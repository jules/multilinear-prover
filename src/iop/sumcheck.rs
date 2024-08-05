//! Sumcheck prover/verifier

use crate::{
    field::{ChallengeField, Field},
    mle::MultilinearExtension,
    pcs::PolynomialCommitmentScheme,
    transcript::Transcript,
    univariate_utils::*,
};

/// Possible errors coming from sumcheck verification.
pub enum SumcheckError {
    MismatchedSum,
}

/// A proof produced by running the Sumcheck protocol. Contains the interpolated polynomials at
/// each variable and the initial claimed sum, which allows a verifier to reduce the initial claim
/// to the final one and check it against (an oracle of) the proven polynomial.
pub struct SumcheckProof<F: Field, E: ChallengeField<F>> {
    pub proofs: Vec<Vec<E>>,
    pub claimed_sum: F,
}

/// Runs the sumcheck prover. Given a list of polynomials and some abstracted transcript, we
/// iteratively perform sumcheck reduction steps, in which we evaluate the polynomials in all of
/// their variables. This also outputs the final evaluation of the polynomial itself so the prover
/// can re-use it for commitment purposes.
///
/// Includes (naive) use of extension field challenges.
// TODO(opt): https://eprint.iacr.org/2024/108.pdf,
// https://people.cs.georgetown.edu/jthaler/small-sumcheck.pdf,
// https://eprint.iacr.org/2024/1210.pdf
pub fn prove<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    polys: &[MultilinearExtension<F>],
    transcript: &mut T,
) -> (SumcheckProof<F, E>, Vec<E>) {
    debug_assert!(polys.iter().all(|p| p.evals.len() == polys[0].evals.len()));

    let n_rounds = polys[0].num_vars();
    let degree = polys.len();
    let mut proofs = Vec::with_capacity(n_rounds);
    let mut challenges = Vec::with_capacity(n_rounds);

    // Perform the polynomial product. This is just entrywise mult.
    // XXX i think this is something that the sumcheck prover actually shouldn't do, and instead we
    // should send some abstraction like a 'CompositionPolynomial' which contains degree
    // information.
    let mut poly = polys[0].clone();
    for p in &polys[1..] {
        poly.evals
            .iter_mut()
            .zip(p.evals.iter())
            .for_each(|(eval, m)| eval.mul_assign(m));
    }

    // In the first round we have no challenge to fix the polynomial with, and we also need to
    // collect the claimed sum from this step; this allows the verifier to reductively check all
    // other claimed sums from just a single field element.
    // So, in this first round, we unroll the logic and do it manually.
    let (coeffs, evals) = sumcheck_step(&poly, degree);
    let claimed_sum = evals.iter().fold(F::ZERO, |mut acc, x| {
        acc.add_assign(x);
        acc
    });
    let coeffs = coeffs.into_iter().map(|c| E::from(c)).collect::<Vec<E>>();
    transcript.observe_witnesses(
        &coeffs
            .iter()
            .flat_map(|c| Into::<Vec<F>>::into(*c))
            .collect::<Vec<F>>(),
    );
    let challenge = transcript.draw_challenge_ext::<E>();
    proofs.push(coeffs);
    challenges.push(challenge);

    // Before the next round, we need to lift to the extension by using the first challenge.
    let mut poly_lifted = poly.fix_variable_ext::<E>(challenge);

    // For the remaining rounds (except for last), we always start by summing the evaluations,
    // interpolating the intermediate polynomial and then generating and fixing a new challenge.
    for i in 0..(n_rounds - 1) {
        let (coeffs, _) = sumcheck_step(&poly_lifted, degree);
        transcript.observe_witnesses(
            &coeffs
                .iter()
                .flat_map(|c| Into::<Vec<F>>::into(*c))
                .collect::<Vec<F>>(),
        );
        let challenge = transcript.draw_challenge_ext::<E>();
        proofs.push(coeffs);
        challenges.push(challenge);

        poly_lifted.fix_variable(challenge); // NOTE this isn't necessary on the very last step
    }

    // Output the final proof and finish.
    (
        SumcheckProof {
            proofs,
            claimed_sum,
        },
        challenges,
    )
}

// In a sumcheck step, we sum the evaluations of a polynomial together and create `d` separate sums
// (where `d` is the degree of the multivariate polynomial), which are then interpolated into
// monomial coefficients.
#[inline(always)]
fn sumcheck_step<F: Field>(poly: &MultilinearExtension<F>, degree: usize) -> (Vec<F>, Vec<F>) {
    let evals = poly.sum_evaluations();

    // Extrapolate any extra `degree - 1` coefficients.
    // Taking into account that any multivariate polynomial can be evaluated at a specific point by
    // the equation (letting x be the point we want to check) p(x, X_1, ..., X_n) = (1 - x) * p(0,
    // X_1, ..., X_n) + x * p(1, X_1, ..., X_n).
    let mut extended_evals = evals.clone();
    for i in 1..degree {
        let mut a = F::from_usize(i + 1);
        let mut one_minus_a = F::ONE;
        one_minus_a.sub_assign(&a);
        a.mul_assign(&evals[1]);
        one_minus_a.mul_assign(&evals[0]);
        a.add_assign(&one_minus_a);
        extended_evals.push(a);
    }

    let mut coeffs = lagrange_interpolation(&extended_evals);
    (coeffs, evals)
}

/// Runs the sumcheck verifier. On being given:
/// - A list of coefficient sets, one per round (or polynomial variable)
/// - An initial claimed sum
/// the verifier can then reductively run the sumcheck protocol and extract a 'final claim',
/// which she can check against (an oracle of) the proven polynomial.
// TODO(opt): the highest coeff can be omitted and should simplify the verification procedure. ref:
// angus gruen's paper https://eprint.iacr.org/2024/108.pdf
pub fn verify<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    proof: SumcheckProof<F, E>,
    transcript: &mut T,
) -> Result<E, SumcheckError> {
    let SumcheckProof {
        proofs,
        mut claimed_sum,
    } = proof;

    // For each step we:
    // - Draw a challenge based on the polynomial coefficients (just as we do in the prover)
    // - Check the polynomial at 0 and 1, to ensure equality with the claimed sum
    // - Reduce the claimed sum by evaluating the polynomial at the challenge point
    let base_field_coeffs = proofs[0].iter().map(|c| c.real_coeff()).collect::<Vec<F>>();
    let mut res = univariate_eval(&base_field_coeffs, F::ZERO);
    res.add_assign(&univariate_eval(&base_field_coeffs, F::ONE));
    if res != claimed_sum {
        return Err(SumcheckError::MismatchedSum);
    }

    transcript.observe_witnesses(
        &proofs[0]
            .iter()
            .flat_map(|e| Into::<Vec<F>>::into(*e))
            .collect::<Vec<F>>(),
    );
    let c = transcript.draw_challenge_ext();
    let mut claimed_sum = univariate_eval(&proofs[0], c);

    // We performed the base field check so now we proceed into the extension field.
    for (i, coeffs) in proofs.into_iter().enumerate().skip(1) {
        let mut res = univariate_eval(&coeffs, E::ZERO);
        res.add_assign(&univariate_eval(&coeffs, E::ONE));
        if res != claimed_sum {
            return Err(SumcheckError::MismatchedSum);
        }

        transcript.observe_witnesses(
            &coeffs
                .iter()
                .flat_map(|e| Into::<Vec<F>>::into(*e))
                .collect::<Vec<F>>(),
        );
        let c = transcript.draw_challenge_ext();
        claimed_sum = univariate_eval(&coeffs, c);
    }

    // Output the final claimed sum. The verifier should use this to check against the proven
    // polynomial.
    Ok(claimed_sum)
}
