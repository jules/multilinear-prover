//! Sumcheck prover/verifier

use crate::{
    field::{ChallengeField, Field},
    polynomial::{MultilinearExtension, MultivariatePolynomial, VirtualPolynomial},
    transcript::Transcript,
    univariate_utils::*,
};
use rayon::prelude::*;

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
/// their variables. This also outputs the final evaluation point constructed by the prover so that
/// we can re-use it for commitment purposes.
///
/// Includes (naive) use of extension field challenges.
// TODO(opt): https://eprint.iacr.org/2024/108.pdf,
// https://people.cs.georgetown.edu/jthaler/small-sumcheck.pdf,
// https://eprint.iacr.org/2024/1210.pdf
pub fn prove<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    poly: &VirtualPolynomial<F>,
    transcript: &mut T,
    precomputed: &[Vec<F>],
) -> (SumcheckProof<F, E>, Vec<E>) {
    let n_rounds = poly.num_vars();
    let degree = poly.degree();
    let mut proofs = Vec::with_capacity(n_rounds);
    let mut challenges = Vec::with_capacity(n_rounds);

    // In the first round we have no challenge to fix the polynomial with, and we also need to
    // collect the claimed sum from this step; this allows the verifier to reductively check all
    // other claimed sums from just a single field element.
    // So, in this first round, we unroll the logic and do it manually.
    let (final_poly, final_evals) = construct_round_poly(poly, degree, precomputed);

    let mut claimed_sum = final_evals[0].clone();
    claimed_sum.add_assign(&final_evals[1]);
    let coeffs = final_poly
        .into_iter()
        .map(|c| E::from(c))
        .collect::<Vec<E>>();
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

    let precomputed = precomputed
        .iter()
        .map(|v| v.iter().map(|e| Into::<E>::into(*e)).collect::<Vec<E>>())
        .collect::<Vec<Vec<E>>>();

    // For the remaining rounds (except for last), we always start by summing the evaluations,
    // interpolating the intermediate polynomial and then generating and fixing a new challenge.
    for _ in 0..(n_rounds - 1) {
        let (final_poly, _) = construct_round_poly(&poly_lifted, degree, &precomputed);
        transcript.observe_witnesses(
            &final_poly
                .iter()
                .flat_map(|c| Into::<Vec<F>>::into(*c))
                .collect::<Vec<F>>(),
        );
        let challenge = transcript.draw_challenge_ext::<E>();
        proofs.push(final_poly);
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

fn construct_round_poly<F: Field>(
    poly: &VirtualPolynomial<F>,
    degree: usize,
    precomputed: &[Vec<F>],
) -> (Vec<F>, Vec<F>) {
    let individual_evals = poly
        .constituents
        .par_iter()
        .map(|p| sumcheck_step(p, degree))
        .collect::<Vec<Vec<F>>>();
    // XXX virtual polynomial constraint should replace this N degree mult
    let final_evals = (0..individual_evals[0].len())
        .into_par_iter()
        .map(|i| {
            (0..individual_evals.len()).fold(F::ONE, |mut acc, j| {
                acc.mul_assign(&individual_evals[j][i]);
                acc
            })
        })
        .collect::<Vec<F>>();
    let final_poly = lagrange_interpolation_with_precompute(&final_evals, precomputed);

    (final_poly, final_evals)
}

#[inline(always)]
fn sumcheck_step<F: Field>(poly: &MultilinearExtension<F>, degree: usize) -> Vec<F> {
    let evals = poly.sum_evaluations();

    // Extrapolate any extra `degree - 1` coefficients.
    // Taking into account that any multivariate polynomial can be evaluated at a specific point by
    // the equation (letting x be the point we want to check) p(x, X_1, ..., X_n) = (1 - x) * p(0,
    // X_1, ..., X_n) + x * p(1, X_1, ..., X_n).
    let mut extended_evals = Vec::with_capacity(degree + 1);
    extended_evals.push(evals.0);
    extended_evals.push(evals.1);
    for i in 1..degree {
        let mut a = F::from_usize(i + 1);
        let mut one_minus_a = F::ONE;
        one_minus_a.sub_assign(&a);
        a.mul_assign(&evals.1);
        one_minus_a.mul_assign(&evals.0);
        a.add_assign(&one_minus_a);
        extended_evals.push(a);
    }

    extended_evals
}

/// Runs the sumcheck verifier. On being given:
/// - A list of coefficient sets, one per round (or polynomial variable)
/// - An initial claimed sum
/// the verifier can then reductively run the sumcheck protocol and extract a 'final claim',
/// which she can check against (an oracle of) the proven polynomial.
// TODO(opt): the highest coeff can be omitted and should simplify the verification procedure. ref:
// angus gruen's paper https://eprint.iacr.org/2024/108.pdf
pub fn verify<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    proof: &SumcheckProof<F, E>,
    transcript: &mut T,
) -> Result<(E, Vec<E>), SumcheckError> {
    let SumcheckProof {
        proofs,
        claimed_sum,
    } = proof;

    let mut challenge_point = Vec::with_capacity(proofs.len());

    // For each step we:
    // - Draw a challenge based on the polynomial coefficients (just as we do in the prover)
    // - Check the polynomial at 0 and 1, to ensure equality with the claimed sum
    // - Reduce the claimed sum by evaluating the polynomial at the challenge point
    let base_field_coeffs = proofs[0].iter().map(|c| c.real_coeff()).collect::<Vec<F>>();
    let mut res = univariate_eval(&base_field_coeffs, F::ZERO);
    res.add_assign(&univariate_eval(&base_field_coeffs, F::ONE));
    if res != *claimed_sum {
        return Err(SumcheckError::MismatchedSum);
    }

    transcript.observe_witnesses(
        &proofs[0]
            .iter()
            .flat_map(|e| Into::<Vec<F>>::into(*e))
            .collect::<Vec<F>>(),
    );
    let c = transcript.draw_challenge_ext();
    challenge_point.push(c);
    let mut claimed_sum = univariate_eval(&proofs[0], c);

    // We performed the base field check so now we proceed into the extension field.
    for coeffs in proofs.into_iter().skip(1) {
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
        challenge_point.push(c);
        claimed_sum = univariate_eval(&coeffs, c);
    }

    // Output the final claimed sum. The verifier should use this to check against the proven
    // polynomial.
    Ok((claimed_sum, challenge_point))
}
