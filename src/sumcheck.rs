//! Sumcheck prover/verifier

use crate::{
    field::{ChallengeField, Field},
    mle::MultilinearExtension,
    pcs::PolynomialCommitmentScheme,
    transcript::Transcript,
};

/// A proof produced by running the Sumcheck protocol. Contains the interpolated polynomials at
/// each variable, the initial claimed sum, a commitment and an opening proof of the polynomial on
/// which the protocol was ran.
pub struct SumcheckProof<
    F: Field,
    T: Transcript<F>,
    E: ChallengeField<F>,
    PCS: PolynomialCommitmentScheme<F, T, E>,
> {
    proofs: Vec<Vec<E>>,
    claimed_sum: F,
    commitment: PCS::Commitment,
    proof: PCS::Proof,
    res: E,
}

/// Runs the sumcheck prover. Given a polynomial and some abstracted transcript, we:
/// - Iteratively perform sumcheck reduction steps, in which we evaluate the polynomial in all of
/// its variables
/// - Commit to the polynomial in full and produce an opening proof at the set of generated
/// challenges
///
/// Two different fields may be specified in case we want to use challenges from a field extension.
///
/// For the generics: F denotes the base field, and E denotes the extension field.
pub fn prove<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, T, E>,
>(
    poly: &MultilinearExtension<F>,
    transcript: &mut T,
) -> SumcheckProof<F, T, E, PCS> {
    let n_rounds = poly.num_vars();
    let mut proofs = Vec::with_capacity(n_rounds);
    let mut challenges = Vec::with_capacity(n_rounds);

    // In the first round we have no challenge to fix the polynomial with, and we also need to
    // collect the claimed sum from this step. This allows the verifier to reductively check all
    // other claimed sums from just a single field element.
    let (coeffs, evals) = sumcheck_step(&poly);
    let claimed_sum = evals.iter().fold(F::ZERO, |mut acc, x| {
        acc.add_assign(x);
        acc
    });
    let coeffs = coeffs.into_iter().map(|c| E::from(c)).collect::<Vec<E>>();
    proofs.push(coeffs);

    // For the next round, we lift to the extension by using the first challenge.
    transcript.observe_witnesses(
        &proofs[0]
            .iter()
            .flat_map(|c| Into::<Vec<F>>::into(*c))
            .collect::<Vec<F>>(),
    );
    let challenge = transcript.draw_challenge_ext::<E>();
    challenges.push(challenge);

    let mut poly_lifted = poly.fix_variable_ext::<E>(challenge);

    // For the remaining rounds, we always start by fixing the polynomial on a challenge
    // element, and then performing sumcheck steps accordingly.
    for i in 0..(n_rounds - 1) {
        let (coeffs, _) = sumcheck_step(&poly_lifted);
        proofs.push(coeffs);

        transcript.observe_witnesses(
            &proofs[i + 1]
                .iter()
                .flat_map(|c| Into::<Vec<F>>::into(*c))
                .collect::<Vec<F>>(),
        );
        let challenge = transcript.draw_challenge_ext::<E>();
        challenges.push(challenge);
        poly_lifted.fix_variable(challenge);
    }

    // Here we need to commit to the full polynomial and pack it into the proof with an evaluation.
    // This provides the verifier with oracle access to the concerning polynomial and lets her
    // check the final summation in the verification procedure.

    debug_assert!(poly_lifted.num_vars() == 0);
    debug_assert!(poly_lifted.evals.len() == 1);
    // Retrieve the final sum at which we open the committed poly.
    let res = poly_lifted.evals[0];

    // Do PCS work now and wrap up proof.
    let commitment = PCS::commit(&[poly.clone()]);
    let proof = PCS::prove(&commitment, challenges, res, transcript);

    SumcheckProof {
        proofs,
        claimed_sum,
        commitment,
        proof,
        res,
    }
}

// In a sumcheck step, we sum the evaluations of a polynomial together and create `d` separate sums
// (where `d` is the degree of the multivariate polynomial), which are then interpolated into
// monomial coefficients.
#[inline(always)]
fn sumcheck_step<F: Field>(poly: &MultilinearExtension<F>) -> (Vec<F>, Vec<F>) {
    let evals = poly.sum_evaluations();
    let coeffs = lagrange_interpolation(&evals);
    (coeffs, evals)
}

/// Runs the sumcheck verifier. On being given:
/// - A list of coefficient sets, one per round (or polynomial variable)
/// - An initial claimed sum
/// - Oracle access to the polynomial
/// - An evaluation proof of the polynomial oracle
/// the verifier can then successfully run the sumcheck protocol and ensure that the proof is
/// correct.
// TODO(opt): the highest coeff can be omitted and should simplify the verification procedure. ref:
// binius
pub fn verify<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, T, E>,
>(
    proof: SumcheckProof<F, T, E, PCS>,
    transcript: &mut T,
) -> bool {
    let SumcheckProof {
        proofs,
        mut claimed_sum,
        commitment,
        proof,
        res,
    } = proof;

    let mut challenges = Vec::with_capacity(proofs.len());
    // For each step we:
    // - Draw a challenge based on the polynomial coefficients (just as we do in the prover)
    // - Check the polynomial at 0 and 1, to ensure equality with the claimed sum
    // - Reduce the claimed sum by evaluating the polynomial at the challenge point
    let base_field_coeffs = proofs[0].iter().map(|c| c.real_coeff()).collect::<Vec<F>>();
    let mut res = univariate_eval(&base_field_coeffs, F::ZERO);
    res.add_assign(&univariate_eval(&base_field_coeffs, F::ONE));
    if res != claimed_sum {
        return false;
    }

    transcript.observe_witnesses(
        &proofs[0]
            .iter()
            .flat_map(|e| Into::<Vec<F>>::into(*e))
            .collect::<Vec<F>>(),
    );
    let c = transcript.draw_challenge_ext();
    challenges.push(c);
    let mut claimed_sum = univariate_eval(&proofs[0], c);

    // We performed the base field check so now we proceed into the extension field.
    for (i, coeffs) in proofs.into_iter().enumerate().skip(1) {
        let mut res = univariate_eval(&coeffs, E::ZERO);
        res.add_assign(&univariate_eval(&coeffs, E::ONE));
        if res != claimed_sum {
            return false;
        }

        transcript.observe_witnesses(
            &coeffs
                .iter()
                .flat_map(|e| Into::<Vec<F>>::into(*e))
                .collect::<Vec<F>>(),
        );
        let c = transcript.draw_challenge_ext();
        challenges.push(c);
        claimed_sum = univariate_eval(&coeffs, c);
    }

    // Finally, check the committed polynomial at the list of challenges.
    PCS::verify(&commitment, challenges, claimed_sum, proof, transcript)
}

// Standard lagrange interpolation, assuming indices for evals are 0, 1, 2, ...
// NOTE: can be sped up if we precompute lagrange coeffs for a given size
fn lagrange_interpolation<F: Field>(evals: &[F]) -> Vec<F> {
    let multiply_polys = |a: &[F], b: &[F]| -> Vec<F> {
        let mut result = vec![F::ZERO; a.len() + b.len() - 1];
        a.iter().enumerate().for_each(|(i, c1)| {
            b.iter().enumerate().for_each(|(j, c2)| {
                let mut m = c1.clone();
                m.mul_assign(c2);
                result[i + j].add_assign(&m);
            });
        });

        result
    };

    let mut polynomial = vec![F::ZERO; evals.len()];
    evals.iter().enumerate().for_each(|(i, eval)| {
        let mut coeffs = vec![F::ONE];
        evals.iter().enumerate().for_each(|(j, _)| {
            if i != j {
                let denom_inv = {
                    let mut i = F::from_usize(i);
                    let j = F::from_usize(j);
                    i.sub_assign(&j);
                    i.inverse().unwrap()
                };

                let mut new_coeffs = vec![F::ZERO; 2];
                new_coeffs[0] = F::from_usize(j);
                new_coeffs[0].negate();
                new_coeffs[0].mul_assign(&denom_inv);
                new_coeffs[1] = denom_inv;
                coeffs = multiply_polys(&coeffs, &new_coeffs);
            }
        });

        coeffs.iter().enumerate().for_each(|(k, c)| {
            let mut res = eval.clone();
            res.mul_assign(c);
            polynomial[k].add_assign(&res);
        });
    });

    polynomial
}

// Simple Horner's method evaluation.
fn univariate_eval<F: Field>(coeffs: &[F], point: F) -> F {
    coeffs
        .iter()
        .enumerate()
        .rev()
        .fold(F::ZERO, |mut acc, (i, coeff)| {
            acc.add_assign(coeff);
            if i != 0 {
                acc.mul_assign(&point);
            }
            acc
        })
}

// Horner's method evaluation that raises into an extension field.
fn univariate_eval_ext<F: Field, E: ChallengeField<F>>(coeffs: &[F], point: E) -> E {
    let mut init = point.clone();
    init.mul_base(&coeffs[coeffs.len() - 1]);

    coeffs
        .iter()
        .enumerate()
        .rev()
        .skip(1)
        .fold(init, |mut acc, (i, coeff)| {
            acc.add_base(coeff);
            if i != 0 {
                acc.mul_assign(&point);
            }
            acc
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::m31::{quartic::M31_4, M31};
    use rand::Rng;
    use std::marker::PhantomData;

    #[test]
    fn test_lagrange_2() {
        let mut evals = vec![M31::default(); 2];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = lagrange_interpolation(&evals);
        assert_eq!(evals[0], univariate_eval(&poly, M31::ZERO));
        assert_eq!(evals[1], univariate_eval(&poly, M31::ONE));
    }

    #[test]
    fn test_lagrange_3() {
        let mut evals = vec![M31::default(); 3];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = lagrange_interpolation(&evals);
        assert_eq!(evals[0], univariate_eval(&poly, M31::ZERO));
        assert_eq!(evals[1], univariate_eval(&poly, M31::ONE));
        assert_eq!(evals[2], univariate_eval(&poly, M31(2)));
    }

    #[test]
    fn test_lagrange_20() {
        let mut evals = vec![M31::default(); 20];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = lagrange_interpolation(&evals);
        for i in 0..20 {
            assert_eq!(evals[i], univariate_eval(&poly, M31::from_usize(i)));
        }
    }

    #[test]
    fn test_eval_2() {
        // f(x) = 5 + 2x
        let poly = vec![M31(5), M31(2)];

        // f(2) = 9
        assert_eq!(M31(9), univariate_eval(&poly, M31(2)));
    }

    #[test]
    fn test_eval_3() {
        // f(x) = 5 + 2x + 3x^2
        let poly = vec![M31(5), M31(2), M31(3)];

        // f(2) = 21
        assert_eq!(M31(21), univariate_eval(&poly, M31(2)));
    }

    #[test]
    fn test_eval_ext() {
        // f(x) = 5 + 2x
        let poly = vec![M31(5), M31(2)];

        // f([2, 0]) = 9
        assert_eq!(
            M31_4::from_usize(9),
            univariate_eval_ext(&poly, M31_4::from_usize(2))
        );
    }

    pub struct MockPCS<F: Field> {
        _marker: PhantomData<F>,
    }

    impl<F: Field> PolynomialCommitmentScheme<F> for MockPCS<F> {
        type Commitment = usize;
        type Proof = usize;

        fn commit(poly: &[MultilinearExtension<F>]) -> Self::Commitment {
            0
        }

        fn open(comm: &Self::Commitment, eval: Vec<F>, result: F) -> Self::Proof {
            0
        }

        fn verify(comm: &Self::Commitment, eval: Vec<F>, result: F, proof: Self::Proof) -> bool {
            true
        }
    }

    pub struct MockTranscript<F: Field> {
        counter: usize,
        _marker: PhantomData<F>,
    }

    impl<F: Field> Transcript<F> for MockTranscript<F> {
        fn draw_challenge(&mut self) -> F {
            self.counter += 1;
            F::from_usize(self.counter)
        }
        fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E {
            self.counter += 1;
            E::from(F::from_usize(self.counter))
        }
        fn observe_witness(&mut self, witness: F) {}
        fn observe_witnesses(&mut self, witness: &[F]) {}
    }

    #[test]
    fn mock_pcs_test() {
        let mut evals = vec![M31::default(); 2u32.pow(20) as usize];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = MultilinearExtension::new(evals);
        let mut transcript = MockTranscript {
            counter: 1,
            _marker: PhantomData::<M31>,
        };
        let proof = prove::<_, M31_4, _, MockPCS<M31_4>>(&poly, &mut transcript);

        let mut transcript = MockTranscript {
            counter: 1,
            _marker: PhantomData::<M31>,
        };
        assert!(verify::<_, M31_4, _, MockPCS<M31_4>>(
            proof,
            &mut transcript
        ));
    }
}
