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
pub struct SumcheckProof<F: Field, PCS: PolynomialCommitmentScheme<F>> {
    proofs: Vec<Vec<F>>,
    claimed_sum: F,
    commitment: PCS::Commitment,
    proof: PCS::Proof,
    res: F,
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
    PCS: PolynomialCommitmentScheme<E>,
>(
    poly: &MultilinearExtension<F>,
    transcript: &mut T,
) -> SumcheckProof<E, PCS> {
    let n_rounds = poly.num_vars();
    let mut proofs = Vec::with_capacity(n_rounds);
    let mut challenges = Vec::with_capacity(n_rounds);
    let mut poly_lifted = poly.lift::<E>();

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

    // For the remaining rounds, we always start by fixing the polynomial on a challenge
    // element, and then performing sumcheck steps accordingly.
    for i in 0..(n_rounds - 1) {
        transcript.observe_witnesses(
            &proofs[i]
                .iter()
                .flat_map(|c| Into::<Vec<F>>::into(*c))
                .collect::<Vec<F>>(),
        );
        let challenge = transcript.draw_challenge_ext::<E>();
        challenges.push(challenge);

        poly_lifted.fix_variable(challenge);
        proofs.push(sumcheck_step(&poly_lifted).0);
    }

    // Here we need to commit to the full polynomial and pack it into the proof with an evaluation.
    // This provides the verifier with oracle access to the concerning polynomial and lets her
    // check the final summation in the verification procedure.

    // Retrieve the final sum at which we open the committed poly.
    transcript.observe_witnesses(
        &proofs[n_rounds - 1]
            .iter()
            .flat_map(|c| Into::<Vec<F>>::into(*c))
            .collect::<Vec<F>>(),
    );
    let challenge = transcript.draw_challenge_ext::<E>();
    challenges.push(challenge);
    poly_lifted.fix_variable(challenge);
    debug_assert!(poly_lifted.num_vars() == 0);
    debug_assert!(poly_lifted.evals.len() == 1);
    let res = poly_lifted.evals[0];

    // Do PCS work here and wrap up proof.
    let commitment = PCS::commit(&poly.lift::<E>());
    let proof = PCS::open(&commitment, challenges, res);

    SumcheckProof {
        proofs,
        claimed_sum: claimed_sum.into(),
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
pub fn verify<F: Field, T: Transcript<F>, PCS: PolynomialCommitmentScheme<F>>(
    proof: SumcheckProof<F, PCS>,
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
    for coeffs in proofs.into_iter() {
        transcript.observe_witnesses(&coeffs);
        let c = transcript.draw_challenge();
        challenges.push(c);

        // NOTE: the first step of this verification should take place in the base field
        let mut res = univariate_eval(&coeffs, F::ZERO);
        res.add_assign(&univariate_eval(&coeffs, F::ONE));
        if res != claimed_sum {
            return false;
        }

        claimed_sum = univariate_eval(&coeffs, c);
    }

    // Finally, check the committed polynomial at the list of challenges.
    PCS::verify(&commitment, challenges, res, proof)
}

// Standard lagrange interpolation, assuming indices for evals are 0, 1, 2, ...
// NOTE: can be sped up if we precompute lagrange coeffs for a given size
fn lagrange_interpolation<F: Field>(evals: &[F]) -> Vec<F> {
    let mut polynomial = vec![F::ZERO; evals.len()];
    evals.iter().enumerate().for_each(|(i, eval)| {
        let denom = {
            (0..evals.len()).fold(F::ONE, |mut acc, j| {
                let mut i = F::from_usize(i);
                let j = F::from_usize(j);
                if i != j {
                    i.sub_assign(&j);
                    acc.mul_assign(&i);
                    acc
                } else {
                    acc
                }
            })
        };

        let denom_inv = denom
            .inverse()
            .expect("lagrange coefficient denominator can not be zero");
        let coeffs = {
            let mut coeffs = vec![F::ZERO; evals.len()];
            coeffs[0] = denom_inv;

            for k in 0..evals.len() {
                if k == i {
                    continue;
                }

                let start = if k < i { k + 1 } else { k };
                let mut new_coeffs = vec![F::ZERO; evals.len()];
                for j in (0..start).rev() {
                    new_coeffs[j + 1].add_assign(&coeffs[j]);
                    let mut i = F::from_usize(i);
                    i.mul_assign(&coeffs[j]);
                    new_coeffs[j].sub_assign(&i);
                }
                coeffs = new_coeffs;
            }

            coeffs
        };

        polynomial.iter_mut().zip(coeffs.iter()).for_each(|(v, c)| {
            let mut eval = *eval;
            eval.mul_assign(c);
            v.add_assign(&eval);
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
            acc.mul_assign(&point.pow(i as u32));
            acc
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::m31::{quartic::M31_4, M31};
    use rand::Rng;
    use std::marker::PhantomData;

    pub struct MockPCS<F: Field> {
        _marker: PhantomData<F>,
    }

    impl<F: Field> PolynomialCommitmentScheme<F> for MockPCS<F> {
        type Commitment = usize;
        type Proof = usize;

        fn commit(poly: &MultilinearExtension<F>) -> Self::Commitment {
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
        let mut evals = vec![M31::default(); 2u32.pow(16) as usize];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = MultilinearExtension::new(evals);
        let mut transcript = MockTranscript {
            counter: 1,
            _marker: PhantomData::<M31>,
        };
        let proof = prove::<_, M31_4, _, MockPCS<M31_4>>(&poly, &mut transcript);
    }
}
