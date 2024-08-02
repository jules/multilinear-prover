//! Sumcheck prover/verifier

use crate::{
    field::{ChallengeField, Field},
    mle::MultilinearExtension,
    pcs::PolynomialCommitmentScheme,
    transcript::Transcript,
    univariate_utils::*,
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
    pcs: PCS,
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
    let commitment = pcs.commit(&[poly.clone()]);
    let proof = pcs.prove(&commitment, &[poly.clone()], challenges, transcript);

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
    pcs: PCS,
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
    pcs.verify(&commitment, challenges, claimed_sum, &proof, transcript)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::m31::{quartic::M31_4, M31};
    use rand::Rng;
    use std::marker::PhantomData;

    #[derive(Default)]
    pub struct MockPCS<F: Field, T: Transcript<F>, E: ChallengeField<F>> {
        _f_marker: PhantomData<F>,
        _t_marker: PhantomData<T>,
        _e_marker: PhantomData<E>,
    }

    impl<F: Field, T: Transcript<F>, E: ChallengeField<F>> PolynomialCommitmentScheme<F, T, E>
        for MockPCS<F, T, E>
    {
        type Commitment = usize;
        type Proof = usize;

        fn commit(&self, poly: &[MultilinearExtension<F>]) -> Self::Commitment {
            0
        }

        fn prove(
            &self,
            comm: &Self::Commitment,
            polys: &[MultilinearExtension<F>],
            eval: Vec<E>,
            transcript: &mut T,
        ) -> Self::Proof {
            0
        }

        fn verify(
            &self,
            comm: &Self::Commitment,
            eval: Vec<E>,
            result: E,
            proof: &Self::Proof,
            transcript: &mut T,
        ) -> bool {
            true
        }
    }

    #[derive(Default)]
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
        fn draw_bits(&mut self, bits: usize) -> usize {
            0
        }
        fn observe_witness(&mut self, witness: F) {}
        fn observe_witnesses(&mut self, witness: &[F]) {}
    }

    // In this test, we only assert that the arithmetic concerning the sumcheck is correct - we
    // don't check for correct functioning of the transcript or the PCS.
    fn mock_pcs_sumcheck<F: Field>(polys: &[MultilinearExtension<F>]) {
        let mut transcript = MockTranscript {
            counter: 1,
            _marker: PhantomData::<M31>,
        };
        let proof = prove::<_, M31_4, _, _>(
            polys,
            &mut transcript,
            MockPCS::<M31, MockTranscript<M31>, M31_4>::default(),
        );

        let mut transcript = MockTranscript {
            counter: 1,
            _marker: PhantomData::<M31>,
        };
        assert!(verify::<
            _,
            M31_4,
            _,
            MockPCS<M31, MockTranscript<M31>, M31_4>,
        >(
            proof,
            &mut transcript,
            MockPCS::<M31, MockTranscript<M31>, M31_4>::default()
        ));
    }

    #[test]
    fn mock_pcs_single_poly_test() {
        let mut evals = vec![M31::default(); 2u32.pow(20) as usize];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = MultilinearExtension::new(evals);
        mock_pcs_sumcheck(&[poly]);
    }

    #[test]
    fn mock_pcs_2_poly_test() {
        mock_pcs_sumcheck(
            &(0..2)
                .map(|_| {
                    let mut evals = vec![M31::default(); 2u32.pow(20) as usize];
                    evals
                        .iter_mut()
                        .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
                    MultilinearExtension::new(evals)
                })
                .collect::<Vec<MultilinearExtension<M31>>>(),
        );
    }

    #[test]
    fn mock_pcs_16_poly_test() {
        mock_pcs_sumcheck(
            &(0..16)
                .map(|_| {
                    let mut evals = vec![M31::default(); 2u32.pow(20) as usize];
                    evals
                        .iter_mut()
                        .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
                    MultilinearExtension::new(evals)
                })
                .collect::<Vec<MultilinearExtension<M31>>>(),
        );
    }
}
