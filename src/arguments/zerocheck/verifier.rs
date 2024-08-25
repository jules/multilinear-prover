//! Wrapper for a full zerocheck verifier, combining the IOP and PCS elements.

use super::prover::ZerocheckProof;
use crate::{
    field::{ChallengeField, Field},
    iop::{sumcheck::SumcheckError, zerocheck},
    pcs::PolynomialCommitmentScheme,
    transcript::Transcript,
    univariate_utils::univariate_eval,
};
use core::marker::PhantomData;

/// A zerocheck verifier, which, on being given a [`ZerocheckProof`], can check correctness of the
/// entire argument (including IOP and PCS elements).
pub struct ZerocheckVerifier<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    pub transcript: T,
    pub pcs: PCS,
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
}

impl<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, E, T>,
    > ZerocheckVerifier<F, E, T, PCS>
{
    /// Creates a new zerocheck verifier for the given transcript and polynomial commitment scheme.
    pub fn new(transcript: T, pcs: PCS) -> Self {
        Self {
            transcript,
            pcs,
            _f_marker: PhantomData::<F>,
            _e_marker: PhantomData::<E>,
        }
    }

    /// Run the zerocheck argument verifier for a given [`ZerocheckProof`].
    pub fn verify(&mut self, proof: &ZerocheckProof<F, E, T, PCS>) -> Result<bool, SumcheckError> {
        // First off, we need to check that the initial claimed sum is indeed zero.
        let base_field_coeffs = proof.zerocheck_proof.proofs[0]
            .iter()
            .map(|c| c.real_coeff())
            .collect::<Vec<F>>();
        let mut res = univariate_eval(&base_field_coeffs, F::ZERO);
        res.add_assign(&univariate_eval(&base_field_coeffs, F::ONE));
        if res != F::ZERO {
            return Ok(false);
        }

        // We run the IOP verifier and extract the final claim and the evaluation of eq.
        let ((final_sum, eval_point), eq) =
            zerocheck::verify(&proof.zerocheck_proof, &mut self.transcript)?;

        // XXX: naively computing the evaluation of eq at challenge point
        let mut eq_lifted = eq.fix_variable_ext(eval_point[0]);
        for v in eval_point.iter().skip(1) {
            eq_lifted.fix_variable(*v);
        }

        let mut evaluations = proof.evaluations.clone();
        evaluations.push(eq_lifted.evals[0]);

        // We ensure that the constraint applied to all evaluations in the proof, multiplied by the
        // evaluation of eq, equals the final claim of the sumcheck verifier. We compute this final
        // value here.
        // XXX DRY
        let final_value = proof
            .products
            .iter()
            .map(|product| {
                let mut a = evaluations[product.0[0]];
                product.0.iter().skip(1).for_each(|i| {
                    a.mul_assign(&evaluations[*i]);
                });

                // TODO exponent

                if product.2 {
                    a.negate();
                }

                a
            })
            .fold(E::ZERO, |mut acc, x| {
                acc.add_assign(&x);
                acc
            });

        // Check that the value holds w.r.t. the sumcheck verifier, and also ensure that all given
        // polynomial evaluations are actually correct w.r.t. the given commitment and evaluation
        // proof.
        Ok(final_sum == final_value
            && self.pcs.verify(
                &proof.commitment,
                &eval_point,
                &proof.evaluations,
                &proof.proof,
                &mut self.transcript,
            ))
    }

    /// We implement a method to relinquish the transcript out of the verifier after it is finished.
    /// This method will consume the prover and leave the transcript to be used for the next
    /// argument.
    pub fn relinquish_transcript(mut self) -> T {
        let transcript = core::mem::take(&mut self.transcript);
        drop(self);
        transcript
    }
}
