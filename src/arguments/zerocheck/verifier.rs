//! Wrapper for a full zerocheck verifier, combining the IOP and PCS elements.

use super::prover::ZerocheckProof;
use crate::{
    field::{ChallengeField, Field},
    iop::{sumcheck::SumcheckError, zerocheck},
    pcs::PolynomialCommitmentScheme,
    polynomial::MultivariatePolynomial,
    transcript::Transcript,
};
use core::marker::PhantomData;

pub struct ZerocheckVerifier<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    transcript: T,
    pcs: PCS,
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
        // First, we run the IOP verifier and extract the final claim and the evaluation of eq.
        // TODO: eq eval and not just recreating it entirely.
        let ((final_sum, eval_point), eq) =
            zerocheck::verify(&proof.zerocheck_proof, &mut self.transcript)?;

        // XXX: naively computing the evaluation of eq at challenge point
        let mut eq_lifted = eq.fix_variable_ext(eval_point[0]);
        for v in eval_point.iter().skip(1) {
            eq_lifted.fix_variable(*v);
        }

        // We ensure that the constraint applied to all evaluations in the proof, multiplied by the
        // evaluation of eq, equals the final claim of the sumcheck verifier. Otherwise, reject.
        // XXX virtual polynomial tooling should replace this N degree multiplication.
        let mut final_value = proof.evaluations.iter().fold(E::ONE, |mut acc, x| {
            acc.mul_assign(x);
            acc
        });
        final_value.mul_assign(&eq_lifted.evals[0]);

        // Check that the value holds w.r.t. the sumcheck verifier, and also ensure that all given
        // polynomial evaluations are actually correct w.r.t. the given commitment and evaluation
        // proof.
        Ok(final_sum != final_value
            && self.pcs.verify(
                &proof.commitment,
                &eval_point,
                &proof.evaluations,
                &proof.proof,
                &mut self.transcript,
            ))
    }

    /// We implement a method to relinquish the transcript out of the prover after it is finished.
    /// This method will consume the prover and leave the transcript to be used for the next
    /// argument.
    pub fn relinquish_transcript(mut self) -> T {
        let transcript = core::mem::take(&mut self.transcript);
        drop(self);
        transcript
    }
}
