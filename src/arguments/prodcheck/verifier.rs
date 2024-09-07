//! A verifier for the prodcheck argument.

use super::prover::ProdcheckProof;
use crate::{
    arguments::zerocheck::verifier::ZerocheckVerifier,
    field::{ChallengeField, Field},
    iop::sumcheck::SumcheckError,
    pcs::PolynomialCommitmentScheme,
    transcript::{IntoObservable, Transcript},
};

/// Verifier for the prodcheck argument. Builds entirely on the [`ZerocheckVerifier`] and uses its
/// internal PCS and transcript to perform any verification needed outside of verifying the
/// zerocheck polynomial created in the prodcheck prover.
pub struct ProdcheckVerifier<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    zerocheck_verifier: ZerocheckVerifier<F, E, T, PCS>,
}

impl<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, E, T>,
    > ProdcheckVerifier<F, E, T, PCS>
{
    /// Creates a new prodcheck verifier by feeding the given arguments into the creation of a
    /// zerocheck verifier.
    pub fn new(transcript: T, pcs: PCS) -> Self {
        Self {
            zerocheck_verifier: ZerocheckVerifier::new(transcript, pcs),
        }
    }

    /// Verifies a [`ProdcheckProof`]. Essentially, does two things:
    /// - Ensures that the commitment to the product polynomial (v) equals 1 at (0, 1, ..., 1)
    /// - Verifies the zerocheck proof of the zerocheck polynomial (h)
    pub fn verify(&mut self, proof: &ProdcheckProof<F, E, T, PCS>) -> Result<bool, SumcheckError> {
        let mut eval_point = vec![E::ONE; proof.num_vars];
        eval_point[0] = E::ZERO;
        let result = self.zerocheck_verifier.verify(&proof.zerocheck_proof)?;

        // Observe commitment into transcript before verifying.
        self.zerocheck_verifier
            .transcript
            .observe_hashes(&proof.commitment.into_observable());

        Ok(result
            && self.zerocheck_verifier.pcs.verify(
                &proof.commitment,
                &eval_point,
                &[E::ONE],
                &proof.proof,
                &mut self.zerocheck_verifier.transcript,
            ))
    }

    /// We implement a method to relinquish the transcript out of the verifier after it is finished.
    /// This method will consume the prover and leave the transcript to be used for the next
    /// argument.
    pub fn relinquish_transcript(mut self) -> T {
        let transcript = core::mem::take(&mut self.zerocheck_verifier.transcript);
        drop(self);
        transcript
    }
}
