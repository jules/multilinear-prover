//! Wrapper for a full logup verifier, combining the IOP and PCS elements.

use crate::{
    arguments::logup::prover::LogUpProof,
    field::{ChallengeField, Field},
    iop::{
        logup,
        sumcheck::{self, SumcheckError},
        zerocheck,
    },
    pcs::PolynomialCommitmentScheme,
    polynomial::{MultilinearExtension, VirtualPolynomial},
    transcript::{IntoObservable, Transcript},
};
use core::marker::PhantomData;
use rayon::prelude::*;

pub struct LogUpVerifier<
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
    > LogUpVerifier<F, E, T, PCS>
{
    pub fn new(transcript: T, pcs: PCS) -> Self {
        Self {
            transcript,
            pcs,
            _f_marker: PhantomData::<F>,
            _e_marker: PhantomData::<E>,
        }
    }

    pub fn verify(&mut self, proof: &LogUpProof<F, E, T, PCS>) -> Result<bool, SumcheckError> {
        // Observe multiplicities.
        self.transcript
            .observe_hashes(&proof.multiplicities_commitment.into_observable());

        let _x = self.transcript.draw_challenge();

        // Observe helpers.
        self.transcript
            .observe_hashes(&proof.helpers_commitment.into_observable());

        let batching_challenges = (0..proof.num_helpers)
            .map(|_| self.transcript.draw_challenge())
            .collect::<Vec<F>>();

        // Draw a list of challenges with which we create the lagrange kernel.
        let mut c = vec![F::ZERO; proof.trace_len];
        c.iter_mut()
            .for_each(|e| *e = self.transcript.draw_challenge());

        // Ensure correct zerocheck.
        if proof.zerocheck_proof.claimed_sum != F::ZERO {
            return Ok(false);
        }

        let (final_claim, eval_point) =
            sumcheck::verify(&proof.zerocheck_proof, &mut self.transcript)?;

        // Ensure correctness of logup constraint.

        // Observe remaining sumcheck constituents commitment.
        self.transcript
            .observe_hashes(&proof.sumcheck_commitment.into_observable());

        // Verify all evaluations.
        Ok(self.pcs.verify(
            &proof.multiplicities_commitment,
            &eval_point,
            &[proof.evaluations[0]],
            &proof.multiplicities_proof,
            &mut self.transcript,
        ) && self.pcs.verify(
            &proof.helpers_commitment,
            &eval_point,
            &proof.evaluations[1..proof.num_helpers + 1],
            &proof.helpers_proof,
            &mut self.transcript,
        ) && self.pcs.verify(
            &proof.sumcheck_commitment,
            &eval_point,
            &proof.evaluations[proof.num_helpers + 1..],
            &proof.sumcheck_proof,
            &mut self.transcript,
        ))
    }
}
