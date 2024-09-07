//! Wrapper for a full logup verifier, combining the IOP and PCS elements.

use crate::{
    arguments::zerocheck::prover::ZerocheckProof,
    field::{ChallengeField, Field},
    iop::{
        logup,
        sumcheck::{self, SumcheckError},
        zerocheck,
    },
    pcs::PolynomialCommitmentScheme,
    polynomial::{MultilinearExtension, VirtualPolynomial},
    transcript::Transcript,
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

    pub fn verify(&mut self, proof: &ZerocheckProof<F, E, T, PCS>) -> Result<bool, SumcheckError> {
        // Observe multiplicities.
        //self.transcript.observe_hashes(
        //    &proof
        //        .multipliticies_commitment
        //        .tree
        //        .elements
        //        .clone()
        //        .into_iter()
        //        .flatten()
        //        .collect::<Vec<[u8; 32]>>(),
        //);

        let _x = self.transcript.draw_challenge();

        // Observe helpers.

        //let batching_challenges = (0..x_plus_columns.len())
        //    .map(|_| self.transcript.draw_challenge())
        //    .collect::<Vec<F>>();

        // Draw a list of challenges with which we create the lagrange kernel.
        //let mut c = vec![F::ZERO; table.len()];
        //c.iter_mut()
        //    .for_each(|e| *e = self.transcript.draw_challenge());

        // Ensure correct zerocheck.
        if proof.zerocheck_proof.claimed_sum != F::ZERO {
            return Ok(false);
        }

        let (final_claim, eval_point) =
            sumcheck::verify(&proof.zerocheck_proof, &mut self.transcript)?;

        // Ensure correctness of logup constraint.

        todo!()
    }
}
