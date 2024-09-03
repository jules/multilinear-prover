//! Wrapper for a full logup verifier, combining the IOP and PCS elements.

use crate::{
    arguments::zerocheck::prover::ZerocheckProof,
    field::{ChallengeField, Field},
    iop::{logup, sumcheck, zerocheck},
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
        // We will have the evaluations laid out as follows.

        let x = self.transcript.draw_challenge();

        // Draw a list of challenges with which we create the lagrange kernel.
        let mut c = vec![F::ZERO; table.len()];
        c.iter_mut()
            .for_each(|e| *e = self.transcript.draw_challenge());

        todo!()
    }
}
