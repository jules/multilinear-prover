//! Multilinear FRI in the extension.

use super::*;
use crate::{
    field::{ChallengeField, Field},
    linear_code::LinearCode,
    merkle_tree::MerkleTree,
    polynomial::mle::MultilinearExtension,
    transcript::Transcript,
};
use core::marker::PhantomData;

pub struct FRIPCS<F: Field, T: Transcript<F>, E: ChallengeField<F>, LC: LinearCode<F>> {
    _f_marker: PhantomData<F>,
    _t_marker: PhantomData<T>,
    _e_marker: PhantomData<E>,
    _lc_marker: PhantomData<LC>,
}

impl<F: Field, T: Transcript<F>, E: ChallengeField<F>, LC: LinearCode<F>>
    PolynomialCommitmentScheme<F, T, E> for FRIPCS<F, T, E, LC>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    type Commitment = usize; // TODO
    type Proof = usize; // TODO

    fn commit(&self, polys: &[VirtualPolynomial<F>], transcript: &mut T) -> Self::Commitment {
        //let mut codewords = vec![];

        // Create codeword
        //let mut codeword = LC::encode(&poly);

        todo!()
    }

    fn prove(
        &self,
        comm: &Self::Commitment,
        polys: &[VirtualPolynomial<F>],
        eval: &[E],
        transcript: &mut T,
    ) -> Self::Proof {
        todo!()
    }

    fn verify(
        &self,
        comm: &Self::Commitment,
        eval: &[E],
        results: &[E],
        proof: &Self::Proof,
        transcript: &mut T,
    ) -> bool {
        todo!()
    }
}
