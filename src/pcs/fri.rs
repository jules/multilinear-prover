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

pub struct FRI<F: Field, E: ChallengeField<F>, T: Transcript<F>, LC: LinearCode<F>> {
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
    _t_marker: PhantomData<T>,
    _lc_marker: PhantomData<LC>,
}

impl<F: Field, E: ChallengeField<F>, T: Transcript<F>, LC: LinearCode<F>>
    PolynomialCommitmentScheme<F, E, T> for FRI<F, E, T, LC>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    type Commitment = usize; // TODO
    type Proof = usize; // TODO

    fn commit(&self, polys: &[MultilinearExtension<F>], transcript: &mut T) -> Self::Commitment {
        //let mut codewords = vec![];

        // Create codeword
        //let mut codeword = LC::encode(&poly);

        todo!()
    }

    fn prove(
        &self,
        comm: &Self::Commitment,
        polys: &[MultilinearExtension<F>],
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
