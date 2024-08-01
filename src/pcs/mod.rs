pub mod tensor_pcs;

use crate::{
    field::{ChallengeField, Field},
    mle::MultilinearExtension,
    transcript::Transcript,
};

pub trait PolynomialCommitmentScheme<F: Field, T: Transcript<F>, E: ChallengeField<F>> {
    type Commitment;
    type Proof;

    fn commit(polys: &[MultilinearExtension<F>]) -> Self::Commitment;
    fn prove(
        comm: &Self::Commitment,
        polys: &[MultilinearExtension<F>],
        eval: Vec<E>,
        result: E,
        transcript: &mut T,
    ) -> Self::Proof;
    fn verify(
        comm: &Self::Commitment,
        eval: Vec<E>,
        result: E,
        proof: Self::Proof,
        transcript: &mut T,
    ) -> bool;
}
