pub mod basefold;
pub mod fri;
pub mod tensor_pcs;
pub use tensor_pcs::TensorPCS;

use crate::{
    field::{ChallengeField, Field},
    polynomial::MultilinearExtension,
    transcript::Transcript,
};

pub trait PolynomialCommitmentScheme<F: Field, E: ChallengeField<F>, T: Transcript<F>> {
    type Commitment;
    type Proof;

    fn commit(&self, polys: &[MultilinearExtension<F>], transcript: &mut T) -> Self::Commitment;
    fn prove(
        &self,
        comm: &Self::Commitment,
        polys: &[MultilinearExtension<F>],
        eval: &[E],
        transcript: &mut T,
    ) -> Self::Proof;
    fn verify(
        &self,
        comm: &Self::Commitment,
        eval: &[E],
        results: &[E],
        proof: &Self::Proof,
        transcript: &mut T,
    ) -> bool;
}
