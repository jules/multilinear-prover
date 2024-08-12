pub mod fri;
pub mod tensor_pcs;

use crate::{
    field::{ChallengeField, Field},
    polynomial::VirtualPolynomial,
    transcript::Transcript,
};

pub trait PolynomialCommitmentScheme<F: Field, T: Transcript<F>, E: ChallengeField<F>> {
    type Commitment;
    type Proof;

    fn commit(&self, polys: &[VirtualPolynomial<F>], transcript: &mut T) -> Self::Commitment;
    fn prove(
        &self,
        comm: &Self::Commitment,
        polys: &[VirtualPolynomial<F>],
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
