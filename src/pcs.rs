use crate::{field::Field, mle::MultilinearExtension};

pub trait PolynomialCommitmentScheme<F: Field> {
    type Commitment;
    type Proof;

    fn commit(poly: MultilinearExtension<F>) -> Self::Commitment;
    fn open(comm: Self::Commitment, eval: Vec<F>, result: F) -> Self::Proof;
    fn verify(comm: Self::Commitment, eval: Vec<F>, result: F, proof: Self::Proof) -> bool;
}
