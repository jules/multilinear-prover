use super::{MultilinearExtension, MultivariatePolynomial};
use crate::field::{ChallengeField, Field};

#[derive(Clone, Debug)]
pub struct VirtualPolynomial<F: Field> {
    degree: usize,
    pub constituents: Vec<MultilinearExtension<F>>,
    pub evals: Vec<F>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn new(evals: Vec<F>, constituents: Vec<MultilinearExtension<F>>) -> Self {
        todo!()
    }
}

impl<F: Field> MultivariatePolynomial<F> for VirtualPolynomial<F> {
    #[inline(always)]
    fn len(&self) -> usize {
        self.evals.len()
    }

    fn degree(&self) -> usize {
        self.degree
    }

    #[inline(always)]
    fn eval(&self, index: usize) -> F {
        self.evals[index]
    }

    #[inline(always)]
    fn eval_mut(&mut self, index: usize) -> &mut F {
        &mut self.evals[index]
    }

    #[inline(always)]
    fn truncate(&mut self, new_len: usize) {
        self.evals.truncate(new_len);
    }

    fn fix_variable_ext<E: ChallengeField<F>>(&self, point: E) -> VirtualPolynomial<E> {
        VirtualPolynomial::<E> {
            degree: self.degree,
            constituents: vec![], // it's likely we don't need this on a lifted polynomial but TBD
            evals: self.fix_variable_ext_internal(point),
        }
    }

    fn add(&self, other: impl MultivariatePolynomial<F>) -> VirtualPolynomial<F> {
        todo!()
    }

    fn mul(&self, other: impl MultivariatePolynomial<F>) -> VirtualPolynomial<F> {
        todo!()
    }
}
