use super::{MultilinearExtension, MultivariatePolynomial};
use crate::field::{ChallengeField, Field};
use core::cmp::max;
use rayon::prelude::*;

/// A representation of multiple polynomials mixed together with addition or multiplication.
// NOTE no bookkeeping on how mixings were performed
#[derive(Clone, Debug)]
pub struct VirtualPolynomial<F: Field> {
    degree: usize,
    pub constituents: Vec<MultilinearExtension<F>>,
    pub evals: Vec<F>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn new(evals: Vec<F>, constituents: Vec<MultilinearExtension<F>>, degree: usize) -> Self {
        Self {
            degree,
            constituents,
            evals,
        }
    }

    pub fn add_assign(&mut self, other: &Self) {
        self.constituents.extend(other.constituents.clone());
        self.evals
            .iter_mut()
            .zip(other.evals.iter())
            .for_each(|(e, o)| {
                e.add_assign(o);
            });
        self.degree = max(self.degree, other.degree);
    }

    pub fn mul_assign(&mut self, other: &Self) {
        self.constituents.extend(other.constituents.clone());
        self.evals
            .iter_mut()
            .zip(other.evals.iter())
            .for_each(|(e, o)| {
                e.mul_assign(o);
            });
        self.degree += other.degree;
    }

    pub fn add_assign_mle(&mut self, other: &MultilinearExtension<F>) {
        self.constituents.push(other.clone());
        self.evals
            .iter_mut()
            .zip(other.evals.iter())
            .for_each(|(e, o)| {
                e.add_assign(o);
            });
    }

    pub fn mul_assign_mle(&mut self, other: &MultilinearExtension<F>) {
        self.constituents.push(other.clone());
        self.evals
            .iter_mut()
            .zip(other.evals.iter())
            .for_each(|(e, o)| {
                e.mul_assign(o);
            });
        self.degree += 1;
    }
}

impl<F: Field> MultivariatePolynomial<F> for VirtualPolynomial<F> {
    #[inline(always)]
    fn len(&self) -> usize {
        self.evals.len()
    }

    #[inline(always)]
    fn degree(&self) -> usize {
        self.degree
    }

    #[inline(always)]
    fn evals(&self) -> &[F] {
        &self.evals
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

    fn fix_variable(&mut self, point: F) {
        self.constituents
            .par_iter_mut()
            .for_each(|mle| mle.fix_variable(point));
    }

    #[inline(always)]
    fn fix_variable_ext<E: ChallengeField<F>>(&self, point: E) -> VirtualPolynomial<E> {
        VirtualPolynomial::<E> {
            degree: self.degree,
            constituents: self
                .constituents
                .par_iter()
                .map(|mle| mle.fix_variable_ext(point))
                .collect::<Vec<MultilinearExtension<E>>>(),
            evals: vec![],
        }
    }
}

impl<F: Field> From<MultilinearExtension<F>> for VirtualPolynomial<F> {
    fn from(value: MultilinearExtension<F>) -> Self {
        Self::new(value.evals.clone(), vec![value], 1)
    }
}
