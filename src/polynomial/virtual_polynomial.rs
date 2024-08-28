use super::MultilinearExtension;
use crate::field::{ChallengeField, Field};
use core::cmp::max;
use rayon::prelude::*;

/// A representation of multiple polynomials mixed together with addition or multiplication.
#[derive(Clone, Debug)]
pub struct VirtualPolynomial<F: Field> {
    degree: usize,
    pub constituents: Vec<MultilinearExtension<F>>,
    // Column indices, negation
    pub products: Vec<(Vec<usize>, bool)>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn new(
        constituents: Vec<MultilinearExtension<F>>,
        degree: usize,
        products: Vec<(Vec<usize>, bool)>,
    ) -> Self {
        Self {
            degree,
            constituents,
            products,
        }
    }

    pub fn add_assign(&mut self, other: &Self) {
        let original_len = self.constituents.len();
        self.constituents.extend(other.constituents.clone());
        self.degree = max(self.degree, other.degree);

        let mut shifted_products = other.products.clone();
        shifted_products.iter_mut().for_each(|product| {
            product.0.iter_mut().for_each(|index| {
                *index = *index + original_len;
            });
        });
        self.products.extend(shifted_products);
    }

    pub fn add_assign_mle(&mut self, other: &MultilinearExtension<F>, negate: bool) {
        self.constituents.push(other.clone());
        self.products
            .push((vec![self.constituents.len() - 1], negate));
    }

    pub fn mul_assign_mle(&mut self, other: &MultilinearExtension<F>) {
        self.constituents.push(other.clone());
        self.degree += 1;
        self.products.iter_mut().for_each(|product| {
            product.0.push(self.constituents.len() - 1);
        });
    }

    pub fn negate_product(&mut self, product_index: usize) {
        self.products[product_index].1 = !self.products[product_index].1;
    }

    #[inline(always)]
    pub fn num_vars(&self) -> usize {
        self.constituents[0].num_vars()
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.constituents[0].len()
    }

    #[inline(always)]
    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn fix_variable(&mut self, point: F) {
        self.constituents
            .par_iter_mut()
            .for_each(|mle| mle.fix_variable(point));
    }

    pub fn fix_variable_ext<E: ChallengeField<F>>(&self, point: E) -> VirtualPolynomial<E> {
        VirtualPolynomial::<E> {
            degree: self.degree,
            constituents: self
                .constituents
                .par_iter()
                .map(|mle| mle.fix_variable_ext(point))
                .collect::<Vec<MultilinearExtension<E>>>(),
            products: self.products.clone(),
        }
    }
}

impl<F: Field> From<MultilinearExtension<F>> for VirtualPolynomial<F> {
    fn from(value: MultilinearExtension<F>) -> Self {
        Self::new(vec![value], 1, vec![(vec![0], false)])
    }
}
