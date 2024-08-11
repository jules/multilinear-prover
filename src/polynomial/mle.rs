use super::{MultivariatePolynomial, VirtualPolynomial};
use crate::field::{ChallengeField, Field};

/// A multilinear polynomial represented in the Lagrange basis.
#[derive(Clone, Debug)]
pub struct MultilinearExtension<F: Field> {
    pub evals: Vec<F>,
}

impl<F: Field> MultilinearExtension<F> {
    pub fn new(evals: Vec<F>) -> Self {
        debug_assert!(evals.len().is_power_of_two());
        Self { evals }
    }
}

impl<F: Field> MultivariatePolynomial<F> for MultilinearExtension<F> {
    #[inline(always)]
    fn len(&self) -> usize {
        self.evals.len()
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

    fn fix_variable_ext<E: ChallengeField<F>>(&self, point: E) -> MultilinearExtension<E> {
        MultilinearExtension::<E>::new(self.fix_variable_ext_internal(point))
    }

    fn add(&self, other: impl MultivariatePolynomial<F>) -> VirtualPolynomial<F> {
        todo!()
    }

    fn mul(&self, other: impl MultivariatePolynomial<F>) -> VirtualPolynomial<F> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{m31::M31, Field};
    use rand::Rng;

    #[test]
    fn test_fixing_consistency() {
        let mut evals = vec![M31::default(); 16];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..5)));
        let mut poly = MultilinearExtension::new(evals);

        let c = M31(2);

        let summed_fixed = {
            let evals = poly.sum_evaluations();
            let mut p0 = evals.0;
            let p1 = evals.1;

            // p0 + c * (p1 - p0)
            let mut rhs = p1;
            rhs.sub_assign(&p0);
            let mut c = c.clone();
            c.mul_assign(&rhs);

            p0.add_assign(&c);
            p0
        };

        let fn_fixed = {
            poly.fix_variable(c);
            let evals = poly.sum_evaluations();
            let mut p0 = evals.0;
            let p1 = evals.1;
            p0.add_assign(&p1);
            p0
        };

        assert_eq!(summed_fixed, fn_fixed);
    }
}
