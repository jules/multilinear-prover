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

    pub fn merge(&self, other: Self) -> Self {
        debug_assert!(self.len() == other.len());
        /*
        let mut new_evals = Vec::with_capacity(self.len() * 2);
        for i in 0..self.len() {
            new_evals.push(self.evals[i]);
            new_evals.push(other.evals[i]);
        }
        Self::new(new_evals)
        */
        let mut evals = self.evals.clone();
        evals.extend(other.evals.clone());
        Self::new(evals)
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.evals.len()
    }

    #[inline(always)]
    pub fn num_vars(&self) -> usize {
        self.len().ilog2() as usize
    }

    #[inline(always)]
    pub fn degree(&self) -> usize {
        1
    }

    /// Sets the first variable of the polynomial to `point`. Used for fixing sumcheck challenges
    /// into the polynomial.
    pub fn fix_variable(&mut self, point: F) {
        for i in 0..(1 << (self.num_vars() - 1)) {
            let mut s = self.evals[(i << 1) + 1];
            s.sub_assign(&self.evals[i << 1]);
            let mut res = point;
            res.mul_assign(&s);
            res.add_assign(&self.evals[i << 1]);
            self.evals[i] = res;
        }

        self.evals.truncate(1 << (self.num_vars() - 1));
    }

    pub fn fix_variable_ext<E: ChallengeField<F>>(&self, point: E) -> MultilinearExtension<E> {
        let mut new_evals = vec![E::ZERO; 1 << (self.num_vars() - 1)];
        for i in 0..(1 << (self.num_vars() - 1)) {
            let mut s = self.evals[(i << 1) + 1];
            s.sub_assign(&self.evals[i << 1]);
            let mut res = point;
            res.mul_base(&s);
            res.add_base(&self.evals[i << 1]);
            new_evals[i] = res;
        }

        MultilinearExtension::<E>::new(new_evals)
    }

    /// Sums all the evaluations of the polynomial at f(0, ...) and f(1, ...).
    pub fn sum_evaluations(&self) -> (F, F) {
        // We want to sum up all the evaluations of the polynomial `f` except for those at the
        // variable at the start. We assume lexicographic ordering, so poly[0] and poly[1] will be
        // f(0, 0, ..., 0) and f(1, 0, ..., 0).
        let mut evals_0 = F::ZERO;
        let mut evals_1 = F::ZERO;
        for i in 0..(1 << (self.num_vars() - 1)) {
            evals_0.add_assign(&self.evals[i << 1]);
            evals_1.add_assign(&self.evals[(i << 1) + 1]);
        }

        (evals_0, evals_1)
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
