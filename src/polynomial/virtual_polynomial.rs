use super::mle::MultilinearExtension;
use crate::field::{ChallengeField, Field};

#[derive(Clone)]
pub struct VirtualPolynomial<F: Field> {
    pub num_vars: usize,
    pub degree: usize,
    pub constituents: Vec<MultilinearExtension<F>>,
    pub evals: Vec<F>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Returns an evaluation on the hypercube.
    #[inline(always)]
    pub fn eval(&self, index: usize) -> F {
        debug_assert!(index < self.evals.len());
        self.evals[index]
    }

    /// Sets the first variable of the polynomial to `point`. Used for fixing sumcheck challenges
    /// into the polynomial.
    pub fn fix_variable(&mut self, point: F) {
        // TODO: par
        // TODO: binius highly optimizes this
        for i in 0..(1 << (self.num_vars - 1)) {
            let mut s = self.evals[(i << 1) + 1];
            s.sub_assign(&self.evals[i << 1]);
            let mut res = point;
            res.mul_assign(&s);
            res.add_assign(&self.evals[i << 1]);
            self.evals[i] = res;
        }

        self.num_vars -= 1;
        self.evals.truncate(1 << self.num_vars);
    }

    pub fn fix_variable_ext<E: ChallengeField<F>>(&self, point: E) -> VirtualPolynomial<E> {
        let mut new_evals = vec![E::ZERO; 1 << (self.num_vars - 1)];
        for i in 0..(1 << (self.num_vars - 1)) {
            let mut s = self.evals[(i << 1) + 1];
            s.sub_assign(&self.evals[i << 1]);
            let mut res = point;
            res.mul_base(&s);
            res.add_base(&self.evals[i << 1]);
            new_evals[i] = res;
        }

        VirtualPolynomial::<E> {
            num_vars: self.num_vars - 1,
            degree: self.degree,
            constituents: vec![],
            evals: new_evals,
        }
    }

    pub fn sum_evaluations(&self) -> Vec<F> {
        // We want to sum up all the evaluations of the polynomial `f` except for those at the
        // variable at the start. We assume lexicographic ordering, so poly[0] and poly[1] will be
        // f(0, 0, ..., 0) and f(1, 0, ..., 0).
        let mut evals_0 = F::ZERO;
        let mut evals_1 = F::ZERO;
        for i in 0..(1 << (self.num_vars() - 1)) {
            evals_0.add_assign(&self.eval(i << 1));
            evals_1.add_assign(&self.eval((i << 1) + 1));
        }

        vec![evals_0, evals_1]
    }
}
