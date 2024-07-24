use crate::field::{ChallengeField, Field};

/// A multilinear polynomial represented in the Lagrange basis.
#[derive(Clone, Debug)]
pub struct MultilinearExtension<F: Field> {
    num_vars: usize,
    pub evals: Vec<F>,
}

impl<F: Field> MultilinearExtension<F> {
    pub fn new(evals: Vec<F>) -> Self {
        debug_assert!(evals.len().is_power_of_two());
        Self {
            num_vars: (evals.len() as f64).log2() as usize,
            evals,
        }
    }

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns an evaluation on the hypercube.
    #[inline(always)]
    pub fn eval(&self, index: usize) -> F {
        debug_assert!(index < self.evals.len());
        self.evals[index]
    }

    /// Lifts the MLE to an extension field. Used for accommodating extension field challenge
    /// evaluations and such.
    pub fn lift<E: ChallengeField<F>>(&self) -> MultilinearExtension<E> {
        let new_evals = self.evals.iter().map(|e| (*e).into()).collect::<Vec<E>>();
        MultilinearExtension::<E>::new(new_evals)
    }

    /// Sets the first variable of the polynomial to `point`. Used for fixing sumcheck challenges
    /// into the polynomial. Supports the use of an extension field challenge.
    pub fn fix_variable(&mut self, point: F) {
        // TODO: par
        // TODO: binius highly optimizes this
        for i in 0..(1 << (self.num_vars - 1)) {
            let mut s = self.evals[(i << 1) + 1];
            s.sub_assign(&self.evals[i << 1]);
            let mut res = point;
            res.mul_assign(&s);
            self.evals[i].add_assign(&res);
        }

        self.num_vars -= 1;
        self.evals.truncate(1 << self.num_vars);
    }

    /// Sums all the evaluations of the polynomial at f(0, ...) and f(1, ...).
    pub fn sum_evaluations(&self) -> Vec<F> {
        // We want to sum up all the evaluations of the polynomial `f` except for those at the
        // variable at the start. We assume lexicographic ordering, so poly[0] and poly[1] will be
        // f(0, 0, ..., 0) and f(1, 0, ..., 0).
        let mut evals_0 = F::ZERO;
        let mut evals_1 = F::ZERO;
        for i in 0..self.num_vars() - 1 {
            evals_0.add_assign(&self.eval(i << 1));
            evals_1.add_assign(&self.eval((i << 1) + 1));
        }

        vec![evals_0, evals_1]
    }
}

pub struct MultilinearComposite<F: Field> {
    polys: Vec<MultilinearExtension<F>>,
    num_vars: usize,
    evals: Vec<F>,
}
