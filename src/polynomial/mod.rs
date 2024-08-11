pub mod mle;
pub use mle::*;
pub mod virtual_polynomial;
pub use virtual_polynomial::*;

use crate::field::{ChallengeField, Field};
use core::fmt::Debug;

/// Mainly takes care of the DRY principle for common functionality in polynomial types.
pub trait MultivariatePolynomial<F: Field>: Clone + Debug {
    fn len(&self) -> usize;

    #[inline(always)]
    fn num_vars(&self) -> usize {
        debug_assert!(self.len().is_power_of_two());
        self.len().ilog2() as usize
    }

    #[inline(always)]
    fn degree(&self) -> usize {
        1
    }

    /// Returns a reference to the underlying evaluations.
    fn evals(&self) -> &[F];

    /// Returns an evaluation on the hypercube.
    fn eval(&self, index: usize) -> F;

    /// Returns a mutable reference to an evaluation on the hypercube.
    fn eval_mut(&mut self, index: usize) -> &mut F;

    /// Truncates the evaluation vector, useful for bookkeeping during sumcheck.
    fn truncate(&mut self, new_len: usize);

    /// Sets the first variable of the polynomial to `point`. Used for fixing sumcheck challenges
    /// into the polynomial.
    fn fix_variable(&mut self, point: F) {
        // TODO: par
        // TODO: binius highly optimizes this
        for i in 0..(1 << (self.num_vars() - 1)) {
            let mut s = self.eval((i << 1) + 1);
            s.sub_assign(&self.eval(i << 1));
            let mut res = point;
            res.mul_assign(&s);
            res.add_assign(&self.eval(i << 1));
            *self.eval_mut(i) = res;
        }

        self.truncate(1 << (self.num_vars() - 1));
    }

    fn fix_variable_ext<E: ChallengeField<F>>(&self, point: E) -> impl MultivariatePolynomial<E>;
    fn fix_variable_ext_internal<E: ChallengeField<F>>(&self, point: E) -> Vec<E> {
        let mut new_evals = vec![E::ZERO; 1 << (self.num_vars() - 1)];
        for i in 0..(1 << (self.num_vars() - 1)) {
            let mut s = self.eval((i << 1) + 1);
            s.sub_assign(&self.eval(i << 1));
            let mut res = point;
            res.mul_base(&s);
            res.add_base(&self.eval(i << 1));
            new_evals[i] = res;
        }

        new_evals
    }

    /// Sums all the evaluations of the polynomial at f(0, ...) and f(1, ...).
    fn sum_evaluations(&self) -> (F, F) {
        // We want to sum up all the evaluations of the polynomial `f` except for those at the
        // variable at the start. We assume lexicographic ordering, so poly[0] and poly[1] will be
        // f(0, 0, ..., 0) and f(1, 0, ..., 0).
        let mut evals_0 = F::ZERO;
        let mut evals_1 = F::ZERO;
        for i in 0..(1 << (self.num_vars() - 1)) {
            evals_0.add_assign(&self.eval(i << 1));
            evals_1.add_assign(&self.eval((i << 1) + 1));
        }

        (evals_0, evals_1)
    }
}
