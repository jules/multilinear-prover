/// A multilinear polynomial represented in the Lagrange basis.
#[derive(Debug)]
pub struct MultilinearPoly<F: Field> {
    num_vars: usize,
    pub evals: Vec<F>,
}

impl<F: Field> MultilinearPoly<F> {
    pub fn new(evals: Vec<F>) -> Self {
        debug_assert!(evals.len().is_power_of_two());
        Self {
            num_vars: f64::log2(evals.len()) as usize,
            evals,
        }
    }

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns an evaluation on the hypercube.
    pub fn eval(&self, index: usize) -> F {
        debug_assert!(index < self.evals.len());
        self.evals[index]
    }

    /// Sets the first variable of the polynomial to `point`. Used for fixing sumcheck challenges
    /// into the polynomial.
    pub fn fix_variable(&self, point: F) -> Self {
        let mut new = MultilinearPoly {
            num_vars: self.num_vars - 1,
            evals: Vec::with_capacity(1 << (self.num_vars - 1)),
        };

        // TODO: par
        for i in 0..(1 << (new.num_vars)) {
            new.evals
                .push(self.evals[i] + (self.evals[(i << 1) + 1] - self.evals[i << 1]) * point);
        }

        new
    }
}
