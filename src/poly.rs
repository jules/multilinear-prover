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

    pub fn eval(&self, index: usize) -> F {
        debug_assert!(index < self.evals.len());
        self.evals[index]
    }
}
