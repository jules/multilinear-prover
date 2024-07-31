pub mod reed_solomon;

use crate::field::Field;

pub trait LinearCode<F: Field> {
    const BLOWUP: usize;

    fn encode(els: &[F]) -> Vec<F>;
}
