use crate::field::Field;

pub trait LinearCode<F: Field> {
    const BLOWUP: u32;

    fn encode(els: &[F]) -> Vec<F>;
}
