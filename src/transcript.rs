use crate::field::Field;

pub trait Transcript<F: Field> {
    fn draw_challenge(&mut self) -> F;
    fn observe_witness(&mut self, witness: F);
    fn observe_witnesses(&mut self, witness: &[F]);
}
