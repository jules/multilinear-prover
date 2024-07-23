use crate::field::{ChallengeField, Field};

pub trait Transcript<F: Field> {
    fn draw_challenge(&mut self) -> F;
    fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E;
    fn observe_witness(&mut self, witness: F);
    fn observe_witnesses(&mut self, witness: &[F]);
}
