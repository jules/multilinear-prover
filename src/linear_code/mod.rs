pub mod reed_solomon;

use crate::field::{ChallengeField, Field};

pub trait LinearCode<F: Field, E: ChallengeField<F>> {
    const BLOWUP: usize;

    fn encode(els: &[F]) -> Vec<F>;
    fn encode_ext(els: &[E]) -> Vec<E>;
}
