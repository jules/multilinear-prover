pub mod reed_solomon;
pub use reed_solomon::ReedSolomonCode;

use crate::field::Field;

pub trait LinearCode<F: Field>: Send + Sync {
    const BLOWUP_BITS: usize;

    fn encode(&self, els: &[F]) -> Vec<F>;
}
