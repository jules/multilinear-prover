use super::*;
use core::marker::PhantomData;

pub struct ReedSolomonCode<F: Field> {
    _marker: PhantomData<F>,
}

impl<F: Field> LinearCode<F> for ReedSolomonCode<F> {
    const BLOWUP: usize = 4;

    fn encode(els: &[F]) -> Vec<F> {
        let mut els = els.clone().to_vec();
        els.extend(els.clone());
        els.extend(els.clone());
        els
    }
}
