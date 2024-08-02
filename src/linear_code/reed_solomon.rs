use super::*;
use core::marker::PhantomData;

pub struct ReedSolomonCode<F: Field, E: ChallengeField<F>> {
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
}

impl<F: Field, E: ChallengeField<F>> LinearCode<F, E> for ReedSolomonCode<F, E> {
    const BLOWUP: usize = 4;

    fn encode(els: &[F]) -> Vec<F> {
        let mut els = els.clone().to_vec();
        els.extend(els.clone());
        els.extend(els.clone());
        els
    }

    fn encode_ext(els: &[E]) -> Vec<E> {
        let mut els = els.clone().to_vec();
        els.extend(els.clone());
        els.extend(els.clone());
        els
    }
}
