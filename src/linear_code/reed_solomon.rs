use super::*;
use crate::univariate_utils::*;
use core::marker::PhantomData;
use rayon::prelude::*;
use std::alloc::Global;

pub struct ReedSolomonCode<F: Field, E: ChallengeField<F>> {
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
}

impl<F: Field, E: ChallengeField<F>> LinearCode<F, E> for ReedSolomonCode<F, E> {
    const BLOWUP: usize = 2;

    // Since we only expect small vectors (the MLE is already turned into a square matrix for which
    // we encode each row separately) we simply multiply with a Vandermonde matrix to get the RS
    // encoding.
    // TODO: FFT WILL MEGA SPEED THIS UP - implement circle group FFT and change the transformation
    // to be like an LDE rather than just a forward transform so that we stay in the base field
    fn encode(els: &[F]) -> Vec<F> {
        (0..(els.len() * Self::BLOWUP))
            .map(|i| univariate_eval(els, F::from_usize(i)))
            .collect::<Vec<F>>()
    }

    // TODO how do we circle group fft with this - i assume we lift into the octal extension?
    fn encode_ext(els: &[E]) -> Vec<E> {
        (0..(els.len() * Self::BLOWUP))
            .map(|i| univariate_eval_ext(els, F::from_usize(i)))
            .collect::<Vec<E>>()
    }
}
