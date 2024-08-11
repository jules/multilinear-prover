use super::*;
use crate::{circle_fft::*, univariate_utils::*};
use core::marker::PhantomData;

pub struct ReedSolomonCode<F: Field, E: ChallengeField<F>> {
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
}

impl<F: Field, E: ChallengeField<F>> LinearCode<F, E> for ReedSolomonCode<F, E> {
    const BLOWUP: usize = 4;

    // Since we only expect small vectors (the MLE is already turned into a square matrix for which
    // we encode each row separately) we simply multiply with a Vandermonde matrix to get the RS
    // encoding.
    fn encode(els: &[F]) -> Vec<F> {
        (0..(els.len() * Self::BLOWUP))
            .map(|i| univariate_eval(els, F::from_usize(i)))
            .collect::<Vec<F>>()
    }

    // TODO: this should be extricated since it's more specific to the tensor pcs
    fn encode_ext(els: &[E]) -> Vec<E> {
        (0..(els.len() * Self::BLOWUP))
            .map(|i| univariate_eval_ext(els, F::from_usize(i)))
            .collect::<Vec<E>>()
    }
}
