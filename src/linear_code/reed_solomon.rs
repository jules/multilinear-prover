use super::*;
use crate::fft::FFT;
use core::marker::PhantomData;

/// Reed-Solomon encoding for polynomial evaluation vectors.
pub struct ReedSolomonCode<F: Field, D: FFT<F>> {
    fft: D,
    _marker: PhantomData<F>,
}

impl<F: Field, D: FFT<F>> ReedSolomonCode<F, D> {
    pub fn new(order_bits: usize) -> Self {
        let fft = D::new(order_bits, Self::BLOWUP_BITS);
        Self {
            fft,
            _marker: PhantomData::<F>,
        }
    }
}

impl<F: Field, D: FFT<F>> LinearCode<F> for ReedSolomonCode<F, D> {
    const BLOWUP_BITS: usize = 2;

    // Perform a simple low-degree extension with an FFT.
    #[inline(always)]
    fn encode(&self, els: &[F]) -> Vec<F> {
        self.fft.fft_extend(&els)
    }
}
