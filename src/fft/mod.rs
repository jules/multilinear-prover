pub mod circle_fft;
pub use circle_fft::CircleFFT;

use crate::field::Field;

pub trait FFT<F: Field>: Send + Sync {
    fn new(order_bits: usize, blowup_bits: usize) -> Self;
    fn fft(&self, coeffs: &[F]) -> Vec<F>;
    fn ifft(&self, coeffs: &[F]) -> Vec<F>;
    fn lde(&self, coeffs: &[F]) -> Vec<F>;
}
