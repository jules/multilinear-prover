//! Mersenne-31 circle group FFT. Used in a variety of polynomial commitment schemes for this
//! project.

use super::FFT;
use crate::field::{
    m31::{complex::M31_2, M31},
    Field, PrimeField, TwoAdicField,
};

const ORDER_SQ: i64 = (M31::ORDER as i64) * (M31::ORDER as i64);
const TWO_ORDER_SQ: i64 = 2 * ORDER_SQ;

pub struct CircleFFT {
    roots: Vec<M31>,
    blowup_roots: Vec<M31>,
    blowup_bits: usize,
}

impl FFT<M31> for CircleFFT {
    fn new(order_bits: usize, blowup_bits: usize) -> Self {
        // We chop off one bit from the order since we are making complex roots of unity, not real
        // ones.
        let roots = precompute_roots(order_bits - 1)
            .iter()
            .flat_map(|v| [v.c0, v.c1])
            .collect::<Vec<M31>>();
        let blowup_roots = precompute_roots(order_bits - 1 + blowup_bits)
            .iter()
            .flat_map(|v| [v.c0, v.c1])
            .collect::<Vec<M31>>();

        Self {
            roots,
            blowup_roots,
            blowup_bits,
        }
    }

    fn fft(&self, coeffs: &[M31]) -> Vec<M31> {
        let mut packed = preprocess(coeffs);
        let packed_roots = preprocess(&self.roots);
        fft(&mut packed, &packed_roots);
        postprocess(&packed)
    }

    fn extend(&self, coeffs: &[M31]) -> Vec<M31> {
        let mut packed = preprocess(coeffs);
        packed.resize(packed.len() << self.blowup_bits, M31_2::ZERO);
        let packed_blowup_roots = preprocess(&self.blowup_roots);
        fft(&mut packed, &packed_blowup_roots);
        postprocess(&packed)
    }

    fn ifft(&self, coeffs: &[M31]) -> Vec<M31> {
        let mut packed = preprocess(coeffs);
        let packed_roots = preprocess(&self.roots);
        ifft(&mut packed, &packed_roots);
        postprocess(&packed)
    }

    fn lde(&self, coeffs: &[M31]) -> Vec<M31> {
        let mut packed = preprocess(coeffs);
        let packed_roots = preprocess(&self.roots);
        ifft(&mut packed, &packed_roots);
        packed.resize(packed.len() << self.blowup_bits, M31_2::ZERO);
        let packed_blowup_roots = preprocess(&self.blowup_roots);
        fft(&mut packed, &packed_blowup_roots);
        postprocess(&packed)
    }
}

fn bitreverse_idx(mut n: u32, l: u32) -> u32 {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }

    r
}

fn bit_reverse(a: &mut [M31_2]) {
    let n = a.len() as u32;
    debug_assert!(n.is_power_of_two());
    let log_n = n.trailing_zeros();

    for k in 0..n {
        let rk = bitreverse_idx(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }
}

/// Precomputes a 2^order_bits sized multiplicative group of complex roots of unity, used for the
/// Circle FFT in M31.
pub fn precompute_roots(order_bits: usize) -> Vec<M31_2> {
    let mut root = M31_2::two_adic_generator();
    root.exp_power_of_2(M31_2::TWO_ADICITY - order_bits);
    let base_root = root.clone();
    debug_assert!(base_root.pow(1 << order_bits) == M31_2::ONE);
    let mut roots = Vec::with_capacity(1 << order_bits);
    for i in 0..(1 << order_bits) {
        if i == (1 << order_bits) - 1 {
            roots.insert(0, root);
        } else {
            roots.push(root);
            root.mul_assign(&base_root);
        }
    }

    debug_assert!(*roots.first().unwrap() == M31_2::ONE);
    roots
}

pub fn preprocess(coeffs: &[M31]) -> Vec<M31_2> {
    coeffs
        .chunks(2)
        .map(|chunk| M31_2 {
            c0: chunk[0],
            c1: chunk[1],
        })
        .collect::<Vec<M31_2>>()
}

pub fn postprocess(coeffs: &[M31_2]) -> Vec<M31> {
    coeffs
        .iter()
        .flat_map(|v| [v.c0, v.c1])
        .collect::<Vec<M31>>()
}

/// Fast Fourier Transform, allowing us to interpret the coefficients of a polynomial, and
/// evaluate the polynomial. Due to usage of the circle group FFT, the resulting vector of
/// coefficients resides in the complex extension of M31.
pub fn fft(coeffs: &mut [M31_2], roots: &[M31_2]) {
    debug_assert!(roots.len() >= coeffs.len());
    let twiddles = &roots[..roots.len() / 2];
    let log_len = coeffs.len().ilog2() as usize;
    bit_reverse(coeffs);

    for i in 0..log_len {
        dit_layer(coeffs, i, log_len - 1 - i, &twiddles);
    }
}

fn dit_layer(coeffs: &mut [M31_2], layer: usize, rev_layer: usize, twiddles: &[M31_2]) {
    let half_block_size = 1 << layer;
    let block_size = 2 * half_block_size;

    for j in (0..coeffs.len()).step_by(block_size) {
        // Unrolled 0 case
        let hi = j;
        let lo = hi + half_block_size;
        let mut sum = coeffs[hi].clone();
        sum.add_assign(&coeffs[lo]);
        let mut diff = coeffs[hi].clone();
        diff.sub_assign(&coeffs[lo]);
        coeffs[hi] = sum;
        coeffs[lo] = diff;

        for i in 1..half_block_size {
            let hi = j + i;
            let lo = hi + half_block_size;
            let twiddle = twiddles[i << rev_layer];

            let (x1, x2) = (coeffs[hi].c0.0 as i64, coeffs[hi].c1.0 as i64);
            let (y1, y2) = (coeffs[lo].c0.0 as i64, coeffs[lo].c1.0 as i64);
            let (w1, w2) = (twiddle.c0.0 as i64, twiddle.c1.0 as i64);
            let z1 = y1 * w1 - y2 * w2;
            let a1 = M31::from_u64_with_reduction((ORDER_SQ + x1 + z1) as u64);
            let b1 = M31::from_u64_with_reduction((ORDER_SQ + x1 - z1) as u64);

            let z2 = y2 * w1 + y1 * w2;
            let a2 = M31::from_u64_with_reduction((x2 + z2) as u64);
            let b2 = M31::from_u64_with_reduction((TWO_ORDER_SQ + x2 - z2) as u64);

            coeffs[hi] = M31_2 { c0: a1, c1: a2 };
            coeffs[lo] = M31_2 { c0: b1, c1: b2 };
        }
    }
}

/// Inverse Fast Fourier transform, which allows us to turn a set of polynomial evaluations (in
/// the extension of M31) into a set of coefficients. The resulting vector will be in the base field
/// of M31.
pub fn ifft(coeffs: &mut [M31_2], roots: &[M31_2]) {
    fft(coeffs, roots);

    // Scale result
    let len_inv = M31_2::from_usize(coeffs.len()).inverse().unwrap();
    coeffs
        .iter_mut()
        .for_each(|coeff| coeff.mul_assign(&len_inv));

    // Postprocessing and unpacking.
    let h = coeffs.len();
    for i in 1..h / 2 {
        coeffs.swap(i, h - i);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::rand_poly;

    #[test]
    fn test_precompute_roots() {
        // compute 2^20 roots.
        let order = 20;
        let roots = precompute_roots(order);
        assert!(roots.len() == 1 << order);
    }

    #[test]
    fn test_fft_roundtrip() {
        let order = 19;
        let roots = precompute_roots(order);
        let poly = rand_poly::<M31>(20);
        let mut preprocessed = preprocess(&poly.evals);
        fft(&mut preprocessed, &roots);
        ifft(&mut preprocessed, &roots);
        let result = postprocess(&preprocessed);
        assert!(
            poly.evals.iter().fold(M31::ZERO, |mut acc, x| {
                acc.add_assign(x);
                acc
            }) == result.iter().fold(M31::ZERO, |mut acc, x| {
                acc.add_assign(x);
                acc
            })
        );
        assert!(poly
            .evals
            .iter()
            .zip(result.iter())
            .all(|(v, eval)| v == eval));
    }

    #[test]
    fn test_ifft_roundtrip() {
        let roots = precompute_roots(19);
        let poly = rand_poly::<M31>(20);
        let mut preprocessed = preprocess(&poly.evals);
        fft(&mut preprocessed, &roots);
        ifft(&mut preprocessed, &roots);
        let result = postprocess(&preprocessed);
        assert!(
            poly.evals.iter().fold(M31::ZERO, |mut acc, x| {
                acc.add_assign(x);
                acc
            }) == result.iter().fold(M31::ZERO, |mut acc, x| {
                acc.add_assign(x);
                acc
            })
        );
        assert!(poly
            .evals
            .iter()
            .zip(result.iter())
            .all(|(v, eval)| v == eval));
    }

    #[test]
    fn test_lde() {
        let roots = precompute_roots(3);
        let blowup_roots = precompute_roots(4);
        let poly = rand_poly::<M31>(4);
        let mut preprocessed = preprocess(&poly.evals);
        ifft(&mut preprocessed, &roots);
        preprocessed.resize(preprocessed.len() << 1, M31_2::ZERO);
        fft(&mut preprocessed, &blowup_roots);
        let result = postprocess(&preprocessed);
        assert!(poly
            .evals
            .iter()
            .all(|v| result.iter().find(|&&r| *v == r).is_some()));
    }
}
