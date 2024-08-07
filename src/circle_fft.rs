//! Mersenne-31 circle group FFT. Used in a variety of polynomial commitment schemes for this
//! project.

// At the moment, I'm foregoing using the AIR compiler work on this - it seems highly specialized
// to the case at hand and for now I should be able to get away with something more rudimentary.
// However, it should be useful for optimization purposes to take optimizations from it later in
// order to reduce prover times (such as cache optimizations).
//
// Most of the code here is taken from Plonky3 https://github.com/Plonky3/Plonky3/blob/main/mersenne-31/.

use crate::field::{
    m31::{complex::M31_2, M31},
    ChallengeField, Field, PrimeField, TwoAdicField,
};

const ORDER_SQ: i64 = (M31::ORDER as i64) * (M31::ORDER as i64);
const TWO_ORDER_SQ: i64 = 2 * ORDER_SQ;

fn bitreverse_idx(mut n: u32, l: u32) -> u32 {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }

    r
}

pub fn bit_reverse(a: &mut [M31_2]) {
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
    for _ in 0..(1 << order_bits) {
        roots.push(root);
        root.mul_assign(&base_root); // NB last one not necessary
    }

    debug_assert!(*roots.last().unwrap() == M31_2::ONE);
    roots
}

/// Fast Fourier Transform, allowing us to interpret the coefficients as evaluations of a
/// polynomial, and interpolating the polynomial in monomial basis. Due to usage of the circle
/// group FFT, the resulting vector of coefficients resides in the complex extension of M31.
pub fn fft(coeffs: &[M31], roots: &[M31_2]) -> Vec<M31_2> {
    debug_assert!(roots.len() * 2 == coeffs.len());

    // Preprocess the coefficients; this means we need to pack them into complex field elements.
    let mut compressed_coeffs = coeffs
        .chunks(2)
        .map(|chunk| M31_2 {
            c0: chunk[0],
            c1: chunk[1],
        })
        .collect::<Vec<M31_2>>();
    debug_assert!(compressed_coeffs.len() * 2 == coeffs.len());

    fft_inner(&mut compressed_coeffs, roots);

    // NB postprocess?
    compressed_coeffs
}

fn fft_inner(coeffs: &mut [M31_2], roots: &[M31_2]) {
    // Now we can perform our transformation.
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

            coeffs[hi] = M31_2::new(vec![a1, a2]);
            coeffs[lo] = M31_2::new(vec![b1, b2]);
        }
    }
}

/// Inverse Fast Fourier transform, which allows us to turn a set of polynomial coefficients (in
/// the extension of M31) into a set of evaluations. The resulting vector will be in the base field
/// of M31.
pub fn ifft(coeffs: &[M31_2], roots: &[M31_2]) -> Vec<M31> {
    let mut inverted_coeffs = coeffs.clone().to_vec();
    fft_inner(&mut inverted_coeffs, roots);

    // Scale result
    let len_inv = M31_2::from_usize(inverted_coeffs.len()).inverse().unwrap();
    inverted_coeffs
        .iter_mut()
        .for_each(|coeff| coeff.mul_assign(&len_inv));

    // Postprocessing and unpacking.
    bit_reverse(&mut inverted_coeffs);
    inverted_coeffs
        .iter()
        .flat_map(|v| [v.c0, v.c1])
        .collect::<Vec<M31>>()
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
        let order = 9;
        let roots = precompute_roots(order);
        let poly = rand_poly::<M31>(2u32.pow(10) as usize);
        let result = ifft(&fft(&poly.evals, &roots), &roots);
        assert!(poly
            .evals
            .iter()
            .zip(result.iter())
            .all(|(v, eval)| v == eval));
    }

    #[test]
    fn test_lde() {
        let order = 9;
        let roots = precompute_roots(order);
        let blowup_order = 10;
        let blowup_roots = precompute_roots(order);
        let poly = rand_poly::<M31>(2u32.pow(10) as usize);
        let result = ifft(&fft(&poly.evals, &roots), &blowup_roots);
        assert!(poly
            .evals
            .iter()
            .zip(result.iter())
            .all(|(v, eval)| v == eval));
    }
}
