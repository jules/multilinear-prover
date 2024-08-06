//! Mersenne-31 circle group FFT. Used in a variety of polynomial commitment schemes for this
//! project.

// At the moment, I'm foregoing using the AIR compiler work on this - it seems highly specialized
// to the case at hand and for now I should be able to get away with something more rudimentary.
// However, it should be useful for optimization purposes to take optimizations from it later in
// order to reduce prover times (such as cache optimizations).

use crate::field::{
    m31::{complex::M31_2, M31},
    Field, TwoAdicField,
};

/// Precomputes an `order` sized multiplicative group of complex roots of unity, used for the
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
    let compressed_coeffs = coeffs
        .chunks(2)
        .map(|chunk| M31_2 {
            c0: chunk[0],
            c1: chunk[1],
        })
        .collect::<Vec<M31_2>>();
    debug_assert!(compressed_coeffs.len() * 2 == coeffs.len());

    // Now we can perform our transformation.
    let twiddles = &roots[..roots.len() / 2];
    panic!()
}

/// Inverse Fast Fourier transform, which allows us to turn a set of polynomial coefficients (in
/// the extension of M31) into a set of evaluations. The resulting vector will be in the base field
/// of M31.
pub fn ifft(coeffs: &[M31_2], roots: &[M31_2]) -> Vec<M31> {
    panic!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::rand_poly;

    #[test]
    fn test_precompute_roots() {
        // compute 2^10 roots.
        let order = 10;
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
