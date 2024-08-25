//! Provides an argument for prodcheck.

pub mod prover;
pub mod verifier;

#[cfg(test)]
mod tests {
    use super::{prover::ProdcheckProver, verifier::ProdcheckVerifier};
    use crate::{
        fft::CircleFFT,
        field::{
            m31::{quartic::M31_4, M31},
            Field,
        },
        linear_code::ReedSolomonCode,
        pcs::TensorPCS,
        polynomial::{MultilinearExtension, VirtualPolynomial},
        test_utils::rand_poly,
        transcript::Blake2sTranscript,
        univariate_utils::precompute_lagrange_coefficients,
    };

    const POLY_SIZE_BITS: u32 = 20;
    const ROOTS_OF_UNITY_BITS: usize = (POLY_SIZE_BITS / 2 + 1) as usize;

    fn make_prover(
        degree: usize,
    ) -> ProdcheckProver<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        // We do degree + 3 here since we need degree + 1 for the sumcheck and one extra for the
        // multiplication with eq, and one for the zerocheck polynomial.
        let lagrange_coefficients = precompute_lagrange_coefficients(degree + 3);
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(100, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        ProdcheckProver::new(transcript, pcs, lagrange_coefficients)
    }

    fn make_verifier() -> ProdcheckVerifier<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(100, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        ProdcheckVerifier::new(transcript, pcs)
    }

    #[test]
    fn test_prove_verify_single_poly() {
        let poly = rand_poly(POLY_SIZE_BITS);
        let mut poly_sorted = poly.clone();
        poly_sorted.evals.sort();

        let mut prover = make_prover(poly.degree());
        let proof = prover.prove(&[poly], &[poly_sorted]);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }
}
