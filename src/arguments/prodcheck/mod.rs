//! Provides an argument for prodcheck.

pub mod prover;
pub mod verifier;

#[cfg(test)]
mod tests {
    use super::{prover::ProdcheckProver, verifier::ProdcheckVerifier};
    use crate::{
        fft::CircleFFT,
        field::m31::{quartic::M31_4, M31},
        linear_code::ReedSolomonCode,
        pcs::TensorPCS,
        polynomial::MultilinearExtension,
        test_utils::rand_poly,
        transcript::Blake2sTranscript,
        univariate_utils::precompute_lagrange_coefficients,
    };

    const POLY_SIZE_BITS: u32 = 20;
    const ROOTS_OF_UNITY_BITS: usize = (POLY_SIZE_BITS / 2 + 1) as usize;

    fn make_prover() -> ProdcheckProver<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        // We do 4 here since we need degree + 1 for the sumcheck and one extra for the
        // multiplication with eq, and one for the zerocheck polynomial. We never increase our
        // sumcheck degree so we can hardcode at 4.
        let lagrange_coefficients = precompute_lagrange_coefficients(4);
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(143, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        ProdcheckProver::new(transcript, pcs, lagrange_coefficients)
    }

    fn make_verifier() -> ProdcheckVerifier<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(143, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        ProdcheckVerifier::new(transcript, pcs)
    }

    #[test]
    fn test_prove_verify_single_poly() {
        let poly = rand_poly(POLY_SIZE_BITS);
        let mut poly_sorted = poly.clone();
        poly_sorted.evals.sort();

        let mut prover = make_prover();
        let proof = prover.prove(&[poly], &[poly_sorted]);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }

    #[test]
    fn test_prove_verify_faulty_product() {
        let poly = rand_poly(POLY_SIZE_BITS);
        let poly_sorted = rand_poly(POLY_SIZE_BITS);

        let mut prover = make_prover();
        let proof = prover.prove(&[poly], &[poly_sorted]);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(!accept);
    }

    #[test]
    fn test_prove_verify_4_poly() {
        let unsorted_polys = (0..4)
            .map(|_| rand_poly(POLY_SIZE_BITS))
            .collect::<Vec<MultilinearExtension<M31>>>();
        let mut sorted_polys = unsorted_polys.clone();
        sorted_polys.iter_mut().for_each(|p| p.evals.sort());

        let mut prover = make_prover();
        let proof = prover.prove(&unsorted_polys, &sorted_polys);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }

    #[test]
    fn test_prove_verify_4_poly_faulty() {
        let unsorted_polys = (0..4)
            .map(|_| rand_poly(POLY_SIZE_BITS))
            .collect::<Vec<MultilinearExtension<M31>>>();
        let mut sorted_polys = unsorted_polys.clone();
        sorted_polys[0] = rand_poly(POLY_SIZE_BITS);
        sorted_polys.iter_mut().for_each(|p| p.evals.sort());

        let mut prover = make_prover();
        let proof = prover.prove(&unsorted_polys, &sorted_polys);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(!accept);
    }
}
