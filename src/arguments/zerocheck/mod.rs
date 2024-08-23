pub mod prover;
pub mod verifier;

#[cfg(test)]
mod tests {
    use super::{prover::ZerocheckProver, verifier::ZerocheckVerifier};
    use crate::{
        fft::CircleFFT,
        field::{
            m31::{quartic::M31_4, M31},
            Field,
        },
        linear_code::ReedSolomonCode,
        pcs::TensorPCS,
        polynomial::{MultilinearExtension, MultivariatePolynomial},
        test_utils::rand_poly,
        transcript::Blake2sTranscript,
        univariate_utils::precompute_lagrange_coefficients,
    };

    const POLY_SIZE_BITS: u32 = 20;
    const ROOTS_OF_UNITY_BITS: usize = (POLY_SIZE_BITS / 2 + 1) as usize;

    fn make_prover(
        degree: usize,
    ) -> ZerocheckProver<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        let lagrange_coefficients = precompute_lagrange_coefficients(degree + 2);
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(100, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        ZerocheckProver::new(transcript, pcs, lagrange_coefficients)
    }

    fn make_verifier() -> ZerocheckVerifier<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(100, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        ZerocheckVerifier::new(transcript, pcs)
    }

    #[test]
    fn test_prove_verify_single_poly() {
        let poly = MultilinearExtension::new(vec![M31::ZERO; 2u32.pow(20) as usize]);
        let mut prover = make_prover(poly.degree());
        let proof = prover.prove(poly.into());

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }
}
