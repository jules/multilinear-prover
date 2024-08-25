//! Provides an argument for zerocheck.

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
        polynomial::{MultilinearExtension, VirtualPolynomial},
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
        // We do degree + 2 here since we need degree + 1 for the sumcheck and one extra for the
        // multiplication with eq.
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
        let poly = MultilinearExtension::new(vec![M31::ZERO; 2u32.pow(POLY_SIZE_BITS) as usize]);

        let mut prover = make_prover(poly.degree());
        let proof = prover.prove(poly.into());

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }

    #[test]
    fn test_prove_mult_gate() {
        let multipliers = (0..2)
            .map(|_| rand_poly(POLY_SIZE_BITS))
            .collect::<Vec<MultilinearExtension<M31>>>();
        let answers = MultilinearExtension::new(
            multipliers[0]
                .evals
                .iter()
                .zip(multipliers[1].evals.iter())
                .map(|(a, b)| {
                    let mut a = a.clone();
                    a.mul_assign(b);
                    a
                })
                .collect::<Vec<M31>>(),
        );

        // a * b - c
        let mut poly: VirtualPolynomial<M31> = multipliers[0].clone().into();
        poly.mul_assign_mle(&multipliers[1]);
        poly.add_assign_mle(&answers, 1, true);

        let mut prover = make_prover(poly.degree());
        let proof = prover.prove(poly);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }

    #[test]
    fn test_prove_mult_gate_fail() {
        let multipliers = (0..2)
            .map(|_| rand_poly(POLY_SIZE_BITS))
            .collect::<Vec<MultilinearExtension<M31>>>();
        let wrong_answers = MultilinearExtension::new(
            multipliers[0]
                .evals
                .iter()
                .zip(multipliers[1].evals.iter())
                .map(|(a, b)| {
                    let mut a = a.clone();
                    a.mul_assign(b);
                    a.add_assign(&M31::ONE);
                    a
                })
                .collect::<Vec<M31>>(),
        );

        // a * b - c
        let mut poly: VirtualPolynomial<M31> = multipliers[0].clone().into();
        poly.mul_assign_mle(&multipliers[1]);
        poly.add_assign_mle(&wrong_answers, 1, true);

        let mut prover = make_prover(poly.degree());
        let proof = prover.prove(poly);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(!accept);
    }

    #[test]
    fn test_prove_verify_4_poly() {
        let mut poly_set = (0..3)
            .map(|_| rand_poly(POLY_SIZE_BITS))
            .collect::<Vec<MultilinearExtension<M31>>>();
        let zero_poly =
            MultilinearExtension::new(vec![M31::ZERO; 2u32.pow(POLY_SIZE_BITS) as usize]);
        poly_set.push(zero_poly);

        let mut poly: VirtualPolynomial<M31> = poly_set[0].clone().into();
        poly_set.iter().skip(1).for_each(|p| {
            poly.mul_assign_mle(p);
        });

        let mut prover = make_prover(poly.degree());
        let proof = prover.prove(poly);

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }
}
