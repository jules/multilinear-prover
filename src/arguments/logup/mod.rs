//! Provides an argument for logup.

pub mod prover;
pub mod verifier;

#[cfg(test)]
mod tests {
    use super::{prover::LogUpProver, verifier::LogUpVerifier};
    use crate::{
        fft::CircleFFT,
        field::{
            m31::{quartic::M31_4, M31},
            Field,
        },
        linear_code::ReedSolomonCode,
        pcs::TensorPCS,
        polynomial::MultilinearExtension,
        test_utils::rand_poly,
        transcript::Blake2sTranscript,
        univariate_utils::precompute_lagrange_coefficients,
    };
    use rand::Rng;

    const POLY_SIZE_BITS: u32 = 4;
    const ROOTS_OF_UNITY_BITS: usize = (POLY_SIZE_BITS / 2 + 1) as usize;

    fn make_prover() -> LogUpProver<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        // We can hardcode at 4 for now.
        let lagrange_coefficients = precompute_lagrange_coefficients(4);
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(143, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        LogUpProver::new(transcript, pcs, lagrange_coefficients)
    }

    fn make_verifier() -> LogUpVerifier<
        M31,
        M31_4,
        Blake2sTranscript<M31>,
        TensorPCS<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>>,
    > {
        let transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::new(143, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));
        LogUpVerifier::new(transcript, pcs)
    }

    #[test]
    fn test_logup() {
        let table = (1..=1 << POLY_SIZE_BITS)
            .map(|i| M31::from_usize(i))
            .collect::<Vec<M31>>();

        let mut multiplicities = vec![M31::ZERO; table.len()];
        let trace_columns = (0..3)
            .map(|_| {
                MultilinearExtension::new(
                    (0..table.len())
                        .map(|_| {
                            let r = rand::thread_rng().gen_range(0..table.len());
                            multiplicities[r].add_assign(&M31::ONE);
                            table[r]
                        })
                        .collect::<Vec<M31>>(),
                )
            })
            .collect::<Vec<MultilinearExtension<M31>>>();

        let mut prover = make_prover();
        let proof = prover.prove(
            &trace_columns,
            &MultilinearExtension::new(table),
            &MultilinearExtension::new(multiplicities),
        );

        let mut verifier = make_verifier();
        let accept = verifier
            .verify(&proof)
            .expect("should not panic after verifying a well-constructed proof");
        assert!(accept);
    }
}
