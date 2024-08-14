#[cfg(test)]
mod tests {
    use crate::{
        fft::CircleFFT,
        field::{
            m31::{quartic::M31_4, M31},
            ChallengeField, Field,
        },
        iop::{prodcheck, sumcheck, zerocheck},
        linear_code::reed_solomon::ReedSolomonCode,
        pcs::{tensor_pcs::TensorPCS, PolynomialCommitmentScheme},
        polynomial::{MultilinearExtension, MultivariatePolynomial, VirtualPolynomial},
        test_utils::{rand_poly, MockTranscript},
        transcript::Transcript,
    };

    const POLY_SIZE_BITS: u32 = 20;
    const ROOTS_OF_UNITY_BITS: usize = 12;

    fn prodcheck_test<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, T, E>,
    >(
        unsorted_columns: &[MultilinearExtension<F>],
        sorted_columns: &[MultilinearExtension<F>],
        transcript_p: &mut T,
        transcript_v: &mut T,
        pcs: PCS,
    ) -> bool {
        // Prover work
        let (sumcheck_claim, eval_point, zero_poly, prod_poly) =
            prodcheck::prove(unsorted_columns, sorted_columns, transcript_p);

        let zero_commitment = pcs.commit(&[zero_poly.clone()], transcript_p);
        let zero_proof = pcs.prove(
            &zero_commitment,
            &[zero_poly.clone()],
            &eval_point,
            transcript_p,
        );

        let prod_commitment = pcs.commit(&[prod_poly.clone().into()], transcript_p);
        let mut prod_eval = vec![E::ONE; prod_poly.num_vars()];
        prod_eval[0] = E::ZERO;
        let prod_proof = pcs.prove(
            &prod_commitment,
            &[prod_poly.clone().into()],
            &prod_eval,
            transcript_p,
        );

        // Verifier work
        if let Ok(final_claim) = zerocheck::verify(sumcheck_claim, transcript_v) {
            if final_claim != E::ZERO {
                return false;
            }
            pcs.verify(
                &zero_commitment,
                &eval_point,
                &[final_claim],
                &zero_proof,
                transcript_v,
            ) && pcs.verify(
                &prod_commitment,
                &prod_eval,
                &[E::ONE],
                &prod_proof,
                transcript_v,
            )
        } else {
            false
        }
    }

    fn zerocheck_test<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, T, E>,
    >(
        poly: &mut VirtualPolynomial<F>,
        transcript_p: &mut T,
        transcript_v: &mut T,
        pcs: PCS,
    ) -> bool {
        // Prover work
        let (sumcheck_claim, eval_point) = zerocheck::prove(poly, transcript_p);

        let commitment = pcs.commit(&[poly.clone()], transcript_p);
        let proof = pcs.prove(&commitment, &[poly.clone()], &eval_point, transcript_p);

        // Verifier work
        if let Ok(final_claim) = zerocheck::verify(sumcheck_claim, transcript_v) {
            if final_claim != E::ZERO {
                return false;
            }
            pcs.verify(
                &commitment,
                &eval_point,
                &[final_claim],
                &proof,
                transcript_v,
            )
        } else {
            false
        }
    }

    fn sumcheck_test<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, T, E>,
    >(
        poly: &VirtualPolynomial<F>,
        transcript_p: &mut T,
        transcript_v: &mut T,
        pcs: PCS,
    ) -> bool {
        // Prover work
        let (sumcheck_claim, eval_point) = sumcheck::prove(poly, transcript_p);
        let commitment = pcs.commit(&[poly.clone()], transcript_p);
        let proof = pcs.prove(&commitment, &[poly.clone()], &eval_point, transcript_p);

        // Verifier work
        if let Ok(final_claim) = sumcheck::verify(sumcheck_claim, transcript_v) {
            pcs.verify(
                &commitment,
                &eval_point,
                &[final_claim],
                &proof,
                transcript_v,
            )
        } else {
            false
        }
    }

    fn tensor_pcs_prodcheck(
        unsorted: &[MultilinearExtension<M31>],
        sorted: &[MultilinearExtension<M31>],
    ) -> bool {
        let tensor_pcs = TensorPCS::new(4, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));

        let mut transcript_p = MockTranscript::default();
        let mut transcript_v = MockTranscript::default();

        prodcheck_test::<
            M31,
            M31_4,
            MockTranscript<M31>,
            TensorPCS<M31, MockTranscript<M31>, M31_4, ReedSolomonCode<M31, CircleFFT>>,
        >(
            unsorted,
            sorted,
            &mut transcript_p,
            &mut transcript_v,
            tensor_pcs,
        )
    }

    fn tensor_pcs_zerocheck(polys: &[MultilinearExtension<M31>]) -> bool {
        let mut poly: VirtualPolynomial<M31> = polys[0].clone().into();
        for p in polys.iter().skip(1) {
            poly.mul_assign_mle(p);
        }

        let tensor_pcs = TensorPCS::new(4, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));

        let mut transcript_p = MockTranscript::default();
        let mut transcript_v = MockTranscript::default();

        zerocheck_test::<
            M31,
            M31_4,
            MockTranscript<M31>,
            TensorPCS<M31, MockTranscript<M31>, M31_4, ReedSolomonCode<M31, CircleFFT>>,
        >(&mut poly, &mut transcript_p, &mut transcript_v, tensor_pcs)
    }

    fn tensor_pcs_sumcheck(polys: &[MultilinearExtension<M31>]) {
        let mut poly: VirtualPolynomial<M31> = polys[0].clone().into();
        for p in polys.iter().skip(1) {
            poly.mul_assign_mle(p);
        }

        let tensor_pcs = TensorPCS::new(4, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));

        let mut transcript_p = MockTranscript::default();
        let mut transcript_v = MockTranscript::default();

        assert!(sumcheck_test::<
            M31,
            M31_4,
            MockTranscript<M31>,
            TensorPCS<M31, MockTranscript<M31>, M31_4, ReedSolomonCode<M31, CircleFFT>>,
        >(
            &poly, &mut transcript_p, &mut transcript_v, tensor_pcs
        ));
    }

    #[test]
    fn tensor_pcs_1_poly_test() {
        tensor_pcs_sumcheck(&[rand_poly(2u32.pow(POLY_SIZE_BITS) as usize)]);
    }

    #[test]
    fn tensor_pcs_64_poly_test() {
        tensor_pcs_sumcheck(
            &(0..64)
                .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
                .collect::<Vec<MultilinearExtension<M31>>>(),
        );
    }

    #[test]
    fn tensor_pcs_1_poly_test_zerocheck_fail() {
        assert!(!tensor_pcs_zerocheck(&[rand_poly(
            2u32.pow(POLY_SIZE_BITS) as usize
        )]));
    }

    #[test]
    fn tensor_pcs_64_poly_test_zerocheck_fail() {
        assert!(!tensor_pcs_zerocheck(
            &(0..64)
                .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
                .collect::<Vec<MultilinearExtension<M31>>>(),
        ));
    }

    #[test]
    fn tensor_pcs_1_poly_test_zerocheck() {
        assert!(tensor_pcs_zerocheck(&[MultilinearExtension::new(
            vec![M31::ZERO; 2u32.pow(POLY_SIZE_BITS) as usize]
        )]));
    }

    #[test]
    fn tensor_pcs_64_poly_test_zerocheck() {
        // Because we currently stand-in the constraint poly for a entrywise mult between all polys, we
        // can create a successful zerocheck by inputting one zero polynomial.
        let mut polys = (0..64)
            .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
            .collect::<Vec<MultilinearExtension<M31>>>();
        polys.push(MultilinearExtension::new(vec![
            M31::ZERO;
            2u32.pow(POLY_SIZE_BITS)
                as usize
        ]));

        assert!(tensor_pcs_zerocheck(&polys));
    }

    // TODO size discrepancies here cause the fft to break because the pcs cant have a power of 2
    // sized square matrix
    #[test]
    fn tensor_pcs_1_poly_test_prodcheck() {
        let poly = rand_poly(2u32.pow(POLY_SIZE_BITS - 1) as usize);
        let mut sorted = poly.clone();
        sorted.evals.sort();
        assert!(tensor_pcs_prodcheck(&[poly], &[sorted]));
    }

    #[test]
    fn tensor_pcs_1_poly_test_prodcheck_fail() {
        let poly = rand_poly(2u32.pow(POLY_SIZE_BITS - 1) as usize);
        let sorted = rand_poly(2u32.pow(POLY_SIZE_BITS - 1) as usize);
        assert!(!tensor_pcs_prodcheck(&[poly], &[sorted]));
    }
}
