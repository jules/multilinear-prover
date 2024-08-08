#[cfg(test)]
mod tests {
    use crate::{
        field::{
            m31::{quartic::M31_4, M31},
            ChallengeField, Field,
        },
        iop::{sumcheck, zerocheck},
        linear_code::reed_solomon::ReedSolomonCode,
        pcs::{tensor_pcs::TensorPCS, PolynomialCommitmentScheme},
        polynomial::{MultilinearExtension, VirtualPolynomial},
        test_utils::{rand_poly, MockTranscript},
        transcript::Transcript,
    };
    use core::marker::PhantomData;

    fn zerocheck_test<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, T, E>,
    >(
        polys: &[MultilinearExtension<F>],
        transcript_p: &mut T,
        transcript_v: &mut T,
        pcs: PCS,
    ) -> bool {
        // Prover work
        let (sumcheck_claim, eval_point) = zerocheck::prove(polys, transcript_p);
        let mut poly = polys[0].clone();
        for p in &polys[1..] {
            poly.evals
                .iter_mut()
                .zip(p.evals.iter())
                .for_each(|(eval, m)| eval.mul_assign(m));
        }
        let commitment = pcs.commit(&[poly.clone()], transcript_p);
        let proof = pcs.prove(&commitment, &[poly], &eval_point, transcript_p);

        // Verifier work
        if let Ok(final_claim) = zerocheck::verify(sumcheck_claim, transcript_v) {
            pcs.verify(&commitment, &eval_point, final_claim, &proof, transcript_v)
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
        polys: &[MultilinearExtension<F>],
        transcript_p: &mut T,
        transcript_v: &mut T,
        pcs: PCS,
    ) -> bool {
        // Prover work
        let (sumcheck_claim, eval_point) = sumcheck::prove(polys, transcript_p);
        let mut poly = polys[0].clone();
        for p in &polys[1..] {
            poly.evals
                .iter_mut()
                .zip(p.evals.iter())
                .for_each(|(eval, m)| eval.mul_assign(m));
        }
        let commitment = pcs.commit(&[poly.clone()], transcript_p);
        let proof = pcs.prove(&commitment, &[poly], &eval_point, transcript_p);

        // Verifier work
        if let Ok(final_claim) = sumcheck::verify(sumcheck_claim, transcript_v) {
            pcs.verify(&commitment, &eval_point, final_claim, &proof, transcript_v)
        } else {
            false
        }
    }

    fn tensor_pcs_sumcheck<F: Field, E: ChallengeField<F>>(polys: &[MultilinearExtension<F>])
    where
        [(); F::NUM_BYTES_IN_REPR]:,
    {
        let tensor_pcs = TensorPCS::new(4);

        let mut transcript_p = MockTranscript::default();
        let mut transcript_v = MockTranscript::default();

        assert!(sumcheck_test::<
            _,
            E,
            _,
            TensorPCS<F, MockTranscript<F>, E, ReedSolomonCode<F, E>>,
        >(
            polys, &mut transcript_p, &mut transcript_v, tensor_pcs
        ));
    }

    #[test]
    fn tensor_pcs_1_poly_test() {
        tensor_pcs_sumcheck::<M31, M31_4>(&[rand_poly(2u32.pow(20) as usize)]);
    }

    #[test]
    fn tensor_pcs_64_poly_test() {
        tensor_pcs_sumcheck::<M31, M31_4>(
            &(0..64)
                .map(|_| rand_poly(2u32.pow(20) as usize))
                .collect::<Vec<MultilinearExtension<M31>>>(),
        );
    }
}
