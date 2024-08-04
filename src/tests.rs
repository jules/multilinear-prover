#[cfg(test)]
mod tests {
    use crate::{
        field::{
            m31::{quartic::M31_4, M31},
            ChallengeField, Field,
        },
        iop::sumcheck,
        linear_code::reed_solomon::ReedSolomonCode,
        mle::MultilinearExtension,
        pcs::{tensor_pcs::TensorPCS, PolynomialCommitmentScheme},
        transcript::Transcript,
    };
    use rand::Rng;
    use std::marker::PhantomData;

    #[derive(Default)]
    pub struct MockTranscript<F: Field> {
        counter: usize,
        _marker: PhantomData<F>,
    }

    impl<F: Field> Transcript<F> for MockTranscript<F> {
        fn draw_challenge(&mut self) -> F {
            self.counter += 1;
            F::from_usize(self.counter)
        }
        fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E {
            self.counter += 1;
            E::new(vec![F::from_usize(self.counter); E::DEGREE])
        }
        fn draw_bits(&mut self, bits: usize) -> usize {
            0
        }
        fn observe_witness(&mut self, witness: F) {}
        fn observe_witnesses(&mut self, witness: &[F]) {}
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
        let commitment = pcs.commit(&[poly.clone()]);
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

        let mut transcript_p = MockTranscript {
            counter: 1,
            _marker: PhantomData::<F>,
        };
        let mut transcript_v = MockTranscript {
            counter: 1,
            _marker: PhantomData::<F>,
        };

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
        let mut evals = vec![M31::default(); 2u32.pow(20) as usize];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = MultilinearExtension::new(evals);
        tensor_pcs_sumcheck::<M31, M31_4>(&[poly]);
    }

    #[test]
    fn tensor_pcs_64_poly_test() {
        tensor_pcs_sumcheck::<M31, M31_4>(
            &(0..64)
                .map(|_| {
                    let mut evals = vec![M31::default(); 2u32.pow(20) as usize];
                    evals
                        .iter_mut()
                        .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
                    MultilinearExtension::new(evals)
                })
                .collect::<Vec<MultilinearExtension<M31>>>(),
        );
    }
}
