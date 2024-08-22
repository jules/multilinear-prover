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
        test_utils::rand_poly,
        transcript::{Blake2sTranscript, Transcript},
        univariate_utils::precompute_lagrange_coefficients,
    };

    const POLY_SIZE_BITS: u32 = 20;
    const ROOTS_OF_UNITY_BITS: usize = (POLY_SIZE_BITS / 2 + 1) as usize;

    //fn prodcheck_test<
    //    F: Field,
    //    E: ChallengeField<F>,
    //    T: Transcript<F>,
    //    PCS: PolynomialCommitmentScheme<F, T, E>,
    //>(
    //    unsorted_columns: &[MultilinearExtension<F>],
    //    sorted_columns: &[MultilinearExtension<F>],
    //    transcript_p: &mut T,
    //    transcript_v: &mut T,
    //    pcs: PCS,
    //) -> bool {
    //    // Prover work
    //    let precomputed = precompute_lagrange_coefficients(4);
    //    let (sumcheck_claim, eval_point, zero_poly, prod_poly) =
    //        prodcheck::prove(unsorted_columns, sorted_columns, transcript_p, &precomputed);

    //    let zero_commitment = pcs.commit(&[zero_poly.clone()], transcript_p);
    //    let zero_proof = pcs.prove(
    //        &zero_commitment,
    //        &[zero_poly.clone()],
    //        &eval_point,
    //        transcript_p,
    //    );

    //    let prod_commitment = pcs.commit(&[prod_poly.clone().into()], transcript_p);
    //    let mut prod_eval = vec![E::ONE; prod_poly.num_vars()];
    //    prod_eval[0] = E::ZERO;
    //    let prod_proof = pcs.prove(
    //        &prod_commitment,
    //        &[prod_poly.clone().into()],
    //        &prod_eval,
    //        transcript_p,
    //    );

    //    // Verifier work
    //    if let Ok(final_claim) = zerocheck::verify(sumcheck_claim, transcript_v) {
    //        if final_claim != E::ZERO {
    //            return false;
    //        }
    //        pcs.verify(
    //            &zero_commitment,
    //            &eval_point,
    //            &[final_claim],
    //            &zero_proof,
    //            transcript_v,
    //        ) && pcs.verify(
    //            &prod_commitment,
    //            &prod_eval,
    //            &[E::ONE],
    //            &prod_proof,
    //            transcript_v,
    //        )
    //    } else {
    //        false
    //    }
    //}

    fn zerocheck_test<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, T, E>,
    >(
        poly: &mut VirtualPolynomial<F>,
        polys: &[MultilinearExtension<F>],
        transcript_p: &mut T,
        transcript_v: &mut T,
        pcs: PCS,
    ) -> bool {
        // Prover work
        let precomputed = precompute_lagrange_coefficients(poly.degree() + 2);
        let now = std::time::Instant::now();
        let (sumcheck_claim, eval_point) = zerocheck::prove(poly, transcript_p, &precomputed);
        let elapsed = std::time::Instant::now();
        println!("zerocheck {:?}", elapsed - now);

        let now = std::time::Instant::now();
        let commitment = pcs.commit(&polys, transcript_p);
        let elapsed = std::time::Instant::now();
        println!("commit {:?}", elapsed - now);
        let evaluations = polys
            .iter()
            .map(|v| {
                let mut v_lifted = v.fix_variable_ext(eval_point[0]);
                for i in 1..20 {
                    v_lifted.fix_variable(eval_point[i]);
                }
                assert!(v_lifted.evals.len() == 1);
                v_lifted.evals[0]
            })
            .collect::<Vec<E>>();
        let now = std::time::Instant::now();
        let proof = pcs.prove(&commitment, &polys, &eval_point, transcript_p);
        let elapsed = std::time::Instant::now();
        println!("prove {:?}", elapsed - now);

        // Verifier work
        if let Ok((final_claim, mut eq)) = zerocheck::verify(sumcheck_claim, transcript_v) {
            let mut final_sum = final_claim.0;
            let challenge_point = final_claim.1;
            assert!(eq.evals == poly.constituents.last().unwrap().evals);
            assert!(challenge_point == eval_point);
            let mut eq_lifted = eq.fix_variable_ext(challenge_point[0]);
            for v in challenge_point.iter().skip(1) {
                eq_lifted.fix_variable(*v);
            }
            assert!(eq_lifted.evals.len() == 1);
            let mut final_value = evaluations.iter().fold(E::ONE, |mut acc, x| {
                acc.mul_assign(x);
                acc
            });
            final_value.mul_assign(&eq_lifted.evals[0]);

            if final_sum != final_value {
                return false;
            }
            pcs.verify(&commitment, &eval_point, &evaluations, &proof, transcript_v)
        } else {
            false
        }
    }

    //fn sumcheck_test<
    //    F: Field,
    //    E: ChallengeField<F>,
    //    T: Transcript<F>,
    //    PCS: PolynomialCommitmentScheme<F, T, E>,
    //>(
    //    poly: &VirtualPolynomial<F>,
    //    transcript_p: &mut T,
    //    transcript_v: &mut T,
    //    pcs: PCS,
    //) -> bool {
    //    // Prover work
    //    let precomputed = precompute_lagrange_coefficients(poly.degree() + 1);
    //    let (sumcheck_claim, eval_point) = sumcheck::prove(poly, transcript_p, &precomputed);
    //    let commitment = pcs.commit(&[poly.clone()], transcript_p);
    //    let proof = pcs.prove(&commitment, &[poly.clone()], &eval_point, transcript_p);

    //    // Verifier work
    //    if let Ok(final_claim) = sumcheck::verify(sumcheck_claim, transcript_v) {
    //        pcs.verify(
    //            &commitment,
    //            &eval_point,
    //            &[final_claim],
    //            &proof,
    //            transcript_v,
    //        )
    //    } else {
    //        false
    //    }
    //}

    //fn tensor_pcs_prodcheck(
    //    unsorted: &[MultilinearExtension<M31>],
    //    sorted: &[MultilinearExtension<M31>],
    //) -> bool {
    //    let tensor_pcs = TensorPCS::new(4, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));

    //    let mut transcript_p = Blake2sTranscript::default();
    //    let mut transcript_v = Blake2sTranscript::default();

    //    prodcheck_test::<
    //        M31,
    //        M31_4,
    //        Blake2sTranscript<M31>,
    //        TensorPCS<M31, Blake2sTranscript<M31>, M31_4, ReedSolomonCode<M31, CircleFFT>>,
    //    >(
    //        unsorted,
    //        sorted,
    //        &mut transcript_p,
    //        &mut transcript_v,
    //        tensor_pcs,
    //    )
    //}

    fn tensor_pcs_zerocheck(polys: &[MultilinearExtension<M31>]) -> bool {
        let mut poly: VirtualPolynomial<M31> = polys[0].clone().into();
        for p in polys.iter().skip(1) {
            poly.mul_assign_mle(p);
        }

        let tensor_pcs = TensorPCS::new(100, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));

        let mut transcript_p = Blake2sTranscript::default();
        let mut transcript_v = Blake2sTranscript::default();

        zerocheck_test::<
            M31,
            M31_4,
            Blake2sTranscript<M31>,
            TensorPCS<M31, Blake2sTranscript<M31>, M31_4, ReedSolomonCode<M31, CircleFFT>>,
        >(
            &mut poly,
            polys,
            &mut transcript_p,
            &mut transcript_v,
            tensor_pcs,
        )
    }

    //fn tensor_pcs_sumcheck(polys: &[MultilinearExtension<M31>]) {
    //    let mut poly: VirtualPolynomial<M31> = polys[0].clone().into();
    //    for p in polys.iter().skip(1) {
    //        poly.mul_assign_mle(p);
    //    }

    //    let tensor_pcs = TensorPCS::new(4, ReedSolomonCode::new(ROOTS_OF_UNITY_BITS));

    //    let mut transcript_p = Blake2sTranscript::default();
    //    let mut transcript_v = Blake2sTranscript::default();

    //    assert!(sumcheck_test::<
    //        M31,
    //        M31_4,
    //        Blake2sTranscript<M31>,
    //        TensorPCS<M31, Blake2sTranscript<M31>, M31_4, ReedSolomonCode<M31, CircleFFT>>,
    //    >(
    //        &poly, &mut transcript_p, &mut transcript_v, tensor_pcs
    //    ));
    //}

    //#[test]
    //fn tensor_pcs_1_poly_test() {
    //    tensor_pcs_sumcheck(&[rand_poly(2u32.pow(POLY_SIZE_BITS) as usize)]);
    //}

    //#[test]
    //fn tensor_pcs_64_poly_test() {
    //    tensor_pcs_sumcheck(
    //        &(0..64)
    //            .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
    //            .collect::<Vec<MultilinearExtension<M31>>>(),
    //    );
    //}

    #[test]
    fn tensor_pcs_1_poly_test_zerocheck_fail() {
        assert!(!tensor_pcs_zerocheck(&[rand_poly(
            2u32.pow(POLY_SIZE_BITS) as usize
        )]));
    }

    //#[test]
    //fn tensor_pcs_64_poly_test_zerocheck_fail() {
    //    assert!(!tensor_pcs_zerocheck(
    //        &(0..64)
    //            .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
    //            .collect::<Vec<MultilinearExtension<M31>>>(),
    //    ));
    //}

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
        let mut polys = (0..63)
            .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
            .collect::<Vec<MultilinearExtension<M31>>>();
        polys.push(MultilinearExtension::new(vec![
            M31::ZERO;
            2u32.pow(POLY_SIZE_BITS)
                as usize
        ]));

        assert!(tensor_pcs_zerocheck(&polys));
    }

    #[test]
    fn tensor_pcs_256() {
        // Because we currently stand-in the constraint poly for a entrywise mult between all polys, we
        // can create a successful zerocheck by inputting one zero polynomial.
        let mut polys = (0..255)
            .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
            .collect::<Vec<MultilinearExtension<M31>>>();
        polys.push(MultilinearExtension::new(vec![
            M31::ZERO;
            2u32.pow(POLY_SIZE_BITS)
                as usize
        ]));

        tensor_pcs_zerocheck(&polys);
    }

    //#[test]
    //fn tensor_pcs_1_poly_test_prodcheck() {
    //    let poly = rand_poly(2u32.pow(POLY_SIZE_BITS) as usize);
    //    let mut sorted = poly.clone();
    //    sorted.evals.sort();
    //    assert!(tensor_pcs_prodcheck(&[poly], &[sorted]));
    //}

    //#[test]
    //fn tensor_pcs_1_poly_test_prodcheck_fail() {
    //    let poly = rand_poly(2u32.pow(POLY_SIZE_BITS) as usize);
    //    let sorted = rand_poly(2u32.pow(POLY_SIZE_BITS) as usize);
    //    assert!(!tensor_pcs_prodcheck(&[poly], &[sorted]));
    //}

    //#[test]
    //fn tensor_pcs_8_poly_test_prodcheck() {
    //    let polys = (0..8)
    //        .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
    //        .collect::<Vec<MultilinearExtension<M31>>>();
    //    let mut sorted = polys.clone();
    //    sorted.iter_mut().for_each(|poly| {
    //        poly.evals.sort();
    //    });
    //    assert!(tensor_pcs_prodcheck(&polys, &sorted));
    //}

    //#[test]
    //fn tensor_pcs_8_poly_test_prodcheck_fail() {
    //    let polys = (0..8)
    //        .map(|_| rand_poly(2u32.pow(POLY_SIZE_BITS) as usize))
    //        .collect::<Vec<MultilinearExtension<M31>>>();
    //    let mut sorted = polys.clone();
    //    sorted[0] = rand_poly(2u32.pow(POLY_SIZE_BITS) as usize); // set one polynomial to be
    //                                                              // different
    //    sorted.iter_mut().for_each(|poly| {
    //        poly.evals.sort();
    //    });
    //    assert!(!tensor_pcs_prodcheck(&polys, &sorted));
    //}
}
