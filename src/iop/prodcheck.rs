//! Product check IOP

use crate::{
    field::{ChallengeField, Field},
    iop::{sumcheck::SumcheckProof, zerocheck},
    polynomial::{MultilinearExtension, MultivariatePolynomial, VirtualPolynomial},
    transcript::Transcript,
};

/// Performs a 'product-check', which proves that one set of polynomials equals another, invariant
/// of their sorting. Essentially, this is done by multiplying the first set together, and then the
/// second set, after which the quotient is computed. The sets equal each other iff s1 / s2 = 1.
/// Unfortunately, checking this isn't quite as simple - much like with a naive sumcheck, a
/// malicious prover can ensure that the elements are set so that the product of s1 does not equal
/// the product of s2, but becomes equal modulo the field order. Therefore, we need to perform some
/// extra trickery to make the claim sound.
pub fn prove<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    unsorted: &[MultilinearExtension<F>],
    sorted: &[MultilinearExtension<F>],
    transcript: &mut T,
) -> (
    SumcheckProof<F, E>,
    Vec<E>,
    VirtualPolynomial<F>,
    MultilinearExtension<F>,
) {
    debug_assert!(unsorted.len() == sorted.len());
    debug_assert!(unsorted
        .iter()
        .skip(1)
        .all(|p| p.num_vars() == unsorted[0].num_vars()));
    debug_assert!(sorted
        .iter()
        .skip(1)
        .all(|p| p.num_vars() == sorted[0].num_vars()));
    debug_assert!(unsorted[0].num_vars() == sorted[0].num_vars());

    // Compute the fractional polynomial.
    let mut nominator = unsorted[0].clone();
    unsorted[1..].iter().for_each(|p| {
        nominator
            .evals
            .iter_mut()
            .zip(p.evals.iter())
            .for_each(|(eval, p_eval)| {
                eval.mul_assign(p_eval);
            });
    });

    let mut denominator = sorted[0].clone();
    sorted[1..].iter().for_each(|p| {
        denominator
            .evals
            .iter_mut()
            .zip(p.evals.iter())
            .for_each(|(eval, p_eval)| {
                eval.mul_assign(p_eval);
            });
    });
    denominator.evals.iter_mut().for_each(|eval| {
        *eval = eval.inverse().unwrap();
    });

    // NOTE: We keep it as a 'MultilinearExtension' since it is supposed to be a multilinear
    // representation.
    let frac_poly = MultilinearExtension::new(
        nominator
            .evals
            .iter()
            .zip(denominator.evals.iter())
            .map(|(n_eval, d_eval)| {
                let mut n_eval = n_eval.clone();
                n_eval.mul_assign(d_eval);
                n_eval
            })
            .collect::<Vec<F>>(),
    );

    // Compute product polynomial. We want to create a polynomial f(X) in \mu + 1 vars where f(x,
    // 0) = frac(x) and f(x, 1) = f(0, x) * f(1, x).
    let mut prod_evals = Vec::with_capacity(frac_poly.len() * 2);
    frac_poly.evals.iter().for_each(|e| prod_evals.push(*e));
    for i in 0..prod_evals.len() {
        let mut lhs = prod_evals[i << 1];
        lhs.mul_assign(&prod_evals[(i << 1) + 1]);
        prod_evals.push(lhs);
    }

    let prod_poly = MultilinearExtension::new(prod_evals);

    // Now we create the zerocheck polynomial and run a zerocheck on it.
    let upper_half_prod = MultilinearExtension::new(prod_poly.evals[frac_poly.len()..].to_vec());
    let mut a = VirtualPolynomial::from(nominator.merge(upper_half_prod));

    let evens = MultilinearExtension::new(
        prod_poly
            .evals
            .iter()
            .step_by(2)
            .map(|e| *e)
            .collect::<Vec<F>>(),
    );
    let mut b = VirtualPolynomial::from(denominator.merge(evens));

    let lower_half_prod = MultilinearExtension::new(prod_poly.evals[..frac_poly.len()].to_vec());
    let odds = MultilinearExtension::new(
        prod_poly
            .evals
            .iter()
            .skip(1)
            .step_by(2)
            .map(|e| *e)
            .collect::<Vec<F>>(),
    );
    let c = VirtualPolynomial::from(lower_half_prod.merge(odds));

    b.mul_assign(&c);
    b.evals.iter_mut().for_each(|e| e.negate());
    a.add_assign(&b);

    let (proof, evals) = zerocheck::prove(&mut a, transcript);
    (proof, evals, a, prod_poly)
}

/// Verify a prodcheck proof.
pub fn verify<F: Field, E: ChallengeField<F>, T: Transcript<F>>() {}
