//! Product check IOP utils. This module does not include any clear prover or verifier since this
//! should be taken from [`zerocheck`]. For a full prodcheck argument, check out
//! [`crate::arguments::prodcheck`].

use crate::{
    field::Field,
    polynomial::{MultilinearExtension, VirtualPolynomial},
};

/// Computes the fractional polynomial, declared as:
/// (unsorted_0 * unsorted_1 * ... * unsorted_n) / (sorted_0 * sorted_1 * .. * sorted_n)
pub fn compute_frac_poly<F: Field>(
    unsorted: &[MultilinearExtension<F>],
    sorted: &[MultilinearExtension<F>],
) -> (
    MultilinearExtension<F>,
    MultilinearExtension<F>,
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

    // NOTE: We keep it as a 'MultilinearExtension' since it is supposed to be a multilinear
    // representation.
    let frac_poly = MultilinearExtension::new(
        nominator
            .evals
            .iter()
            .zip(denominator.evals.iter())
            .map(|(n_eval, d_eval)| {
                let mut n_eval = n_eval.clone();
                n_eval.mul_assign(&d_eval.inverse().unwrap());
                n_eval
            })
            .collect::<Vec<F>>(),
    );

    (frac_poly, nominator, denominator)
}

/// Computes the 'product polynomial'. We want to create a polynomial f(X) in \mu + 1 vars where f(x,
/// 0) = frac(x) and f(x, 1) = f(0, x) * f(1, x).
pub fn compute_v<F: Field>(frac_poly: MultilinearExtension<F>) -> MultilinearExtension<F> {
    let mut prod_evals = Vec::with_capacity(frac_poly.len() * 2);
    frac_poly.evals.iter().for_each(|e| prod_evals.push(*e));
    for i in 0..(prod_evals.len() - 1) {
        let mut lhs = prod_evals[i << 1];
        lhs.mul_assign(&prod_evals[(i << 1) + 1]);
        prod_evals.push(lhs);
    }

    // The final product should equal one.
    debug_assert!(*prod_evals.last().unwrap() == F::ONE);

    // Last evaluation should be zero.
    prod_evals.push(F::ZERO);

    MultilinearExtension::new(prod_evals)
}

/// Computes the zerocheck polynomial for the prodcheck.
pub fn compute_h<F: Field>(
    prod_poly: &MultilinearExtension<F>,
    nominator: MultilinearExtension<F>,
    denominator: MultilinearExtension<F>,
) -> VirtualPolynomial<F> {
    // Now we create the zerocheck polynomial and run a zerocheck on it.
    let upper_half_prod =
        MultilinearExtension::new(prod_poly.evals[prod_poly.len() / 2..].to_vec());
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

    let lower_half_prod =
        MultilinearExtension::new(prod_poly.evals[..prod_poly.len() / 2].to_vec());
    let odds = MultilinearExtension::new(
        prod_poly
            .evals
            .iter()
            .skip(1)
            .step_by(2)
            .map(|e| *e)
            .collect::<Vec<F>>(),
    );
    let c = lower_half_prod.merge(odds);

    // a - b * c
    b.mul_assign_mle(&c);
    b.negate_product(0);
    a.add_assign(&b);
    a
}
