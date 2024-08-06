//! Zerocheck prover/verifier

use super::sumcheck::{self, SumcheckError, SumcheckProof};
use crate::{
    field::{ChallengeField, Field},
    mle::MultilinearExtension,
    transcript::Transcript,
};

/// Runs the zerocheck prover. This is essentially just a sumcheck, but prior to running the
/// sumcheck protocol, first multiplies the concerning polynomial with an `eq` polynomial which
/// randomly zeroes out all points except for one - allowing us to check for a constraint
/// polynomial being zero-everywhere on the cube without potential cheating by making all
/// evaluations sum to zero, for instance.
pub fn prove<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    polys: &[MultilinearExtension<F>],
    transcript: &mut T,
) -> (SumcheckProof<F, E>, Vec<E>) {
    // Draw a list of challenges with which we create the `eq` polynomial.
    let mut c = vec![F::ZERO; polys[0].num_vars()];
    c.iter_mut().for_each(|e| *e = transcript.draw_challenge());
    // For the purposes of creating `eq`, we will also create a vector of `1 - challenge`.
    let one_minus_c = c
        .iter()
        .map(|chal| {
            let mut one_minus_chal = F::ONE;
            one_minus_chal.sub_assign(chal);
            one_minus_chal
        })
        .collect::<Vec<F>>();

    // Create the `eq` polynomial. Depending on when a variable is 0 or 1, we either input (1 - r)
    // or r as the term for the multiplication. Do this for 2^num_vars.
    let mut eq_evals = vec![F::ZERO; polys[0].evals.len()];
    eq_evals.iter_mut().enumerate().for_each(|(i, eval)| {
        *eval = c.iter().enumerate().fold(F::ZERO, |mut acc, (j, chal)| {
            let flip = ((i >> j) ^ 1) != 0;
            if flip {
                acc.mul_assign(&one_minus_c[j]);
            } else {
                acc.mul_assign(&chal);
            }
            acc
        });
    });

    // Perform the polynomial product. This is just entrywise mult.
    // XXX i think this is something that the sumcheck prover actually shouldn't do, and instead we
    // should send some abstraction like a 'CompositionPolynomial' which contains degree
    // information.
    let mut poly = polys[0].clone();
    for p in &polys[1..] {
        poly.evals
            .iter_mut()
            .zip(p.evals.iter())
            .for_each(|(eval, m)| eval.mul_assign(m));
    }

    // Now multiply with our constraint poly.
    poly.evals
        .iter_mut()
        .zip(eq_evals.iter())
        .for_each(|(eval, eq)| {
            eval.mul_assign(eq);
        });

    // Finally, just run the sumcheck prover and return the needed information.
    sumcheck::prove(&[poly], transcript)
}

/// Runs the zerocheck verifier.
pub fn verify<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    proof: SumcheckProof<F, E>,
    eval: Vec<E>,
    transcript: &mut T,
) -> Result<E, SumcheckError> {
    // Draw a list of challenges with which we create the `eq` polynomial.
    let mut c = vec![F::ZERO; eval.len()];
    c.iter_mut().for_each(|e| *e = transcript.draw_challenge());

    // Get a claimed sum for this proof.
    let mut claim = sumcheck::verify(proof, transcript)?;

    // Evaluate `eq` at `eval`.
    let mut res = E::ONE;
    c.iter().zip(eval.iter()).for_each(|(c, e)| {
        let mut ce = e.clone();
        ce.mul_base(c);
        ce.add_assign(&ce.clone());
        ce.sub_base(c);
        ce.sub_assign(e);
        ce.add_assign(&E::ONE);
        res.mul_assign(&ce);
    });

    // Output claim / res
    let res_inverse = res.inverse().expect("res should have an inverse");
    claim.mul_assign(&res_inverse);
    Ok(claim)
}
