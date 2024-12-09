//! Zerocheck prover/verifier.

use super::sumcheck::{self, SumcheckError, SumcheckProof};
use crate::{
    field::{ChallengeField, Field},
    polynomial::{MultilinearExtension, VirtualPolynomial},
    transcript::Transcript,
};

/// Runs the zerocheck prover. This is essentially just a sumcheck, but prior to running the
/// sumcheck protocol, first multiplies the concerning polynomial with an `eq` polynomial which
/// randomly zeroes out all points except for one - allowing us to check for a constraint
/// polynomial being zero-everywhere on the cube without potential cheating by making all
/// evaluations sum to zero, for instance.
pub fn prove<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    poly: &mut VirtualPolynomial<F>,
    transcript: &mut T,
    precomputed: &[Vec<F>],
) -> (SumcheckProof<F, E>, Vec<E>) {
    // Draw a list of challenges with which we create the `eq` polynomial.
    let mut c = vec![F::ZERO; poly.num_vars()];
    c.iter_mut().for_each(|e| *e = transcript.draw_challenge());

    // Multiply with our polynomial.
    poly.mul_assign_mle(&compute_eq(c, poly.num_vars() as u32));

    // Finally, just run the sumcheck prover and return the needed information.
    sumcheck::prove(&poly, transcript, precomputed)
}

/// Runs the zerocheck verifier.
pub fn verify<F: Field, E: ChallengeField<F>, T: Transcript<F>>(
    proof: &SumcheckProof<F, E>,
    transcript: &mut T,
) -> Result<((E, Vec<E>), MultilinearExtension<F>), SumcheckError> {
    // Draw a list of challenges with which we create the `eq` polynomial.
    let mut c = vec![F::ZERO; proof.proofs.len()];
    c.iter_mut().for_each(|e| *e = transcript.draw_challenge());

    // TODO: eq eval and not just recreating it entirely.
    // Return a claimed sum for this proof.
    Ok((
        sumcheck::verify(proof, transcript)?,
        compute_eq(c, proof.proofs.len() as u32),
    ))
}

pub fn compute_eq<F: Field>(challenges: Vec<F>, num_vars: u32) -> MultilinearExtension<F> {
    // For the purposes of creating `eq`, we will also create a vector of `1 - challenge`.
    let one_minus_c = challenges
        .iter()
        .map(|chal| {
            let mut one_minus_chal = F::ONE;
            one_minus_chal.sub_assign(chal);
            one_minus_chal
        })
        .collect::<Vec<F>>();

    // Create the `eq` polynomial. Depending on when a variable is 0 or 1, we either input (1 - r)
    // or r as the term for the multiplication.
    let mut eq_evals = vec![F::ZERO; 2u32.pow(num_vars) as usize];
    eq_evals.iter_mut().enumerate().for_each(|(i, eval)| {
        *eval = challenges
            .iter()
            .enumerate()
            .fold(F::ONE, |mut acc, (j, chal)| {
                let flip = ((i >> j) & 1) == 0;
                if flip {
                    acc.mul_assign(&one_minus_c[j]);
                } else {
                    acc.mul_assign(&chal);
                }
                acc
            });
    });

    MultilinearExtension::new(eq_evals)
}
