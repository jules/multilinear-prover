use crate::{field::Field, mle::MultilinearExtension, transcript::Transcript};

pub enum SumcheckError {
    ExpectedChallenge,
    ChallengeOnFirstRound,
}

pub struct SumcheckProof<F: Field> {
    proofs: Vec<Vec<F>>,
    claimed_sum: F,
}

pub fn prove<F: Field, T: Transcript<F>>(
    poly: MultilinearExtension<F>,
    transcript: &mut T,
) -> Result<SumcheckProof<F>, SumcheckError> {
    let n_rounds = poly.num_vars();
    let mut proofs = Vec::with_capacity(n_rounds);

    let mut poly = poly.clone();
    // first round
    let (evals, mut challenge) = sumcheck_step(&poly, transcript, &mut proofs);
    // Reuse the work to set our claimed sum.
    let claimed_sum = evals.clone().into_iter().fold(F::ZERO, |acc, x| acc + x);

    // the rest
    for _ in 0..(n_rounds - 1) {
        poly.fix_variable(challenge);
        (_, challenge) = sumcheck_step(&poly, transcript, &mut proofs);
    }

    // NOTE: need to do the last summation here right?

    Ok(SumcheckProof {
        proofs,
        claimed_sum,
    })
}

fn sumcheck_step<F: Field, T: Transcript<F>>(
    poly: &MultilinearExtension<F>,
    transcript: &mut T,
    proofs: &mut Vec<Vec<F>>,
) -> (Vec<F>, F) {
    let evals = poly.sum_evaluations();
    let coeffs = lagrange_interpolation(&evals);

    transcript.observe_witnesses(&coeffs);
    proofs.push(coeffs);
    (evals, transcript.draw_challenge())
}

pub fn verify<F: Field, T: Transcript<F>>(
    proof: SumcheckProof<F>,
    transcript: &mut T,
) -> Result<bool, SumcheckError> {
    let SumcheckProof {
        proofs,
        mut claimed_sum,
    } = proof;

    let mut challenges = Vec::with_capacity(proofs.len());
    for p in proofs {
        transcript.observe_witnesses(&p);
        let c = transcript.draw_challenge();
        challenges.push(c);

        let res = univariate_eval(&p, F::ZERO) + univariate_eval(&p, F::ONE);
        if res != claimed_sum {
            return Ok(false);
        }

        claimed_sum = univariate_eval(&p, c);
    }

    // check final eval on committed poly

    Ok(true)
}

// Standard lagrange interpolation, assuming indices for evals are 0, 1, 2, ...
// NOTE: can be sped up if we precompute lagrange coeffs for a given size
fn lagrange_interpolation<F: Field>(evals: &[F]) -> Vec<F> {
    let mut polynomial = vec![F::ZERO; evals.len()];
    evals.iter().enumerate().for_each(|(i, eval)| {
        let denom = {
            (0..evals.len()).fold(F::ONE, |acc, j| {
                let i = F::from_usize(i);
                let j = F::from_usize(j);
                if i != j {
                    acc * (i - j)
                } else {
                    acc
                }
            })
        };

        let denom_inv = denom
            .inverse()
            .expect("lagrange coefficient denominator can not be zero");
        let coeffs = {
            let mut coeffs = vec![F::ZERO; evals.len()];
            coeffs[0] = denom_inv;

            for k in 0..evals.len() {
                if k == i {
                    continue;
                }

                let start = if k < i { k + 1 } else { k };
                let mut new_coeffs = vec![F::ZERO; evals.len()];
                for j in (0..start).rev() {
                    new_coeffs[j + 1] += coeffs[j];
                    new_coeffs[j] -= F::from_usize(i) * coeffs[j];
                }
                coeffs = new_coeffs;
            }

            coeffs
        };

        polynomial
            .iter_mut()
            .zip(coeffs.iter())
            .for_each(|(v, c)| *v += *eval * c);
    });

    polynomial
}

// Simple Horner's method evaluation.
fn univariate_eval<F: Field>(coeffs: &[F], point: F) -> F {
    coeffs
        .iter()
        .enumerate()
        .rev()
        .fold(F::ZERO, |acc, (i, coeff)| {
            (acc + coeff) * point.pow(i as u32)
        })
}
