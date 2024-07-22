use crate::{field::Field, mle::MultilinearExtension, transcript::Transcript};

pub enum SumcheckError {
    ExpectedChallenge,
    ChallengeOnFirstRound,
}

pub struct SumcheckProof<F: Field> {
    proofs: Vec<Vec<F>>,
    claimed_sums: Vec<F>,
}

pub fn prove<F: Field, T: Transcript<F>>(
    poly: &MultilinearExtension<F>,
    transcript: &mut T,
) -> Result<SumcheckProof<F>, SumcheckError> {
    let n_rounds = poly.num_vars();
    let mut proofs = Vec::with_capacity(n_rounds);
    let mut claimed_sums = Vec::with_capacity(n_rounds);
    let mut poly = poly.clone();

    for i in 0..(n_rounds - 1) {
        let (coeffs, sum) = sumcheck_step(&poly, transcript);
        transcript.observe_witnesses(&coeffs);
        let challenge = transcript.draw_challenge();

        proofs.push(coeffs);
        claimed_sums.push(sum);
        poly.fix_variable(challenge);
    }

    // at the end, we only need to figure out the last poly and sum
    let (coeffs, sum) = sumcheck_step(&poly, transcript);
    proofs.push(coeffs);
    claimed_sums.push(sum);

    Ok(SumcheckProof {
        proofs,
        claimed_sums,
    })
}

#[inline(always)]
fn sumcheck_step<F: Field, T: Transcript<F>>(
    poly: &MultilinearExtension<F>,
    transcript: &mut T,
) -> (Vec<F>, F) {
    let evals = poly.sum_evaluations();
    let coeffs = lagrange_interpolation(&evals);
    (coeffs, evals.iter().fold(F::ZERO, |acc, x| acc + x))
}

pub fn verify<F: Field, T: Transcript<F>>(
    proof: SumcheckProof<F>,
    transcript: &mut T,
) -> Result<bool, SumcheckError> {
    let SumcheckProof {
        proofs,
        claimed_sums,
    } = proof;

    let mut claim = None;
    let mut challenges = Vec::with_capacity(proofs.len());
    for (coeffs, sum) in proofs.iter().zip(claimed_sums) {
        if let Some(c) = claim {
            if c != sum {
                return Ok(false);
            }
        }

        transcript.observe_witnesses(&coeffs);
        let c = transcript.draw_challenge();
        challenges.push(c);

        let res = univariate_eval(&coeffs, F::ZERO) + univariate_eval(&coeffs, F::ONE);
        if res != sum {
            return Ok(false);
        }

        claim = Some(univariate_eval(&coeffs, c));
    }

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
