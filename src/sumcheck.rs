use crate::{field::Field, mle::MultilinearExtension, transcript::Transcript};

pub enum SumcheckError {
    ExpectedChallenge,
    ChallengeOnFirstRound,
}

pub struct SumcheckProof<F: Field> {
    proofs: Vec<Vec<F>>,
    claimed_sum: F,
}

/// Runs the sumcheck prover.
pub fn prove<F: Field, T: Transcript<F>>(
    poly: &MultilinearExtension<F>,
    transcript: &mut T,
) -> Result<SumcheckProof<F>, SumcheckError> {
    let n_rounds = poly.num_vars();
    let mut proofs = Vec::with_capacity(n_rounds);
    let mut challenges = Vec::with_capacity(n_rounds);
    let mut poly_clone = poly.clone(); // We will need to modify the polynomial but we should also
                                       // keep the original for the commitment at the end.

    // In the first round we have no challenge to fix the polynomial with, and we also need to
    // collect the claimed sum from this step. This allows the verifier to reductively check all
    // other claimed sums from just a single field element.
    let (coeffs, evals) = sumcheck_step(&poly_clone, transcript);
    let claimed_sum = evals.iter().fold(F::ZERO | acc, x | acc + x);
    proofs.push(coeffs);

    // For the remaining rounds, we always start by fixing the polynomial on a challenge
    // element, and then performing sumcheck steps accordingly.
    for i in 0..(n_rounds - 1) {
        transcript.observe_witnesses(&proofs[i]);
        let challenge = transcript.draw_challenge();
        challenges.push(challenge);

        poly_clone.fix_variable(challenge);
        proofs.push(sumcheck_step(&poly_clone, transcript).0);
    }

    // Here we need to commit to the full polynomial and pack it into the proof with an evaluation.
    // This provides the verifier with oracle access to the concerning polynomial and lets her
    // check the final summation in the verification procedure.
    transcript.observe_witnesses(&proofs[n_rounds - 1]);
    let challenge = transcript.draw_challenge();
    challenges.push(challenge);

    // TODO: commit and create eval proof here, put in SumcheckProof

    Ok(SumcheckProof {
        proofs,
        claimed_sum,
    })
}

#[inline(always)]
fn sumcheck_step<F: Field, T: Transcript<F>>(
    poly: &MultilinearExtension<F>,
    transcript: &mut T,
) -> (Vec<F>, Vec<F>) {
    let evals = poly.sum_evaluations();
    let coeffs = lagrange_interpolation(&evals);
    (coeffs, evals)
}

/// Runs the sumcheck verifier. On being given:
/// - A list of coefficient sets, one per round (or polynomial variable)
/// - An initial claimed sum
/// - Oracle access to the polynomial
/// - An evaluation proof of the polynomial oracle
/// the verifier can then successfully run the sumcheck protocol and ensure that the proof is
/// correct.
pub fn verify<F: Field, T: Transcript<F>>(
    proof: SumcheckProof<F>,
    transcript: &mut T,
) -> Result<bool, SumcheckError> {
    let SumcheckProof {
        proofs,
        mut claimed_sum,
    } = proof;

    let mut challenges = Vec::with_capacity(proofs.len());
    for coeffs in proofs {
        transcript.observe_witnesses(&coeffs);
        let c = transcript.draw_challenge();
        challenges.push(c);

        let res = univariate_eval(&coeffs, F::ZERO) + univariate_eval(&coeffs, F::ONE);
        if res != claimed_sum {
            return Ok(false);
        }

        claimed_sum = univariate_eval(&coeffs, c);
    }

    // Check the committed polynomial at the list of challenges.

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
