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
    let mut challenge: Option<F> = None;
    let n_rounds = poly.num_vars();
    let mut proofs = Vec::with_capacity(n_rounds);

    let mut poly = poly.clone();
    let mut claimed_sum = None;
    for round in 0..n_rounds {
        let evals = sum_polynomial(&mut challenge, &mut poly, round + 1)?;
        if round == 0 {
            // Reuse the work to set our claimed sum.
            claimed_sum = Some(evals.clone().into_iter().fold(F::ZERO, |acc, x| acc + x));
        }
        let coeffs = lagrange_interpolation(evals);

        transcript.observe_witnesses(&coeffs);
        proofs.push(coeffs);
        challenge = Some(transcript.draw_challenge());
    }

    Ok(SumcheckProof {
        proofs,
        claimed_sum: claimed_sum.unwrap(),
    })
}

fn sum_polynomial<F: Field>(
    challenge: &mut Option<F>,
    poly: &mut MultilinearExtension<F>,
    round: usize,
) -> Result<Vec<F>, SumcheckError> {
    // If we have some challenge, we first need to fix the polynomial on that challenge.
    if let Some(c) = challenge {
        if round == 1 {
            return Err(SumcheckError::ChallengeOnFirstRound);
        }

        *poly = poly.fix_variable(*c);
    } else {
        if round > 1 {
            return Err(SumcheckError::ExpectedChallenge);
        }
    }

    // Here is where we get into the meat of the sumcheck.
    // We want to sum up all the evaluations of the polynomial `f` except for those at the variable
    // at `round`. We assume lexicographic ordering, so poly[0] and poly[1] will be f(0, 0, ..., 0)
    // and f(1, 0, ..., 0).
    let mut evals_0 = F::ZERO;
    let mut evals_1 = F::ZERO;
    for i in 0..poly.num_vars() - 1 {
        evals_0 += poly.eval(i << 1);
        evals_1 += poly.eval((i << 1) + 1);
    }

    Ok(vec![evals_0, evals_1])
}

pub fn verify<F: Field, T: Transcript<F>>(
    proof: SumcheckProof<F>,
    transcript: T,
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

        let res = univariate_eval(p, 0) + univariate_eval(p, 1);
        if res != claimed_sum {
            return Ok(false);
        }

        claimed_sum = univariate_eval(p, c);
    }

    // check final eval on committed poly

    Ok(true)
}

// Standard lagrange interpolation, assuming indices for evals are 0, 1, 2, ...
// NOTE: can be sped up if we precompute lagrange coeffs for a given size
fn lagrange_interpolation<F: Field>(evals: Vec<F>) -> Vec<F> {
    evals
        .iter()
        .enumerate()
        .map(|(i, eval)| {
            let denom = {
                (0..evals.len())
                    .map(|j| {
                        let i = F::from_usize(i);
                        let j = F::from_usize(j);
                        if i == j {
                            F::ONE
                        } else {
                            i - j
                        }
                    })
                    .fold(F::ONE, |acc, x| acc * x)
            };

            let num = F::ONE; // TODO
            let denom_inv = denom
                .inverse()
                .expect("lagrange coefficient denominator can not be zero");

            *eval * num * denom_inv
        })
        .collect::<Vec<F>>()
}
