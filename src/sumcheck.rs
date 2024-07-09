use crate::{field::Field, mle::MultilinearExtension, transcript::Transcript};

pub enum SumcheckError {
    ExpectedChallenge,
    ChallengeOnFirstRound,
}

pub struct SumcheckProof<F: Field> {
    challenges: Vec<F>,
    proofs: Vec<Vec<F>>,
}

pub fn prove<F: Field, T: Transcript<F>>(
    poly: MultilinearExtension<F>,
    claimed_sum: F,
    transcript: &mut T,
) -> Result<SumcheckProof, SumcheckError> {
    let mut challenge: Option<F> = None;
    let n_rounds = poly.num_vars();
    transcript.witness_field_element(F::from_u64_unchecked(n_rounds));
    let mut challenges = Vec::with_capacity(n_rounds);
    let mut proofs = Vec::with_capacity(n_rounds);

    let mut poly = poly.clone();
    let mut claimed_sum = claimed_sum.clone();
    for round in 0..n_rounds {
        let evals = do_sumcheck_round(&mut challenge, &mut poly, &mut claimed_sum, round + 1)?;
        transcript.observe_witnesses(&evals);
        proofs.push(evals);
        challenge = Some(transcript.draw_challenge());
        challenges.push(challenge);
    }

    Ok(SumcheckProof { challenges, proofs })
}

fn do_sumcheck_round<F: Field>(
    challenge: &mut Option<F>,
    poly: &mut MultilinearExtension<F>,
    sum: &mut F,
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

    // Here is where we do the meat of the sumcheck.
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
    proof: SumcheckProof,
    transcript: T,
) -> Result<bool, SumcheckError> {
    // interpolate the lines on each sumcheck round and check low degree
    Ok(false)
}
