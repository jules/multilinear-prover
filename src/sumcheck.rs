pub struct SumcheckProof<F: Field> {
    challenges: Vec<F>,
    proofs: Vec<Vec<F>>,
}

fn do_sumcheck<F: Field, T: Transcript<F>>(
    poly: MultilinearPoly<F>,
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
    for _ in 0..n_rounds {
        let evals = do_sumcheck_round(&mut challenge, &mut poly, &mut claimed_sum, &mut transcript);
        transcript.witness_field_elements(&evals);
        proofs.push(evals);
        challenge = Some(transcript.get_challenge());
        challenges.push(challenge);
    }

    Ok(SumcheckProof { challenges, proofs })
}

fn do_sumcheck_round<F: Field, T: Transcript<F>>(
    challenge: &mut Option<F>,
    poly: &mut MultilinearPoly<F>,
    sum: &mut F,
    transcript: &mut T,
) -> Vec<F> {
    // If we have some challenge, we first need to fix the polynomial on that challenge.
    if let Some(c) = challenge {
        // TODO: check round
        poly = poly.fix_variable(c);
    }

    // Here is where we do the meat of the sumcheck.
}
