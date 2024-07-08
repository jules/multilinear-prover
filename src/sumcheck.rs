pub struct SumcheckProof<F: Field, MP: MultilinearPoly<F>> {
    poly: MP,
    sum: F,
}

fn do_sumcheck<F: Field, MP: MultilinearPoly<F>>(
    poly: MP,
    claimed_sum: F,
) -> Result<SumcheckProof, SumcheckError> {
    let mut challenge: Option<F> = None;
    let n_rounds = poly.num_vars();

    let mut poly = poly.clone();
    let mut claimed_sum = claimed_sum.clone();
    for _ in 0..n_rounds {
        do_sumcheck_round(&mut challenge, &mut poly, &mut claimed_sum);
    }

    Ok(SumcheckProof {
        poly,
        sum: claimed_sum,
    })
}
