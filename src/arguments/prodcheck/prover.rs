//! A prover for the prodcheck argument.

use crate::{
    arguments::zerocheck::prover::{ZerocheckProof, ZerocheckProver},
    field::{ChallengeField, Field},
    iop::prodcheck,
    pcs::PolynomialCommitmentScheme,
    polynomial::MultilinearExtension,
    transcript::Transcript,
};

pub struct ProdcheckProver<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    zerocheck_prover: ZerocheckProver<F, E, T, PCS>,
}

pub struct ProdcheckProof<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    pub zerocheck_proof: ZerocheckProof<F, E, T, PCS>,
    pub commitment: PCS::Commitment,
    pub proof: PCS::Proof,
    pub num_vars: usize, // actually, this can come from a verifier key that gives circuit length
}

impl<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, E, T>,
    > ProdcheckProver<F, E, T, PCS>
{
    /// Creates a new prodcheck prover by feeding the given arguments into the creation of a
    /// zerocheck prover.
    pub fn new(transcript: T, pcs: PCS, lagrange_coefficients: Vec<Vec<F>>) -> Self {
        Self {
            zerocheck_prover: ZerocheckProver::new(transcript, pcs, lagrange_coefficients),
        }
    }

    pub fn prove(
        &mut self,
        unsorted_columns: &[MultilinearExtension<F>],
        sorted_columns: &[MultilinearExtension<F>],
    ) -> ProdcheckProof<F, E, T, PCS> {
        let (frac_poly, nominator, denominator) =
            prodcheck::compute_frac_poly(unsorted_columns, sorted_columns);
        let mut v = prodcheck::compute_v(frac_poly);
        let h = prodcheck::compute_h(&v, nominator, denominator);
        let zerocheck_proof = self.zerocheck_prover.prove(h);

        // We also need a commitment and evaluation proof of v(0, 1, ..., 1).
        let mut eval_point = vec![E::ONE; v.num_vars()];
        eval_point[0] = E::ZERO;

        let num_vars = v.num_vars();
        let v = [v];
        let commitment = self
            .zerocheck_prover
            .pcs
            .commit(&v, &mut self.zerocheck_prover.transcript);

        let proof = self.zerocheck_prover.pcs.prove(
            &commitment,
            &v,
            &eval_point,
            &mut self.zerocheck_prover.transcript,
        );

        ProdcheckProof {
            zerocheck_proof,
            commitment,
            proof,
            num_vars,
        }
    }

    /// We implement a method to relinquish the transcript out of the prover after it is finished.
    /// This method will consume the prover and leave the transcript to be used for the next
    /// argument.
    pub fn relinquish_transcript(mut self) -> T {
        let transcript = core::mem::take(&mut self.zerocheck_prover.transcript);
        drop(self);
        transcript
    }
}
