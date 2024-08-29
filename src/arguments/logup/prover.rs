//! Wrapper for a full logup prover, combining the IOP and PCS elements.

use crate::{
    arguments::zerocheck::prover::{ZerocheckProof, ZerocheckProver},
    field::{ChallengeField, Field},
    iop::logup,
    pcs::PolynomialCommitmentScheme,
    polynomial::MultilinearExtension,
    transcript::Transcript,
};

/// A prover for the LogUp batch-column lookup protocol. Internally, consists only of a zerocheck
/// argument prover, as this is used near the end of the proof to establish correct relations
/// between the trace columns and the table.
pub struct LogUpProver<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    zerocheck_prover: ZerocheckProver<F, E, T, PCS>,
}

impl<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, E, T>,
    > LogUpProver<F, E, T, PCS>
{
    /// Creates a new LogUp prover by feeding the given arguments into the creation of a
    /// zerocheck prover.
    pub fn new(transcript: T, pcs: PCS, lagrange_coefficients: Vec<Vec<F>>) -> Self {
        Self {
            zerocheck_prover: ZerocheckProver::new(transcript, pcs, lagrange_coefficients),
        }
    }

    /// Proves the lookup relation with support for batch-column lookups which is particularly
    /// efficient in the VM scenario where many columns are, for instance, subject to one and the
    /// same range check.
    pub fn prove(
        &mut self,
        trace_columns: &[MultilinearExtension<F>],
        table: &MultilinearExtension<F>,
    ) -> ZerocheckProof<F, E, T, PCS> {
        let multiplicities = logup::compute_multiplicities(trace_columns, table);

        self.zerocheck_prover
            .transcript
            .observe_witnesses(&multiplicities.evals);
        // XXX do we lift into extension here?
        let x = self.zerocheck_prover.transcript.draw_challenge();
        let helpers = logup::compute_helper_columns(trace_columns, table, &multiplicities, x);
        for helper in helpers {
            self.zerocheck_prover
                .transcript
                .observe_witnesses(&helper.evals);
        }

        // TODO zerocheck here
        todo!()
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
