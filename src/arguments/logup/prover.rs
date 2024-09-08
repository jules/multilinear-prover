//! Wrapper for a full logup prover, combining the IOP and PCS elements.

use crate::{
    field::{ChallengeField, Field},
    iop::{
        logup,
        sumcheck::{self, SumcheckProof},
        zerocheck,
    },
    pcs::PolynomialCommitmentScheme,
    polynomial::{MultilinearExtension, VirtualPolynomial},
    transcript::{IntoObservable, Transcript},
};
use core::marker::PhantomData;
use rayon::prelude::*;

/// A prover for the LogUp batch-column lookup protocol.
pub struct LogUpProver<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    lagrange_coefficients: Vec<Vec<F>>,
    pub transcript: T,
    pub pcs: PCS,
    _e_marker: PhantomData<E>,
}

pub struct LogUpProof<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    pub zerocheck_proof: SumcheckProof<F, E>,
    pub multiplicities_commitment: PCS::Commitment,
    pub multiplicities_proof: PCS::Proof,
    pub helpers_commitment: PCS::Commitment,
    pub helpers_proof: PCS::Proof,
    pub sumcheck_commitment: PCS::Commitment,
    pub sumcheck_proof: PCS::Proof,
    pub evaluations: Vec<E>,
    pub num_helpers: usize,
    pub num_vars: usize,
}

impl<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, E, T>,
    > LogUpProver<F, E, T, PCS>
{
    pub fn new(transcript: T, pcs: PCS, lagrange_coefficients: Vec<Vec<F>>) -> Self {
        Self {
            lagrange_coefficients,
            transcript,
            pcs,
            _e_marker: PhantomData::<E>,
        }
    }

    /// Proves the lookup relation with support for batch-column lookups which is particularly
    /// efficient in the VM scenario where many columns are, for instance, subject to one and the
    /// same range check.
    pub fn prove(
        &mut self,
        trace_columns: &[MultilinearExtension<F>],
        table: &MultilinearExtension<F>,
        multiplicities: &MultilinearExtension<F>,
    ) -> LogUpProof<F, E, T, PCS> {
        debug_assert!(self.lagrange_coefficients.len() == 4);

        let multiplicities_commitment = self
            .pcs
            .commit(&[multiplicities.clone()], &mut self.transcript);
        self.transcript
            .observe_hashes(&multiplicities_commitment.into_observable());

        // XXX do we lift into extension here?
        let x = self.transcript.draw_challenge();

        let x_plus_table = MultilinearExtension::new(
            table
                .evals
                .iter()
                .map(|value| {
                    let mut value = value.clone();
                    value.add_assign(&x);
                    value
                })
                .collect::<Vec<F>>(),
        );

        // Compute x + column for all trace columns.
        let x_plus_columns = trace_columns
            .par_iter()
            .map(|column| {
                MultilinearExtension::new(
                    column
                        .evals
                        .iter()
                        .map(|value| {
                            let mut value = value.clone();
                            value.add_assign(&x);
                            value
                        })
                        .collect::<Vec<F>>(),
                )
            })
            .collect::<Vec<MultilinearExtension<F>>>();

        let helpers =
            logup::compute_helper_columns(&x_plus_columns, &x_plus_table, &multiplicities);

        // Commit and observe.
        let helpers_commitment = self.pcs.commit(&helpers, &mut self.transcript);
        self.transcript
            .observe_hashes(&helpers_commitment.into_observable());

        let batching_challenges = (0..helpers.len())
            .map(|_| self.transcript.draw_challenge())
            .collect::<Vec<F>>();

        // Draw a list of challenges with which we create the lagrange kernel.
        let mut c = vec![F::ZERO; table.num_vars()];
        c.iter_mut()
            .for_each(|e| *e = self.transcript.draw_challenge());

        // Compute lagrange kernel of the hypercube.
        let lagrange_kernel = {
            let mut evals = vec![F::ZERO; table.len()];
            let one_over_2_to_n = F::from_usize(table.len()).inverse().unwrap();

            evals.iter_mut().enumerate().for_each(|(i, eval)| {
                let mut i = i.clone();
                let mut result = F::ONE;
                (0..table.num_vars()).for_each(|j| {
                    if i & 1 == 1 {
                        let mut one_plus_c = c[j].clone();
                        one_plus_c.add_assign(&F::ONE);
                        result.mul_assign(&one_plus_c);
                    }

                    i >>= 1;
                });

                result.mul_assign(&one_over_2_to_n);
                *eval = result;
            });

            MultilinearExtension::new(evals)
        };

        // We construct a virtual polynomial that attests to the correct lookup relation.
        let mut rho = vec![x_plus_table.clone()];
        rho.extend(x_plus_columns.clone());

        let mut rhs = Vec::with_capacity(x_plus_columns.len() + 1);
        rhs.push(multiplicities.clone());
        let mut negative_one = F::ONE;
        negative_one.negate();
        rhs.extend(
            (0..x_plus_columns.len())
                .map(|_| MultilinearExtension::new(vec![negative_one; table.len()]))
                .collect::<Vec<MultilinearExtension<F>>>(),
        );

        let sumchecks = helpers
            .par_iter()
            .enumerate()
            .map(|(i, helper)| {
                let scaled_kernel = MultilinearExtension::new(
                    lagrange_kernel
                        .evals
                        .iter()
                        .map(|v| {
                            let mut v = v.clone();
                            v.mul_assign(&batching_challenges[i]);
                            v
                        })
                        .collect::<Vec<F>>(),
                );

                let mut sumcheck = VirtualPolynomial::from(helper.clone());
                sumcheck.mul_assign_mle(&rho[i]);
                sumcheck.add_assign_mle(&rhs[i], true);
                sumcheck.mul_assign_mle(&scaled_kernel);
                sumcheck.add_assign_mle(helper, false);
                sumcheck
            })
            .collect::<Vec<VirtualPolynomial<F>>>();

        let mut sumcheck = sumchecks[0].clone();
        sumchecks.iter().skip(1).for_each(|s| {
            sumcheck.add_assign(s);
        });

        let (zerocheck_proof, eval_point) =
            sumcheck::prove(&sumcheck, &mut self.transcript, &self.lagrange_coefficients);

        let multiplicities_evaluation = {
            let mut m_lifted = multiplicities.fix_variable_ext(eval_point[0]);
            (1..eval_point.len()).for_each(|i| {
                m_lifted.fix_variable(eval_point[i]);
            });

            m_lifted.evals[0]
        };

        let helpers_evaluations = helpers
            .par_iter()
            .map(|p| {
                let mut p_lifted = p.fix_variable_ext(eval_point[0]);
                (1..eval_point.len()).for_each(|i| {
                    p_lifted.fix_variable(eval_point[i]);
                });

                debug_assert!(p_lifted.evals.len() == 1);
                p_lifted.evals[0]
            })
            .collect::<Vec<E>>();

        let rho_evaluations = rho
            .par_iter()
            .map(|p| {
                let mut p_lifted = p.fix_variable_ext(eval_point[0]);
                (1..eval_point.len()).for_each(|i| {
                    p_lifted.fix_variable(eval_point[i]);
                });

                debug_assert!(p_lifted.evals.len() == 1);
                p_lifted.evals[0]
            })
            .collect::<Vec<E>>();

        let mut evaluations = vec![multiplicities_evaluation];
        evaluations.extend(helpers_evaluations);
        evaluations.extend(rho_evaluations);

        let sumcheck_commitment = self.pcs.commit(&rho, &mut self.transcript);
        self.transcript
            .observe_hashes(&sumcheck_commitment.into_observable());

        let multiplicities_proof = self.pcs.prove(
            &multiplicities_commitment,
            &[multiplicities.clone()],
            &eval_point,
            &mut self.transcript,
        );

        let helpers_proof = self.pcs.prove(
            &helpers_commitment,
            &helpers,
            &eval_point,
            &mut self.transcript,
        );

        let sumcheck_proof = self.pcs.prove(
            &sumcheck_commitment,
            &rho,
            &eval_point,
            &mut self.transcript,
        );

        LogUpProof {
            zerocheck_proof,
            multiplicities_commitment,
            multiplicities_proof,
            helpers_commitment,
            helpers_proof,
            sumcheck_commitment,
            sumcheck_proof,
            evaluations,
            num_helpers: helpers.len(),
            num_vars: table.num_vars(),
        }
    }

    /// We implement a method to relinquish the transcript out of the prover after it is finished.
    /// This method will consume the prover and leave the transcript to be used for the next
    /// argument.
    pub fn relinquish_transcript(mut self) -> T {
        let transcript = core::mem::take(&mut self.transcript);
        drop(self);
        transcript
    }
}
