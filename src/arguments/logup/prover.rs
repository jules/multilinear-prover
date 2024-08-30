//! Wrapper for a full logup prover, combining the IOP and PCS elements.

use crate::{
    arguments::zerocheck::prover::{ZerocheckProof, ZerocheckProver},
    field::{ChallengeField, Field},
    iop::{logup, sumcheck, zerocheck},
    pcs::PolynomialCommitmentScheme,
    polynomial::{MultilinearExtension, VirtualPolynomial},
    transcript::Transcript,
};
use core::marker::PhantomData;
use rayon::prelude::*;

/// A prover for the LogUp batch-column lookup protocol. Internally, consists only of a zerocheck
/// argument prover, as this is used near the end of the proof to establish correct relations
/// between the trace columns and the table.
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
    ) -> ZerocheckProof<F, E, T, PCS> {
        debug_assert!(self.lagrange_coefficients.len() == 4);

        let multiplicities = logup::compute_multiplicities(trace_columns, table);

        self.transcript.observe_witnesses(&multiplicities.evals);
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
            .iter()
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

        let helpers = logup::compute_helper_columns(&x_plus_columns, table, &multiplicities, x);
        helpers.iter().for_each(|helper| {
            self.transcript.observe_witnesses(&helper.evals);
        });

        // XXX should probably also be extension?
        let batching_challenges = (0..x_plus_columns.len())
            .map(|_| self.transcript.draw_challenge())
            .collect::<Vec<F>>();

        let mut rho = vec![x_plus_table.clone()];
        rho.extend(x_plus_columns.clone());

        let rho_product = MultilinearExtension::new(
            x_plus_table
                .evals
                .iter()
                .enumerate()
                .map(|(i, value)| {
                    let mut value = value.clone();
                    x_plus_columns.iter().for_each(|column| {
                        value.mul_assign(&column.evals[i]);
                    });
                    value
                })
                .collect::<Vec<F>>(),
        );

        // Draw a list of challenges with which we create the lagrange kernel.
        let mut c = vec![F::ZERO; table.len()];
        c.iter_mut()
            .for_each(|e| *e = self.transcript.draw_challenge());

        // Compute lagrange kernel of the hypercube.
        // NOTE: currently just using eq, it seems to serve a similar purpose?
        let lagrange_kernel = zerocheck::compute_eq(c, table.num_vars() as u32);

        // We construct a virtual polynomial that attests to the correct lookup relation.
        let mut negated_multiplicites = multiplicities.evals.clone();
        negated_multiplicites.iter_mut().for_each(|value| {
            value.negate();
        });
        let negated_multiplicities = MultilinearExtension::new(negated_multiplicites);

        let sumchecks = helpers
            .iter()
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

                let mut product_left = VirtualPolynomial::from(helper.clone());
                product_left.mul_assign_mle(&rho[i]);

                let mut product_right = VirtualPolynomial::from(helper.clone());
                product_right.mul_assign_mle(&rho_product);

                // right hand side TODO

                product_left.add_assign(&product_right);
                product_left.mul_assign_mle(&scaled_kernel);
                product_left.add_assign_mle(helper, false);
                product_left
            })
            .collect::<Vec<VirtualPolynomial<F>>>();

        let mut sumcheck = sumchecks[0].clone();
        sumchecks.iter().skip(1).for_each(|s| {
            sumcheck.add_assign(s);
        });

        let (zerocheck_proof, eval_point) =
            sumcheck::prove(&sumcheck, &mut self.transcript, &self.lagrange_coefficients);

        let evaluations = sumcheck
            .constituents
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

        let commitment = self
            .pcs
            .commit(&sumcheck.constituents, &mut self.transcript);

        let proof = self.pcs.prove(
            &commitment,
            &sumcheck.constituents,
            &eval_point,
            &mut self.transcript,
        );

        ZerocheckProof {
            zerocheck_proof,
            evaluations,
            commitment,
            proof,
            products: sumcheck.products.clone(),
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
