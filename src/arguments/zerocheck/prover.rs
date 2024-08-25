//! Wrapper for a full zerocheck prover, combining the IOP and PCS elements.

use crate::{
    field::{ChallengeField, Field},
    iop::{sumcheck::SumcheckProof, zerocheck},
    pcs::PolynomialCommitmentScheme,
    polynomial::VirtualPolynomial,
    transcript::Transcript,
};
use core::marker::PhantomData;
use rayon::prelude::*;

/// A zerocheck prover, which, on being given a virtual polynomial, constructs a zerocheck proof
/// and a commitment of the contents. This prover specifically convinces a verifier that, for a
/// given trace and a constraint polynomial, each row of the trace evaluates zero-everywhere on the
/// hypercube. In other words, this prover argues that all simple constraints (excluding memory
/// arguments and lookup arguments) are correctly applied.
pub struct ZerocheckProver<
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

/// A full zerocheck proof, which includes the oracle messages for the zerocheck IOP, as well as
/// the commitment and evaluation proof to all the constituents of the virtual polynomial for which
/// a zerocheck was proven.
pub struct ZerocheckProof<
    F: Field,
    E: ChallengeField<F>,
    T: Transcript<F>,
    PCS: PolynomialCommitmentScheme<F, E, T>,
> {
    pub zerocheck_proof: SumcheckProof<F, E>,
    pub evaluations: Vec<E>,
    pub commitment: PCS::Commitment,
    pub proof: PCS::Proof,
    pub products: Vec<(Vec<usize>, usize, bool)>,
}

impl<
        F: Field,
        E: ChallengeField<F>,
        T: Transcript<F>,
        PCS: PolynomialCommitmentScheme<F, E, T>,
    > ZerocheckProver<F, E, T, PCS>
{
    /// Creates a new zerocheck prover with the given transcript, polynomial commitment scheme and
    /// a set of lagrange coefficients used in the zerocheck IOP.
    pub fn new(transcript: T, pcs: PCS, lagrange_coefficients: Vec<Vec<F>>) -> Self {
        Self {
            lagrange_coefficients,
            transcript,
            pcs,
            _e_marker: PhantomData::<E>,
        }
    }

    /// Create a zerocheck proof for the given [`VirtualPolynomial`].
    pub fn prove(&mut self, mut poly: VirtualPolynomial<F>) -> ZerocheckProof<F, E, T, PCS> {
        debug_assert!(self.lagrange_coefficients.len() == poly.degree() + 2);

        // We first need to collect how many multilinears we have so that we only make evaluations
        // of these and not include the eq polynomial.
        let n_polys = poly.constituents.len();

        // Perform the IOP.
        let (zerocheck_proof, eval_point) =
            zerocheck::prove(&mut poly, &mut self.transcript, &self.lagrange_coefficients);

        // With a generated evaluation point in hand, we can create the full polynomial
        // evaluations, batch commitment and proof.
        let evaluations = poly
            .constituents
            .par_iter()
            .take(n_polys)
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
            .commit(&poly.constituents[..n_polys], &mut self.transcript);

        let proof = self.pcs.prove(
            &commitment,
            &poly.constituents[..n_polys],
            &eval_point,
            &mut self.transcript,
        );

        ZerocheckProof {
            zerocheck_proof,
            evaluations,
            commitment,
            proof,
            products: poly.products.clone(),
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
