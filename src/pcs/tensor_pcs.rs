//! Tensor PCS from [DP23] https://eprint.iacr.org/2023/630.pdf.

use crate::{
    field::{ChallengeField, Field},
    linear_code::LinearCode,
    merkle_tree::{verify_path, MerkleTree},
    pcs::PolynomialCommitmentScheme,
    polynomial::mle::MultilinearExtension,
    transcript::Transcript,
};
use blake2::{Blake2s256, Digest};
use core::marker::PhantomData;
use rayon::prelude::*;

pub struct TensorPCS<F: Field, T: Transcript<F>, E: ChallengeField<F>, LC: LinearCode<F, E>> {
    n_test_queries: usize,
    _f_marker: PhantomData<F>,
    _t_marker: PhantomData<T>,
    _e_marker: PhantomData<E>,
    _lc_marker: PhantomData<LC>,
}

impl<F: Field, T: Transcript<F>, E: ChallengeField<F>, LC: LinearCode<F, E>>
    TensorPCS<F, T, E, LC>
{
    pub fn new(n_test_queries: usize) -> Self {
        Self {
            n_test_queries,
            _f_marker: PhantomData::<F>,
            _t_marker: PhantomData::<T>,
            _e_marker: PhantomData::<E>,
            _lc_marker: PhantomData::<LC>,
        }
    }
}

impl<F: Field, LC: LinearCode<F, E>, T: Transcript<F>, E: ChallengeField<F>>
    PolynomialCommitmentScheme<F, T, E> for TensorPCS<F, T, E, LC>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    type Commitment = (MerkleTree<F>, Vec<Vec<F>>);
    type Proof = (Vec<E>, Vec<(Vec<Vec<[u8; 32]>>, Vec<Vec<F>>)>);

    fn commit(&self, polys: &[MultilinearExtension<F>], _transcript: &mut T) -> Self::Commitment {
        debug_assert!(polys.iter().all(|p| p.evals.len() == polys[0].evals.len()));

        // Turn the polys into m x m matrices row-wise, then encode row-wise.
        let log_size = polys[0].evals.len().isqrt();
        let matrices = polys
            .iter()
            .map(|poly| {
                poly.evals
                    .par_chunks(log_size)
                    .flat_map(|chunk| LC::encode(chunk))
                    .collect::<Vec<F>>()
            })
            .collect::<Vec<Vec<F>>>();

        // Create the column hashes for each matrix.
        let row_size = log_size * LC::BLOWUP;
        let leaves = matrices
            .iter()
            .map(|matrix| {
                (0..row_size)
                    .map(|i| {
                        let mut hasher = Blake2s256::new();
                        (0..log_size).for_each(|j| {
                            hasher.update(matrix[i + j * row_size].to_le_bytes().to_vec());
                        });
                        <[u8; 32]>::from(hasher.finalize())
                    })
                    .collect::<Vec<[u8; 32]>>()
            })
            .collect::<Vec<Vec<[u8; 32]>>>();

        // Build a merkle tree out of our encoded matrices. We populate the base layer with a hash
        // of each column and then work up.
        // The merkle tree impl takes care of the layer-on-layer hashing.
        (
            MerkleTree::new(leaves, row_size as usize, log_size),
            matrices,
        )
    }

    // for now, we assume commitment was observed by transcript
    fn prove(
        &self,
        comm: &Self::Commitment,
        polys: &[MultilinearExtension<F>],
        eval: &[E],
        transcript: &mut T,
    ) -> Self::Proof {
        debug_assert!(eval.len() == polys[0].num_vars());

        // We want to compute the tensor product expansion of the top half of the evaluation
        // variables. We then compute the vector-matrix product of this tensor product expansion
        // paired with each polynomial.
        let outer_expansion = tensor_product_expansion(&eval[eval.len() / 2..]);

        // Calculate the t' for each t (poly) with the vector-matrix product.
        let log_size = polys[0].evals.len().isqrt();
        let mut t_primes = vec![vec![E::ZERO; log_size]; polys.len()];
        t_primes.iter_mut().enumerate().for_each(|(i, t_prime)| {
            outer_expansion.iter().enumerate().for_each(|(j, tensor)| {
                t_prime.iter_mut().enumerate().for_each(|(k, entry)| {
                    let mut tensor = tensor.clone();
                    tensor.mul_base(&polys[i].evals[j * log_size + k]);
                    entry.add_assign(&tensor);
                });
            });
        });

        // Mix all t' into one using challenges.
        let n_challenges = comm.1.len().ilog2() as usize;
        let t_prime: Vec<E> = if n_challenges > 0 {
            //let mut challenges = Vec::with_capacity(n_challenges);
            //for _ in 0..n_challenges {
            //    challenges.push(transcript.draw_challenge_ext());
            //}

            // MIX
            t_primes.into_iter().flatten().collect() // TODO
        } else {
            t_primes.into_iter().flatten().collect()
        };

        // Enter t_prime into the transcript.
        transcript.observe_witnesses(
            &t_prime
                .iter()
                .flat_map(|e| Into::<Vec<F>>::into(*e))
                .collect::<Vec<F>>(),
        );

        // Compute merkle proofs for n rows in the committed matrix. Extract columns to include
        // in the proof.
        // XXX figure out optimal n
        let row_size = (log_size * LC::BLOWUP) as usize;
        (
            t_prime,
            (0..self.n_test_queries)
                .map(|_| {
                    // Sample a random index for a given row.
                    let index = transcript.draw_bits(row_size);
                    let path = comm.0.get_proof(index);
                    // Extract a column per polynomial that we use in the proof.
                    let cols = comm
                        .1
                        .iter()
                        .map(|matrix| {
                            (0..log_size)
                                .map(|i| matrix[i * row_size])
                                .collect::<Vec<F>>()
                        })
                        .collect::<Vec<Vec<F>>>();
                    (path, cols)
                })
                .collect::<Vec<_>>(),
        )
    }

    fn verify(
        &self,
        comm: &Self::Commitment,
        eval: &[E],
        result: E,
        proof: &Self::Proof,
        transcript: &mut T,
    ) -> bool {
        let inner_expansion = tensor_product_expansion(&eval[..eval.len() / 2]);
        let outer_expansion = tensor_product_expansion(&eval[eval.len() / 2..]);

        // Enter t_prime into the transcript.
        transcript.observe_witnesses(
            &proof
                .0
                .iter()
                .flat_map(|e| Into::<Vec<F>>::into(*e))
                .collect::<Vec<F>>(),
        );

        // Encode t_prime.
        let enc_t_prime = LC::encode_ext(&proof.0);

        // Ensure that merkle paths and column evaluations are correct.
        let log_size = proof.0.len();
        let row_size = (log_size * LC::BLOWUP) as usize;
        for i in 0..self.n_test_queries {
            let index = transcript.draw_bits(row_size);
            if !verify_val::<F, E>(&proof.1[i].1, &outer_expansion, enc_t_prime[index])
                || !verify_path(&proof.1[i].0)
            {
                return false;
            }
        }

        // Ensure that t_prime times inner tensor expansion equals result.
        inner_expansion
            .iter()
            .zip(proof.0.iter())
            .fold(E::ZERO, |mut acc, (t, e)| {
                let mut t = t.clone();
                t.mul_assign(e);
                acc.add_assign(&t);
                acc
            })
            == result
    }
}

fn verify_val<F: Field, E: ChallengeField<F>>(columns: &[Vec<F>], tensor: &[E], eval: E) -> bool {
    let evals = columns
        .iter()
        .map(|column| {
            tensor
                .iter()
                .zip(column.iter())
                .fold(E::ZERO, |mut acc, (t, e)| {
                    let mut t = t.clone();
                    t.mul_base(e);
                    acc.add_assign(&t);
                    acc
                })
        })
        .collect::<Vec<E>>();

    // Mix if we have multiple polys.
    let final_eval = if evals.len() > 1 {
        evals.iter().fold(E::ZERO, |mut acc, e| {
            // MIX TODO
            acc.add_assign(e);
            acc
        })
    } else {
        evals[0]
    };

    final_eval == eval
}

fn tensor_product_expansion<F: Field, E: ChallengeField<F>>(eval: &[E]) -> Vec<E> {
    let mut buf = vec![E::ONE];
    for e in eval {
        let middle = buf.len();
        buf.resize(middle << 1, E::ZERO);
        let (left, right) = buf.split_at_mut(middle);
        left.iter_mut().zip(right.iter_mut()).for_each(|(l, r)| {
            let mut prod = l.clone();
            prod.mul_assign(&e);
            l.sub_assign(&prod);
            *r = prod;
        });
    }

    buf
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        field::{
            m31::{quartic::M31_4, M31},
            Field,
        },
        linear_code::reed_solomon::ReedSolomonCode,
        polynomial::mle::MultilinearExtension,
        test_utils::{rand_poly, MockTranscript},
        transcript::Transcript,
    };
    use rand::Rng;

    #[test]
    fn commit_prove_verify_single_poly() {
        let poly = rand_poly(2u32.pow(20) as usize);

        let pcs = TensorPCS::<M31, MockTranscript<M31>, M31_4, ReedSolomonCode<M31, M31_4>> {
            n_test_queries: 4,
            _f_marker: PhantomData::<M31>,
            _t_marker: PhantomData::<MockTranscript<M31>>,
            _e_marker: PhantomData::<M31_4>,
            _lc_marker: PhantomData::<ReedSolomonCode<M31, M31_4>>,
        };
        let commitment = pcs.commit(&[poly.clone()], &mut MockTranscript::default());

        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let proof = pcs.prove(
            &commitment,
            &[poly.clone()],
            &eval,
            &mut MockTranscript::default(),
        );

        let mut res = poly.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res.fix_variable(*e);
        });
        assert!(pcs.verify(
            &commitment,
            &eval,
            res.evals[0],
            &proof,
            &mut MockTranscript::default()
        ));
    }
}
