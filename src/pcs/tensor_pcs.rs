//! Tensor PCS from [DP23] https://eprint.iacr.org/2023/1784.pdf, section 3.3 construction 3.7.

use crate::{
    field::{ChallengeField, Field},
    linear_code::LinearCode,
    merkle_tree::MerkleTree,
    pcs::PolynomialCommitmentScheme,
    polynomial::MultilinearExtension,
    transcript::{IntoObservable, Transcript},
};
use blake2::{Blake2s256, Digest};
use core::marker::PhantomData;
use rayon::prelude::*;

/// An implementation of Diamond and Posen's tensor polynomial commitment scheme, which is linear
/// in prover complexity and has easily configurable soundness.
pub struct TensorPCS<F: Field, E: ChallengeField<F>, T: Transcript<F>, LC: LinearCode<F>> {
    n_test_queries: usize,
    code: LC,
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
    _t_marker: PhantomData<T>,
}

/// A commitment made with the tensor polynomial commitment scheme.
pub struct TensorCommitment<F: Field> {
    // Merkle tree where each leaf `n` is the hash of the column at `n` in each polynomial matrix.
    tree: MerkleTree,
    // Inputted polynomials reshaped as matrices and encoded row-wise with a linear code.
    encoded_matrices: Vec<Vec<F>>,
}

impl<F: Field> IntoObservable for TensorCommitment<F> {
    fn into_observable(&self) -> Vec<&[u8; 32]> {
        self.tree.to_hashes()
    }
}

/// An evaluation proof made with the tensor polynomial commitment scheme.
pub struct TensorProof<F: Field, E: ChallengeField<F>> {
    // The mixed product of all polynomials committed to.
    t_prime: Vec<E>,
    // Authentication paths for the query indices.
    authentication_paths: Vec<Vec<[[u8; 32]; 2]>>,
    // The columns of the encoded matrices at the query indices.
    columns: Vec<Vec<Vec<F>>>,
}

impl<F: Field, E: ChallengeField<F>, T: Transcript<F>, LC: LinearCode<F>> TensorPCS<F, E, T, LC> {
    pub fn new(n_test_queries: usize, code: LC) -> Self {
        Self {
            n_test_queries,
            code,
            _f_marker: PhantomData::<F>,
            _e_marker: PhantomData::<E>,
            _t_marker: PhantomData::<T>,
        }
    }
}

impl<F: Field, E: ChallengeField<F>, T: Transcript<F>, LC: LinearCode<F>> TensorPCS<F, E, T, LC> {
    fn encode_ext(&self, coeffs: &[E]) -> Vec<E> {
        // We unpack the coeffs into E::DEGREE different vectors and compute low-degree extensions
        // of each, after which we zip them back up into extension field elements.
        let unpacked = coeffs
            .iter()
            .map(|e| Into::<Vec<F>>::into(*e))
            .collect::<Vec<Vec<F>>>();
        let unpacked = (0..E::DEGREE)
            .map(|i| unpacked.iter().map(|el| el[i]).collect::<Vec<F>>())
            .collect::<Vec<Vec<F>>>();
        let mut encoded_columns = vec![];
        for column in unpacked {
            encoded_columns.push(self.code.encode(&column));
        }

        (0..encoded_columns[0].len())
            .map(|i| {
                let mut recoded_el = Vec::with_capacity(E::DEGREE);
                (0..E::DEGREE).for_each(|j| {
                    recoded_el.push(encoded_columns[j][i]);
                });
                E::new(recoded_el)
            })
            .collect::<Vec<E>>()
    }
}

impl<F: Field, E: ChallengeField<F>, T: Transcript<F>, LC: LinearCode<F>>
    PolynomialCommitmentScheme<F, E, T> for TensorPCS<F, E, T, LC>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    type Commitment = TensorCommitment<F>;
    type Proof = TensorProof<F, E>;

    fn commit(&self, polys: &[MultilinearExtension<F>], transcript: &mut T) -> Self::Commitment {
        debug_assert!(polys.iter().all(|p| p.len() == polys[0].len()));

        // Turn the polys into matrices, then encode row-wise.
        // We take the next-power-of-two as the number of rows.
        let log_size = polys[0].len().isqrt().next_power_of_two();
        let encoded_matrices = polys
            .par_iter()
            .map(|poly| {
                poly.evals
                    .chunks(log_size)
                    .flat_map(|chunk| self.code.encode(chunk))
                    .collect::<Vec<F>>()
            })
            .collect::<Vec<Vec<F>>>();

        // Create the column hashes for each matrix.
        let row_size = log_size << LC::BLOWUP_BITS;
        let depth = polys[0].len() / row_size;
        let leaves = encoded_matrices
            .par_iter()
            .map(|matrix| {
                (0..row_size)
                    .into_iter()
                    .map(|i| {
                        let mut hasher = Blake2s256::new();
                        (0..depth).for_each(|j| {
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
        TensorCommitment {
            tree: MerkleTree::new(leaves),
            encoded_matrices,
        }
    }

    // We assume that the merkle tree has been observed into the transcript by now.
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
        let log_size = polys[0].len().isqrt().next_power_of_two();
        let outer_expansion = tensor_product_expansion(&eval[(log_size.ilog2() as usize)..]);

        // Calculate the t' for each t (poly) with the vector-matrix product.
        let mut t_primes = vec![vec![E::ZERO; log_size]; polys.len()];
        t_primes
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, t_prime)| {
                outer_expansion.iter().enumerate().for_each(|(j, tensor)| {
                    t_prime.iter_mut().enumerate().for_each(|(k, entry)| {
                        let mut tensor = tensor.clone();
                        tensor.mul_base(&polys[i].evals[j * log_size + k]);
                        entry.add_assign(&tensor);
                    });
                });
            });

        // Mix all t' into one using challenges.
        let n_challenges = (comm.encoded_matrices.len() as f32).log2().ceil() as usize;
        let t_prime: Vec<E> = if n_challenges > 0 {
            let mut challenges = Vec::with_capacity(n_challenges);
            for _ in 0..n_challenges {
                challenges.push(transcript.draw_challenge_ext());
            }

            let expansion = tensor_product_expansion(&challenges)[..t_primes.len()].to_vec();
            expansion.iter().enumerate().for_each(|(i, coeff)| {
                t_primes[i].iter_mut().for_each(|e| e.mul_assign(coeff));
            });

            t_primes[0]
                .clone()
                .into_par_iter()
                .enumerate()
                .map(|(i, e)| {
                    (1..t_primes.len()).fold(e, |mut acc, j| {
                        acc.add_assign(&t_primes[j][i]);
                        acc
                    })
                })
                .collect::<Vec<E>>()
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
        let row_size = (log_size << LC::BLOWUP_BITS) as usize;
        let row_size_bits = row_size.ilog2() as usize;
        let depth = comm.encoded_matrices[0].len() / row_size;
        let (authentication_paths, columns): (Vec<Vec<[[u8; 32]; 2]>>, Vec<Vec<Vec<F>>>) = (0
            ..self.n_test_queries)
            .into_iter()
            .map(|_| {
                // Sample a random index for a given column.
                let index = transcript.draw_bits(row_size_bits);
                let path = comm.tree.get_proof(index);
                // Extract a column per polynomial that we use in the proof.
                let cols = comm
                    .encoded_matrices
                    .par_iter()
                    .map(|matrix| {
                        (0..depth)
                            .map(|i| matrix[i * row_size + index])
                            .collect::<Vec<F>>()
                    })
                    .collect::<Vec<Vec<F>>>();
                (path, cols)
            })
            .collect::<Vec<_>>()
            .into_iter()
            .unzip();

        TensorProof {
            t_prime,
            authentication_paths,
            columns,
        }
    }

    // We assume similarly that the verifier has observed the commitment in question into the
    // transcript.
    fn verify(
        &self,
        comm: &Self::Commitment,
        eval: &[E],
        results: &[E],
        proof: &Self::Proof,
        transcript: &mut T,
    ) -> bool {
        let log_size = proof.t_prime.len();
        let inner_expansion = tensor_product_expansion(&eval[..(log_size.ilog2() as usize)]);
        let outer_expansion = tensor_product_expansion(&eval[(log_size.ilog2() as usize)..]);

        // Create batch value.
        let n_challenges = (comm.encoded_matrices.len() as f32).log2().ceil() as usize;
        let mut expansion = vec![];
        let result = if n_challenges > 0 {
            let mut challenges: Vec<E> = Vec::with_capacity(n_challenges);
            for _ in 0..n_challenges {
                challenges.push(transcript.draw_challenge_ext());
            }

            expansion = tensor_product_expansion(&challenges)[..results.len()].to_vec();

            // Return the inner product of the expanded challenges and the results.
            expansion
                .clone()
                .into_iter()
                .zip(results.iter())
                .map(|(a, b)| {
                    let mut a = a.clone();
                    a.mul_assign(b);
                    a
                })
                .fold(E::ZERO, |mut acc, x| {
                    acc.add_assign(&x);
                    acc
                })
        } else {
            results[0].clone()
        };

        // Enter t_prime into the transcript.
        transcript.observe_witnesses(
            &proof
                .t_prime
                .iter()
                .flat_map(|e| Into::<Vec<F>>::into(*e))
                .collect::<Vec<F>>(),
        );

        // Encode t_prime.
        let enc_t_prime = self.encode_ext(&proof.t_prime);

        // Ensure that merkle paths and column evaluations are correct.
        for i in 0..self.n_test_queries {
            let index = transcript.draw_bits(enc_t_prime.len().ilog2() as usize);
            if !verify_val(
                &proof.columns[i],
                &outer_expansion,
                enc_t_prime[index],
                &expansion,
            ) || !comm.tree.verify_path(&proof.authentication_paths[i])
            {
                return false;
            }
        }

        // Ensure that t_prime times inner tensor expansion equals result.
        inner_expansion
            .iter()
            .zip(proof.t_prime.iter())
            .fold(E::ZERO, |mut acc, (t, e)| {
                let mut t = t.clone();
                t.mul_assign(e);
                acc.add_assign(&t);
                acc
            })
            == result
    }
}

fn verify_val<F: Field, E: ChallengeField<F>>(
    columns: &[Vec<F>],
    tensor: &[E],
    eval: E,
    mixing_coeffs: &[E],
) -> bool {
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
        debug_assert!(evals.len() == mixing_coeffs.len());
        mixing_coeffs
            .into_iter()
            .zip(evals.iter())
            .map(|(a, b)| {
                let mut a = a.clone();
                a.mul_assign(b);
                a
            })
            .fold(E::ZERO, |mut acc, x| {
                acc.add_assign(&x);
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
        fft::CircleFFT,
        field::{
            m31::{quartic::M31_4, M31},
            Field,
        },
        linear_code::reed_solomon::ReedSolomonCode,
        test_utils::rand_poly,
        transcript::Blake2sTranscript,
    };
    use rand::Rng;

    const POLY_SIZE_BITS: u32 = 20;
    const ROOTS_OF_UNITY_BITS: usize = (POLY_SIZE_BITS / 2 + 1) as usize;
    const N_QUERIES: usize = 143; // roughly 96 bits of security with an RS code that has a
                                  // blowup of 4; see binius

    #[test]
    fn commit_prove_verify_single_poly() {
        let poly = rand_poly(POLY_SIZE_BITS);

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let commitment = pcs.commit(&[poly.clone().into()], &mut prover_transcript);
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let proof = pcs.prove(
            &commitment,
            &[poly.clone().into()],
            &eval,
            &mut prover_transcript,
        );

        let mut res = poly.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res.fix_variable(*e);
        });

        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(pcs.verify(
            &commitment,
            &eval,
            &[res.evals[0]],
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn commit_prove_verify_single_poly_wrong_eval() {
        let poly = rand_poly(POLY_SIZE_BITS);

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let commitment = pcs.commit(&[poly.clone().into()], &mut prover_transcript);
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let proof = pcs.prove(
            &commitment,
            &[poly.clone().into()],
            &eval,
            &mut prover_transcript,
        );

        let mut res = poly.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res.fix_variable(*e);
        });
        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(!pcs.verify(
            &commitment,
            &eval,
            &[res.evals[0]],
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn commit_prove_verify_single_poly_wrong_commitment() {
        let poly = rand_poly(POLY_SIZE_BITS);

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let commitment = pcs.commit(&[poly.clone().into()], &mut Blake2sTranscript::default());
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let fake_poly = rand_poly(POLY_SIZE_BITS);
        let fake_commitment = pcs.commit(&[fake_poly.clone().into()], &mut prover_transcript);
        let proof = pcs.prove(
            &fake_commitment,
            &[poly.clone().into()],
            &eval,
            &mut prover_transcript,
        );

        let mut res = poly.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res.fix_variable(*e);
        });
        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(!pcs.verify(
            &commitment,
            &eval,
            &[res.evals[0]],
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn commit_prove_verify_single_poly_wrong_proof() {
        let poly = rand_poly(POLY_SIZE_BITS);

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let commitment = pcs.commit(&[poly.clone().into()], &mut prover_transcript);
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let fake_poly = rand_poly(POLY_SIZE_BITS);
        let proof = pcs.prove(
            &commitment,
            &[fake_poly.clone().into()],
            &eval,
            &mut prover_transcript,
        );

        let mut res = poly.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res.fix_variable(*e);
        });
        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(!pcs.verify(
            &commitment,
            &eval,
            &[res.evals[0]],
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn commit_prove_verify_2_poly() {
        let poly = rand_poly(POLY_SIZE_BITS);
        let poly_2 = rand_poly(POLY_SIZE_BITS);

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let commitment = pcs.commit(&[poly.clone(), poly_2.clone()], &mut prover_transcript);
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let proof = pcs.prove(
            &commitment,
            &[poly.clone(), poly_2.clone()],
            &eval,
            &mut prover_transcript,
        );

        let mut res = poly.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res.fix_variable(*e);
        });
        let mut res_2 = poly_2.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res_2.fix_variable(*e);
        });
        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(pcs.verify(
            &commitment,
            &eval,
            &[res.evals[0], res_2.evals[0]],
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn commit_prove_verify_3_poly() {
        let poly = rand_poly(POLY_SIZE_BITS);
        let poly_2 = rand_poly(POLY_SIZE_BITS);
        let poly_3 = rand_poly(POLY_SIZE_BITS);

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let commitment = pcs.commit(
            &[poly.clone(), poly_2.clone(), poly_3.clone()],
            &mut prover_transcript,
        );
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); poly.num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let proof = pcs.prove(
            &commitment,
            &[poly.clone(), poly_2.clone(), poly_3.clone()],
            &eval,
            &mut prover_transcript,
        );

        let mut res = poly.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res.fix_variable(*e);
        });
        let mut res_2 = poly_2.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res_2.fix_variable(*e);
        });
        let mut res_3 = poly_3.fix_variable_ext(eval[0]);
        eval.iter().skip(1).for_each(|e| {
            res_3.fix_variable(*e);
        });
        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(pcs.verify(
            &commitment,
            &eval,
            &[res.evals[0], res_2.evals[0], res_3.evals[0]],
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn commit_prove_verify_32_poly() {
        let polys = (0..32)
            .map(|_| rand_poly(POLY_SIZE_BITS))
            .collect::<Vec<MultilinearExtension<M31>>>();

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let now = std::time::Instant::now();
        let commitment = pcs.commit(&polys, &mut prover_transcript);
        let elapsed = std::time::Instant::now();
        println!("commitment takes {:?}", elapsed - now);
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); polys[0].num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let now = std::time::Instant::now();
        let proof = pcs.prove(&commitment, &polys, &eval, &mut prover_transcript);
        let elapsed = std::time::Instant::now();
        println!("proof takes {:?}", elapsed - now);

        let results = (0..32)
            .map(|i| {
                let mut res = polys[i].fix_variable_ext(eval[0]);
                eval.iter().skip(1).for_each(|e| {
                    res.fix_variable(*e);
                });

                res.evals[0]
            })
            .collect::<Vec<M31_4>>();
        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(pcs.verify(
            &commitment,
            &eval,
            &results,
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn commit_prove_verify_256_poly() {
        let polys = (0..256)
            .map(|_| rand_poly(POLY_SIZE_BITS))
            .collect::<Vec<MultilinearExtension<M31>>>();

        let mut prover_transcript = Blake2sTranscript::default();
        let pcs = TensorPCS::<M31, M31_4, Blake2sTranscript<M31>, ReedSolomonCode<M31, CircleFFT>> {
            n_test_queries: N_QUERIES,
            code: ReedSolomonCode::new(ROOTS_OF_UNITY_BITS),
            _f_marker: PhantomData::<M31>,
            _e_marker: PhantomData::<M31_4>,
            _t_marker: PhantomData::<Blake2sTranscript<M31>>,
        };
        let now = std::time::Instant::now();
        let commitment = pcs.commit(&polys, &mut prover_transcript);
        let elapsed = std::time::Instant::now();
        println!("commitment takes {:?}", elapsed - now);
        prover_transcript.observe_hashes(&commitment.into_observable());

        let mut eval = vec![M31_4::default(); polys[0].num_vars()];
        eval.iter_mut().for_each(|e| {
            *e = M31_4::from_usize(rand::thread_rng().gen_range(0..M31::ORDER) as usize);
        });
        let now = std::time::Instant::now();
        let proof = pcs.prove(&commitment, &polys, &eval, &mut prover_transcript);
        let elapsed = std::time::Instant::now();
        println!("proof takes {:?}", elapsed - now);

        let results = (0..256)
            .map(|i| {
                let mut res = polys[i].fix_variable_ext(eval[0]);
                eval.iter().skip(1).for_each(|e| {
                    res.fix_variable(*e);
                });

                res.evals[0]
            })
            .collect::<Vec<M31_4>>();
        let mut verifier_transcript = Blake2sTranscript::default();
        verifier_transcript.observe_hashes(&commitment.into_observable());
        assert!(pcs.verify(
            &commitment,
            &eval,
            &results,
            &proof,
            &mut verifier_transcript
        ));
    }
}
