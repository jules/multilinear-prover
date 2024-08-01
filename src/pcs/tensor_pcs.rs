use crate::{
    field::{ChallengeField, Field},
    linear_code::LinearCode,
    merkle_tree::MerkleTree,
    mle::MultilinearExtension,
    pcs::PolynomialCommitmentScheme,
    transcript::Transcript,
};
use blake2::{Blake2s256, Digest};
use core::{
    hash::{Hash, Hasher},
    marker::PhantomData,
};

/// [DP23]: https://eprint.iacr.org/2023/630.pdf
pub struct TensorPCS<F: Field, LC: LinearCode<F>> {
    _f_marker: PhantomData<F>,
    _lc_marker: PhantomData<LC>,
}

impl<F: Field, LC: LinearCode<F>, T: Transcript<F>, E: ChallengeField<F>>
    PolynomialCommitmentScheme<F, T, E> for TensorPCS<F, LC>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    type Commitment = ([u8; 32], Vec<Vec<F>>);
    type Proof = Vec<[u8; 32]>;

    fn commit(polys: &[MultilinearExtension<F>]) -> Self::Commitment {
        // Turn the polys into m x m matrices, then encode row-wise.
        // XXX: ensure same length
        let log_size = polys[0].evals.len() >> 1;
        let matrices = polys
            .iter()
            .map(|poly| {
                poly.evals
                    .chunks(log_size as usize)
                    .flat_map(|chunk| LC::encode(chunk))
                    .collect::<Vec<F>>()
            })
            .collect::<Vec<Vec<F>>>();

        // Build a merkle tree out of our encoded matrices. We populate the base layer with a hash
        // of each row and then work up.
        let row_size = log_size * LC::BLOWUP;

        // The merkle tree impl takes care of the layer-on-layer hashing.
        let tree = MerkleTree::new(matrices.clone(), row_size as usize);

        (tree.root(), matrices)
    }

    // for now, we assume commitment was observed by transcript
    fn prove(
        comm: &Self::Commitment,
        polys: &[MultilinearExtension<F>],
        eval: Vec<E>,
        result: E,
        transcript: &mut T,
    ) -> Self::Proof {
        debug_assert!(eval.len() == polys[0].num_vars());

        let mut challenges = Vec::with_capacity(comm.1.len());
        for _ in 0..comm.1.len() {
            challenges.push(transcript.draw_challenge_ext());
        }

        // We want to compute the tensor product expansion of the top half of the evaluation
        // variables. We then compute the vector-matrix product of this tensor product expansion
        // paired with each polynomial.
        // XXX HIGHLY UNOPTIMIZED
        let expansion = {
            let mut buf = vec![E::ONE];
            for e in (eval[eval.len() / 2..]) {
                let buf_1 = buf
                    .iter()
                    .map(|v| {
                        let mut prod = v;
                        prod.mul_assign(&e);
                        let mut s = v;
                        s.sub_assign(&prod);
                        s
                    })
                    .collect::<Vec<E>>();
                let buf_2 = buf
                    .iter()
                    .map(|v| {
                        let mut prod = v;
                        prod.mul_assign(&e);
                        prod
                    })
                    .collect::<Vec<E>>();
                buf = buf_1.iter().chain(buf_2.iter()).collect();
            }

            buf
        };

        let log_size = polys[0].evals.len() >> 1;

        // Calculate the t' for each t (poly) with the vector-matrix product.
        let t_primes = polys.iter().zip(challenges.iter()).map(|(poly, c)| {
            poly.evals.chunks(log_size as usize).flat_map(|chunk|

        // Compute merkle proofs for every row in the committed matrix. Extract columns to include
        // in the proof.
        
    }

    fn verify(
        comm: &Self::Commitment,
        eval: Vec<E>,
        result: E,
        proof: Self::Proof,
        transcript: &mut T,
    ) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        field::m31::M31, linear_code::reed_solomon::ReedSolomonCode, mle::MultilinearExtension,
    };
    use rand::Rng;

    #[test]
    fn commit_single_poly() {
        let mut evals = vec![M31::default(); 2u32.pow(20) as usize];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = MultilinearExtension::new(evals);

        let commitment = TensorPCS::<M31, ReedSolomonCode<M31>>::commit(&[poly]);
    }
}
