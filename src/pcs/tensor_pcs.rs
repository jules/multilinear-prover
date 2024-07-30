use crate::{
    field::Field, linear_code::LinearCode, merkle_tree::MerkleTree, mle::MultilinearExtension,
    pcs::PolynomialCommitmentScheme,
};
use blake2::{Blake2s256, Digest};
use core::{
    hash::{Hash, Hasher},
    marker::PhantomData,
};

pub struct TensorPCS<F: Field, LC: LinearCode<F>> {
    _f_marker: PhantomData<F>,
    _lc_marker: PhantomData<LC>,
}

impl<F: Field, LC: LinearCode<F>> PolynomialCommitmentScheme<F> for TensorPCS<F, LC>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    type Commitment = ([u8; 32], Vec<Vec<F>>);
    type Proof = Vec<[u8; 32]>;

    fn commit(polys: &[MultilinearExtension<F>]) -> Self::Commitment {
        // Turn the polys into m x m matrices, then encode row-wise.
        // XXX: ensure same length
        let log_size = polys.len().ilog2();
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
        let leaves = matrices
            .iter()
            .flat_map(|matrix| {
                matrix
                    .chunks(row_size as usize)
                    .map(|row| {
                        let mut hasher = Blake2s256::new();
                        row.iter()
                            .for_each(|el| hasher.update(el.to_le_bytes().to_vec()));
                        <[u8; 32]>::from(hasher.finalize())
                    })
                    .collect::<Vec<[u8; 32]>>()
            })
            .collect::<Vec<[u8; 32]>>();

        // The merkle tree impl takes care of the layer-on-layer hashing.
        let tree = MerkleTree::new(leaves);

        (tree.root(), matrices)
    }

    fn open(comm: &Self::Commitment, eval: Vec<F>, result: F) -> Self::Proof {
        panic!()
    }

    fn verify(comm: &Self::Commitment, eval: Vec<F>, result: F, proof: Self::Proof) -> bool {
        false
    }
}
