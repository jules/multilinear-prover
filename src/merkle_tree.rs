use crate::field::Field;
use blake2::{Blake2s256, Digest};
use core::marker::PhantomData;

pub struct MerkleTree<F: Field> {
    root: [u8; 32],
    elements: Vec<Vec<[u8; 32]>>, // bottom layer first
    _marker: PhantomData<F>,
}

impl<F: Field> MerkleTree<F>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    /// Constructs a new MerkleTree from matrices of elements. Each row is hashed and then 2-arity
    /// hashing up to the root is performed.
    pub fn new(mut leaves: Vec<Vec<[u8; 32]>>, row_size: usize, col_size: usize) -> Self {
        debug_assert!(leaves[0].len().is_power_of_two());
        // Create initial leaves. If we only have one layer, do nothing. Otherwise, we hash the
        // i'th entry of each layer together to make a single leaf.
        let compressed_leaves = if leaves.len() > 1 {
            (0..col_size)
                .map(|i| {
                    let mut hasher = Blake2s256::new();
                    leaves.iter().for_each(|layer| {
                        hasher.update(layer[i]);
                    });
                    <[u8; 32]>::from(hasher.finalize())
                })
                .collect::<Vec<[u8; 32]>>()
        } else {
            leaves.pop().unwrap()
        };

        // Perform hash-chaining up to root.
        let mut size = compressed_leaves.len();
        let mut layers = vec![compressed_leaves];
        loop {
            if size == 1 {
                break;
            }

            let last = layers.last().unwrap();
            let new_layer = last
                .chunks(2)
                .map(|c| {
                    let mut hasher = Blake2s256::new();
                    hasher.update(c[0]);
                    hasher.update(c[1]);
                    <[u8; 32]>::from(hasher.finalize())
                })
                .collect::<Vec<[u8; 32]>>();
            layers.push(new_layer);
            size >>= 1;
        }

        Self {
            root: layers.last().unwrap()[0],
            elements: layers,
            _marker: PhantomData::<F>,
        }
    }

    pub fn get_proof(&self, mut index: usize) -> Vec<Vec<[u8; 32]>> {
        let mut path = Vec::with_capacity(self.elements.len());
        let mut base_layer = vec![];
        base_layer.push(self.elements[0][index]); // leaf
        base_layer.push(self.elements[0][index ^ 1]); // neighboring element
        path.push(base_layer);
        for i in 1..(self.elements.len() - 1) {
            let mut layer = vec![];
            index >>= 1;
            layer.push(self.elements[i][index]);
            layer.push(self.elements[i][index ^ 1]);
            path.push(layer);
        }
        path.push(vec![self.elements[self.elements.len() - 1][0]]);

        path
    }

    pub fn root(&self) -> [u8; 32] {
        self.root
    }
}

pub fn verify_path(path: Vec<Vec<[u8; 32]>>) -> bool {
    for i in 0..(path.len() - 1) {
        let mut hasher = Blake2s256::new();
        path[i].iter().for_each(|e| hasher.update(e));
        let hash = <[u8; 32]>::from(hasher.finalize());
        if path[i + 1].iter().all(|e| *e != hash) {
            return false;
        }
    }

    true
}
