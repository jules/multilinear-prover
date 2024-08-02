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
    pub fn new(elements: Vec<Vec<F>>, row_size: usize, col_size: usize) -> Self {
        debug_assert!(elements[0].len().is_power_of_two());
        // Create initial leaves.
        let leaves = elements
            .iter()
            .flat_map(|matrix| {
                (0..row_size)
                    .map(|i| {
                        let mut hasher = Blake2s256::new();
                        (0..col_size).for_each(|j| {
                            hasher.update(matrix[i + j * row_size].to_le_bytes().to_vec());
                        });
                        <[u8; 32]>::from(hasher.finalize())
                    })
                    .collect::<Vec<[u8; 32]>>()
            })
            .collect::<Vec<[u8; 32]>>();

        // Perform hash-chaining up to root.
        let mut size = leaves.len();
        let mut layers = vec![leaves];
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

    pub fn get_proof(&self, mut index: usize) -> Vec<[u8; 32]> {
        let mut path = Vec::with_capacity(self.elements.len());
        path.push(self.elements[0][index]); // leaf
        for i in 1..self.elements.len() {
            index >>= 1;
            path.push(self.elements[i][index]);
        }

        path
    }

    pub fn root(&self) -> [u8; 32] {
        self.root
    }
}
