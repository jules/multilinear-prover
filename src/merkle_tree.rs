use crate::field::Field;
use blake2::{Blake2s256, Digest};
use core::marker::PhantomData;

pub struct MerkleTree<F: Field> {
    root: [u8; 32],
    elements: Vec<Vec<[u8; 32]>>,
    _marker: PhantomData<F>,
}

impl<F: Field> MerkleTree<F>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    pub fn new(elements: Vec<Vec<F>>, row_size: usize) -> Self {
        debug_assert!(elements[0].len().is_power_of_two());
        // Create initial leaves.
        let leaves = elements
            .iter()
            .flat_map(|matrix| {
                matrix
                    .chunks(row_size)
                    .map(|row| {
                        let mut hasher = Blake2s256::new();
                        row.iter()
                            .for_each(|el| hasher.update(el.to_le_bytes().to_vec()));
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

    pub fn root(&self) -> [u8; 32] {
        self.root
    }
}
