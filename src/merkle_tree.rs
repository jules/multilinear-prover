use crate::field::m31::M31;
use blake2::{Blake2s256, Digest};

pub struct MerkleTree {
    root: [u8; 32],
    elements: Vec<Vec<[u8; 32]>>,
}

impl MerkleTree {
    pub fn new(leaves: Vec<[u8; 32]>) -> Self {
        debug_assert!(leaves.len().is_power_of_two());
        // Perform hash-chaining up to root.
        let mut layers = vec![leaves.clone()];
        let mut size = leaves.len() >> 1;
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
        }
    }

    pub fn root(&self) -> [u8; 32] {
        self.root
    }
}
