use blake2::{Blake2s256, Digest};
use rayon::prelude::*;

pub struct MerkleTree {
    root: [u8; 32],
    pub elements: Vec<Vec<[u8; 32]>>, // bottom layer first
}

impl MerkleTree {
    /// Constructs a new MerkleTree from matrices of elements. Each row is hashed and then 2-arity
    /// hashing up to the root is performed.
    pub fn new(mut leaves: Vec<Vec<[u8; 32]>>) -> Self {
        debug_assert!(leaves[0].len().is_power_of_two());
        debug_assert!(leaves.iter().skip(1).all(|a| a.len() == leaves[0].len()));

        // Create initial leaves. If we only have one layer, do nothing. Otherwise, we hash the
        // i'th entry of each layer together to make a single leaf.
        let compressed_leaves = if leaves.len() > 1 {
            (0..leaves[0].len())
                .into_par_iter()
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
                .par_chunks(2)
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

    pub fn get_proof(&self, mut index: usize) -> Vec<[[u8; 32]; 2]> {
        let mut path = Vec::with_capacity(self.elements.len() - 1);
        let mut base_layer = [[0u8; 32]; 2];
        base_layer[index & 1] = self.elements[0][index]; // leaf
        base_layer[(index ^ 1) & 1] = self.elements[0][index ^ 1]; // neighboring element
        path.push(base_layer);
        for i in 1..(self.elements.len() - 1) {
            let mut layer = [[0u8; 32]; 2];
            index >>= 1;
            layer[index & 1] = self.elements[i][index];
            layer[(index ^ 1) & 1] = self.elements[i][index ^ 1];
            path.push(layer);
        }

        path
    }

    pub fn root(&self) -> [u8; 32] {
        self.root
    }

    pub fn verify_path(&self, path: &[[[u8; 32]; 2]]) -> bool {
        // Ensure correct authentication path.
        for i in 0..(path.len() - 1) {
            let mut hasher = Blake2s256::new();
            path[i].iter().for_each(|e| hasher.update(e));
            let hash = <[u8; 32]>::from(hasher.finalize());
            if path[i + 1].iter().all(|e| *e != hash) {
                return false;
            }
        }

        // Check if this authentication path belongs to this tree.
        let mut hasher = Blake2s256::new();
        path[path.len() - 1].iter().for_each(|e| hasher.update(e));
        let hash = <[u8; 32]>::from(hasher.finalize());
        hash == self.root
    }

    pub fn to_hashes(&self) -> Vec<&[u8; 32]> {
        self.elements.iter().flatten().collect::<Vec<&[u8; 32]>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    const TREE_SIZE: usize = 1024;

    fn make_tree() -> MerkleTree {
        let leaves = (0..TREE_SIZE)
            .map(|_| rand::thread_rng().gen::<[u8; 32]>())
            .collect::<Vec<[u8; 32]>>();
        MerkleTree::new(vec![leaves])
    }

    #[test]
    fn test_merkle_proof() {
        let tree = make_tree();

        let index = rand::thread_rng().gen_range(0..TREE_SIZE);
        let proof = tree.get_proof(index);
        assert!(tree.verify_path(&proof));
    }

    #[test]
    fn test_merkle_proof_wrong_tree() {
        let tree = make_tree();

        let index = rand::thread_rng().gen_range(0..TREE_SIZE);
        let proof = tree.get_proof(index);

        let tree = make_tree();
        assert!(!tree.verify_path(&proof));
    }
}
