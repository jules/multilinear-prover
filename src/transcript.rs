use crate::field::{ChallengeField, Field};
use blake2::{Blake2s256, Digest};
use core::marker::PhantomData;

pub trait Transcript<F: Field>: Default + Send + Sync {
    fn draw_challenge(&mut self) -> F;
    fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E;
    fn draw_bits(&mut self, bits: usize) -> usize;
    fn observe_witness(&mut self, witness: F);
    fn observe_witnesses(&mut self, witness: &[F]);
    fn observe_hashes(&mut self, hashes: &[&[u8; 32]]);
}

pub trait IntoObservable {
    fn into_observable(&self) -> Vec<&[u8; 32]>;
}

#[derive(Default)]
pub struct Blake2sTranscript<F: Field> {
    hasher: Blake2s256,
    state: Vec<u8>,
    _marker: PhantomData<F>,
}

impl<F: Field> Blake2sTranscript<F> {
    // We wrap this into a convenient method since we need to do some state trickery.
    fn output_hash(&mut self) -> [u8; 32] {
        self.hasher.update(self.state.clone());
        let output = <[u8; 32]>::from(self.hasher.finalize_reset());
        self.state.extend(output.to_vec());
        output
    }
}

impl<F: Field> Transcript<F> for Blake2sTranscript<F>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    fn draw_challenge(&mut self) -> F {
        let output = self.output_hash();
        let mut arr = [0u8; 8];
        arr[..F::NUM_BYTES_IN_REPR].copy_from_slice(&output[..F::NUM_BYTES_IN_REPR]);
        let n = usize::from_le_bytes(arr);
        F::from_usize(n)
    }

    fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E {
        debug_assert!(F::NUM_BYTES_IN_REPR * E::DEGREE <= 32); // we only extract 256 bits so can't
                                                               // make much more
        let output = self.output_hash();
        let els = output
            .chunks(F::NUM_BYTES_IN_REPR)
            .take(E::DEGREE)
            .map(|chunk| {
                let mut arr = [0u8; 8];
                arr[..F::NUM_BYTES_IN_REPR].copy_from_slice(&chunk);
                let n = usize::from_le_bytes(arr);
                F::from_usize(n)
            })
            .collect::<Vec<F>>();

        E::new(els)
    }

    fn draw_bits(&mut self, bits: usize) -> usize {
        debug_assert!(bits <= 64); // need to fit in a usize
        let bytes = bits / 8;
        let leftover = bits % 8;
        let output = self.output_hash();
        let mut v = output[..bytes + 1].to_vec();

        // Flip superfluous bits in last byte.
        v[bytes] &= (0..leftover).fold(0u8, |mut acc, i| {
            if i != 0 {
                acc <<= 1;
            }

            acc + 1
        });

        // Now fill out the array and cast to usize.
        let mut arr = [0u8; 8];
        arr[..v.len()].copy_from_slice(&v);
        usize::from_le_bytes(arr)
    }

    fn observe_witness(&mut self, witness: F) {
        let bytes = witness.to_le_bytes().to_vec();
        self.state.extend(bytes.clone());
    }

    fn observe_witnesses(&mut self, witness: &[F]) {
        self.state.reserve(witness.len() * F::NUM_BYTES_IN_REPR);
        for w in witness {
            let bytes = w.to_le_bytes().to_vec();
            self.state.extend(bytes.clone());
        }
    }

    fn observe_hashes(&mut self, hashes: &[&[u8; 32]]) {
        self.state.reserve(hashes.len() * 32);
        for h in hashes {
            self.state.extend(h.clone());
        }
    }
}
