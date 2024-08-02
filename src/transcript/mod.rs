use crate::field::{ChallengeField, Field, PrimeField};
use blake2::{Blake2s256, Digest};
use core::marker::PhantomData;

pub trait Transcript<F: Field> {
    fn draw_challenge(&mut self) -> F;
    fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E;
    fn draw_bits(&mut self, bits: usize) -> usize;
    fn observe_witness(&mut self, witness: F);
    fn observe_witnesses(&mut self, witness: &[F]);
}

pub struct Blake2sTranscript<P: PrimeField> {
    hasher: Blake2s256,
    _marker: PhantomData<P>,
}

impl<F: Field, P: PrimeField> Transcript<F> for Blake2sTranscript<P>
where
    [(); F::NUM_BYTES_IN_REPR]:,
{
    fn draw_challenge(&mut self) -> F {
        let output = <[u8; 32]>::from(self.hasher.finalize_reset());
        let mut arr = [0u8; 8];
        arr.copy_from_slice(&output[..F::NUM_BYTES_IN_REPR]);
        let n = usize::from_le_bytes(arr);
        F::from_usize(n)
    }

    fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E {
        debug_assert!(F::NUM_BYTES_IN_REPR * E::DEGREE <= 32); // we only extract 256 bits so can't
                                                               // make much more
        let output = <[u8; 32]>::from(self.hasher.finalize_reset());
        let els = output
            .chunks(F::NUM_BYTES_IN_REPR)
            .map(|chunk| {
                let mut arr = [0u8; 8];
                arr.copy_from_slice(&chunk);
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
        let output = <[u8; 32]>::from(self.hasher.finalize_reset());
        let mut v = output[..bytes + 1].to_vec();

        // Flip redundant bits in last byte.
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
        self.hasher.update(witness.to_le_bytes().to_vec());
    }

    fn observe_witnesses(&mut self, witness: &[F]) {
        for w in witness {
            self.hasher.update(w.to_le_bytes().to_vec());
        }
    }
}
