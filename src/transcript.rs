// also taken from air compiler for convenience

use crate::field::{Field, FieldExtension, PrimeField};
use crate::field_utils::AlgebraicRoundFunction;
use crate::utils::*;
use derivative::Derivative;
use unroll::unroll_for_loops;

#[derive(Derivative)]
#[derivative(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AbsorbtionMode {
    Addition = 0,
    Overwrite, // https://keccak.team/files/CSF-0.1.pdf - e.g. some info about this mode
}

pub trait Transcript<F: PrimeField>: Send + Sync + std::fmt::Debug {
    const IS_ALGEBRAIC: bool = true;

    fn new() -> Self;
    fn witness_field_element(&mut self, field_el: F);
    fn get_challenge(&mut self) -> F;

    fn get_multiple_challenges_fixed<const N: usize>(&mut self) -> [F; N] {
        let mut result = [F::ZERO; N];
        for dst in result.iter_mut() {
            *dst = self.get_challenge();
        }

        result
    }

    fn get_multiple_challenges<B: GoodAllocator>(&mut self, num_challenges: usize) -> Vec<F, B> {
        let mut result = Vec::with_capacity_in(num_challenges, B::default());
        for _ in 0..num_challenges {
            let chal = self.get_challenge();
            result.push(chal);
        }

        result
    }

    fn get_challenge_bytes(&mut self, _num_bytes: usize) -> Vec<u8> {
        if Self::IS_ALGEBRAIC {
            unimplemented!("Should not be called on algebraic transcripts")
        } else {
            todo!("Should be implemented if non-algebraic")
        }
    }

    fn get_challenge_in_ext<T: FieldExtension<F>>(&mut self) -> T
    where
        [(); T::DEGREE]:,
    {
        let inner = self.get_multiple_challenges_fixed::<{ T::DEGREE }>();
        T::from_coeffs_in_base(&inner)
    }

    fn get_multiple_challenges_in_ext_fixed<T: FieldExtension<F> + Field, const N: usize>(
        &mut self,
    ) -> [T; N]
    where
        [(); T::DEGREE]:,
    {
        let mut result = [T::ZERO; N];
        for dst in result.iter_mut() {
            *dst = self.get_challenge_in_ext();
        }

        result
    }

    fn get_multiple_challenges_in_ext<T: FieldExtension<F>, B: GoodAllocator>(
        &mut self,
        num_challenges: usize,
    ) -> Vec<T, B>
    where
        [(); T::DEGREE]:,
    {
        let mut result = Vec::with_capacity_in(num_challenges, B::default());
        for _ in 0..num_challenges {
            let chal = self.get_challenge_in_ext();
            result.push(chal);
        }

        result
    }

    fn witness_field_elements(&mut self, field_els: &[F]) {
        field_els
            .iter()
            .for_each(|current| self.witness_field_element(*current))
    }

    fn witness_field_element_in_ext<T: FieldExtension<F>>(&mut self, field_els: T)
    where
        [(); T::DEGREE]:,
    {
        self.witness_field_elements(field_els.coeffs_in_base())
    }

    fn witness_field_elements_in_ext<T: FieldExtension<F>>(&mut self, field_els: &[T])
    where
        [(); T::DEGREE]:,
    {
        for elem in field_els.iter() {
            self.witness_field_elements(elem.coeffs_in_base())
        }
    }

    fn witness_merkle_tree_cap<H: CompressionFunction<F>>(&mut self, cap: &[H::Output]) {
        cap.iter()
            .for_each(|root| self.witness_field_elements(root.as_ref()))
    }
}

// algebraic_transcript = algebraic_round_function + defined rate and capacity + buffers for already allocated challenges and elements to commit
// e.g. consider the following situation: Rate = 2 but we iteratively squeeze two challenges - so there is actually no need to roll round function in between
// built with coherence with the SAFE API as defined in https://eprint.iacr.org/2023/522
#[derive(Derivative)]
#[derivative(Debug)]
pub struct AlgebraicSpongeBasedTranscript<
    F: PrimeField,
    const RATE: usize,
    const CAPACITY: usize,
    R: AlgebraicRoundFunction<F, { RATE + CAPACITY }>,
> {
    input_buffer: [F; RATE],
    filled: usize,
    squeezed: usize,
    state: [F; RATE + CAPACITY],
    #[derivative(Debug = "ignore")]
    round_function: R,
    mode: AbsorbtionMode,
}

impl<
        F: PrimeField,
        const RATE: usize,
        const CAPACITY: usize,
        R: AlgebraicRoundFunction<F, { RATE + CAPACITY }>,
    > Transcript<F> for AlgebraicSpongeBasedTranscript<F, RATE, CAPACITY, R>
where
    R: Default,
{
    fn new() -> Self {
        Self {
            input_buffer: [F::ZERO; RATE],
            round_function: R::default(),
            filled: 0,
            squeezed: 0,
            state: [F::ZERO; { RATE + CAPACITY }],
            mode: AbsorbtionMode::Addition,
        }
    }

    fn witness_field_element(&mut self, field_el: F) {
        self.witness_field_elements(&[field_el])
    }

    #[unroll_for_loops]
    fn witness_field_elements(&mut self, values: &[F]) {
        let mut rest = values;
        if rest.len() <= RATE - self.filled {
            self.input_buffer[self.filled..self.filled + rest.len()].copy_from_slice(values);
            self.filled += values.len();
            return;
        }

        // quickly fill the unfilled, then do full absorbtion rounds, then save the rest
        if rest.len() > RATE - self.filled {
            let (to_use, other) = rest.split_at(RATE - self.filled);
            rest = other;
            for (dst, src) in self.state[..RATE]
                .iter_mut()
                .zip(self.input_buffer[..self.filled].iter().chain(to_use.iter()))
            {
                match self.mode {
                    AbsorbtionMode::Addition => {
                        dst.add_assign(src);
                    }
                    AbsorbtionMode::Overwrite => {
                        *dst = *src;
                    }
                }
            }
            self.round_function.round_function(&mut self.state);
            self.filled = 0;
        }

        assert!(self.filled == 0);
        let mut it = rest.chunks_exact(RATE);
        for chunk in &mut it {
            for (dst, src) in self.state[..RATE].iter_mut().zip(chunk.iter()) {
                match self.mode {
                    AbsorbtionMode::Addition => {
                        dst.add_assign(src);
                    }
                    AbsorbtionMode::Overwrite => {
                        *dst = *src;
                    }
                }
            }

            self.round_function.round_function(&mut self.state);
        }

        let tail = it.remainder();
        self.filled += tail.len();
        for (dst, src) in self.input_buffer.iter_mut().zip(tail.iter()) {
            *dst = *src; // save only for a future
        }

        // TODO: check this place in Boojum - probably a bug if we forgot to clear vec of available challenges
        self.squeezed = RATE;
    }

    fn get_challenge(&mut self) -> F {
        if self.squeezed == RATE {
            self.round_function.round_function(&mut self.state);
            self.filled = 0;
            self.squeezed = 0;
        }

        let elem_to_squeeze = self.state[self.squeezed];
        self.squeezed += 1;
        elem_to_squeeze
    }
}

#[derive(Derivative)]
#[derivative(Clone, Debug)]
pub struct Blake2sTranscript {
    inner: blake2::Blake2s256,
    buffer: Vec<u8>,
    available_challenge_bytes: Vec<u8>,
}

use blake2::Digest;

use super::merkle_tree::CompressionFunction;

impl<F: PrimeField> Transcript<F> for Blake2sTranscript {
    const IS_ALGEBRAIC: bool = false;

    fn new() -> Self {
        Self {
            inner: blake2::Blake2s256::new(),
            buffer: Vec::with_capacity(64),
            available_challenge_bytes: Vec::with_capacity(32),
        }
    }

    fn witness_field_element(&mut self, field_el: F) {
        self.available_challenge_bytes.clear();
        self.buffer.extend(field_el.as_u64_reduced().to_le_bytes());
    }

    fn witness_field_elements(&mut self, field_els: &[F]) {
        self.available_challenge_bytes.clear();
        for el in field_els.iter() {
            self.buffer.extend(el.as_u64_reduced().to_le_bytes());
        }
    }

    fn get_challenge(&mut self) -> F {
        if self.buffer.is_empty() == false {
            self.inner.update(&self.buffer);
            self.buffer.clear();
            self.available_challenge_bytes.clear();

            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);
        }

        if self.available_challenge_bytes.is_empty() == false {
            assert!(self.available_challenge_bytes.len() % 8 == 0);
            let mut buffer = [0u8; 8];
            for (dst, src) in buffer
                .iter_mut()
                .zip(self.available_challenge_bytes.drain(..8))
            {
                *dst = src;
            }
            let as_u64 = u64::from_le_bytes(buffer);
            let challenge = F::from_u64_with_reduction(as_u64);

            challenge
        } else {
            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);

            assert!(self.available_challenge_bytes.is_empty() == false);
            Transcript::<F>::get_challenge(self)
        }
    }

    fn get_challenge_bytes(&mut self, num_bytes: usize) -> Vec<u8> {
        if self.buffer.is_empty() == false {
            self.inner.update(&self.buffer);
            self.buffer.clear();
            self.available_challenge_bytes.clear();

            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);
        }

        if self.available_challenge_bytes.len() >= num_bytes {
            let result: Vec<u8> = self.available_challenge_bytes.drain(..num_bytes).collect();

            result
        } else {
            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);

            assert!(self.available_challenge_bytes.is_empty() == false);
            Transcript::<F>::get_challenge_bytes(self, num_bytes)
        }
    }
}

#[derive(Derivative)]
#[derivative(Clone, Debug)]
pub struct Keccak256Transcript {
    inner: sha3::Keccak256,
    buffer: Vec<u8>,
    available_challenge_bytes: Vec<u8>,
}

impl<F: PrimeField> Transcript<F> for Keccak256Transcript {
    const IS_ALGEBRAIC: bool = false;

    fn new() -> Self {
        Self {
            inner: sha3::Keccak256::new(),
            buffer: Vec::with_capacity(64),
            available_challenge_bytes: Vec::with_capacity(32),
        }
    }

    fn witness_field_element(&mut self, field_el: F) {
        self.available_challenge_bytes.clear();
        self.buffer.extend(field_el.as_u64_reduced().to_le_bytes());
    }

    fn witness_field_elements(&mut self, field_els: &[F]) {
        self.available_challenge_bytes.clear();
        for el in field_els.iter() {
            self.buffer.extend(el.as_u64_reduced().to_le_bytes());
        }
    }

    fn get_challenge(&mut self) -> F {
        if self.buffer.is_empty() == false {
            self.inner.update(&self.buffer);
            self.buffer.clear();
            self.available_challenge_bytes.clear();

            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);
        }

        if self.available_challenge_bytes.is_empty() == false {
            assert!(self.available_challenge_bytes.len() % 8 == 0);
            let mut buffer = [0u8; 8];
            for (dst, src) in buffer
                .iter_mut()
                .zip(self.available_challenge_bytes.drain(..8))
            {
                *dst = src;
            }
            let as_u64 = u64::from_le_bytes(buffer);
            let challenge = F::from_u64_with_reduction(as_u64);

            challenge
        } else {
            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);

            assert!(self.available_challenge_bytes.is_empty() == false);
            Transcript::<F>::get_challenge(self)
        }
    }

    fn get_challenge_bytes(&mut self, num_bytes: usize) -> Vec<u8> {
        if self.buffer.is_empty() == false {
            self.inner.update(&self.buffer);
            self.buffer.clear();
            self.available_challenge_bytes.clear();

            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);
        }

        if self.available_challenge_bytes.len() >= num_bytes {
            let result: Vec<u8> = self.available_challenge_bytes.drain(..num_bytes).collect();

            result
        } else {
            // reseed
            let mut output = [0u8; 32];
            let raw_output = self.inner.finalize_reset();
            output[..].copy_from_slice(raw_output.as_slice());

            self.inner.update(output);
            self.available_challenge_bytes.extend(output);

            assert!(self.available_challenge_bytes.is_empty() == false);
            Transcript::<F>::get_challenge_bytes(self, num_bytes)
        }
    }
}

pub struct BoolsBuffer {
    pub available: Vec<bool>,
    pub max_needed: usize,
}

impl BoolsBuffer {
    pub fn new(max_needed_bits: usize) -> Self {
        BoolsBuffer {
            available: vec![],
            max_needed: max_needed_bits,
        }
    }

    pub fn get_bits<F: PrimeField, T: Transcript<F>>(
        &mut self,
        transcript: &mut T,
        num_bits: usize,
    ) -> Vec<bool> {
        if self.available.len() >= num_bits {
            let give: Vec<_> = self.available.drain(..num_bits).collect();

            give
        } else {
            if T::IS_ALGEBRAIC {
                // we use some heuristics to only take part of the bits from field
                // element to get better uniformity of query indexes
                let bits_avaiable = F::CHAR_BITS - self.max_needed;

                // get 1 field element from transcript
                let field_el = transcript.get_challenge();
                let field_el_u64 = field_el.as_u64_reduced();
                let t = [field_el_u64];
                let mut lsb_iterator = LSBIterator::new(&t);
                for _ in 0..bits_avaiable {
                    let bit = lsb_iterator.next().unwrap();
                    self.available.push(bit);
                }
            } else {
                // we assume that it produced BYTES and those are uniform
                let bytes: [u8; 8] = transcript
                    .get_challenge_bytes(8)
                    .try_into()
                    .expect("length must match");
                let as_u64 = u64::from_le_bytes(bytes);
                let t = [as_u64];
                let mut lsb_iterator = LSBIterator::new(&t);
                for _ in 0..64 {
                    let bit = lsb_iterator.next().unwrap();
                    self.available.push(bit);
                }
            }

            self.get_bits(transcript, num_bits)
        }
    }
}
