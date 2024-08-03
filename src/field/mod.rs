pub mod m31;

use core::{
    cmp::{Eq, PartialEq},
    default::Default,
    fmt::{Debug, Display},
    hash::Hash,
    marker::{Send, Sync},
};

pub trait Field:
    'static
    + Clone
    + Copy
    + Default
    + Display
    + Debug
    + Hash
    + PartialEq
    + Eq
    + Send
    + Sync
    + Default
    + Sized
{
    const ZERO: Self;
    const ONE: Self;
    type CharField = Self;
    const NUM_BYTES_IN_REPR: usize;

    fn from_usize(v: usize) -> Self;

    fn is_zero(&self) -> bool;
    fn inverse(&self) -> Option<Self>;

    fn add_assign(&mut self, other: &Self);
    fn sub_assign(&mut self, other: &Self);
    fn mul_assign(&mut self, other: &Self);
    fn negate(&mut self);
    fn square(&mut self);
    fn double(&mut self);
    fn div2(&mut self);

    fn pow(&self, mut exp: u32) -> Self {
        let mut base = *self;
        let mut result = Self::ONE;
        while exp > 0 {
            if exp % 2 == 1 {
                result.mul_assign(&base);
            }

            exp >>= 1;
            base.square();
        }

        result
    }

    fn exp_power_of_2(&mut self, power_log: usize) {
        for _ in 0..power_log {
            self.square();
        }
    }

    fn to_le_bytes(self) -> [u8; Self::NUM_BYTES_IN_REPR];
}

pub trait PrimeField: Field {
    const TWO: Self;
    const MINUS_ONE: Self;

    const CHAR_BITS: usize;
    const CHARACTERISTICS: u64;

    fn as_u64(self) -> u64;
    fn from_u64_unchecked(value: u64) -> Self;
    fn from_u64_with_reduction(value: u64) -> Self;
    fn from_u64(value: u64) -> Option<Self>;
    fn as_u64_reduced(&self) -> u64;

    fn as_boolean(&self) -> bool;

    fn from_boolean(flag: bool) -> Self {
        if flag {
            Self::ONE
        } else {
            Self::ZERO
        }
    }

    fn increment_unchecked(&'_ mut self);
}

pub trait ChallengeField<F: Field>: Field + Into<Vec<F>> + From<F> {
    const DEGREE: usize;

    fn new(values: Vec<F>) -> Self;
    fn add_base(&mut self, other: &F);
    fn sub_base(&mut self, other: &F);
    fn mul_base(&mut self, other: &F);
    fn real_coeff(&self) -> F;
}

pub trait TwoAdicField: Field {
    /// The number of factors of two in this field's multiplicative group.
    const TWO_ADICITY: usize;

    /// Returns a generator of the multiplicative group of order `2^bits`.
    /// Assumes `bits < TWO_ADICITY`, otherwise the result is undefined.
    /// all functions here except for two_adic_generator should not even exist
    #[must_use]
    fn two_adic_generator() -> Self;

    #[must_use]
    fn two_adic_group_order() -> usize;
}
