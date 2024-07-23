pub mod m31;

use core::{
    cmp::{Eq, PartialEq},
    default::Default,
    fmt::{Debug, Display},
    hash::Hash,
    marker::{Send, Sync},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
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
    + Neg
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + Sized
{
    const ZERO: Self;
    const ONE: Self;
    type CharField = Self;

    fn from_usize(v: usize) -> Self;

    // zero check
    fn is_zero(&self) -> bool;
    fn inverse(&self) -> Option<Self>;

    fn square(&self) -> Self;
    fn double(&self) -> Self;
    fn div2(&self) -> Self;

    fn pow(&self, mut exp: u32) -> Self {
        let mut base = *self;
        let mut result = Self::ONE;
        while exp > 0 {
            if exp % 2 == 1 {
                result.mul_assign(&base);
            }

            exp >>= 1;
            base = base.square();
        }

        result
    }

    fn exp_power_of_2(&mut self, power_log: usize) {
        for _ in 0..power_log {
            self.square();
        }
    }
}

pub trait PrimeField: Field {
    const TWO: Self;
    const MINUS_ONE: Self;
    const NUM_BYTES_IN_REPR: usize;

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

    fn to_le_bytes(self) -> [u8; Self::NUM_BYTES_IN_REPR];

    fn increment_unchecked(&'_ mut self);
}

pub trait ChallengeField<F: Field>:
    Field
    + Add<F, Output = Self>
    + Sub<F, Output = Self>
    + Mul<F, Output = Self>
    + Into<Vec<F>>
    + From<F>
{
    const DEGREE: usize;

    fn new(values: [F; Self::DEGREE]) -> Self;
}
