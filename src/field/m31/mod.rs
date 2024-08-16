pub mod complex;
pub mod quartic;

use crate::field::*;
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
};

/// The prime field `F_p` where `p = 2^31 - 1`.
#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
#[repr(transparent)]
pub struct M31(pub u32);

impl M31 {
    pub const ORDER: u32 = (1 << 31) - 1;
    pub const MSBITMASK: u32 = 1 << 31;

    pub const fn new(value: u32) -> Self {
        debug_assert!(value < Self::ORDER);
        Self(value)
    }

    #[inline(always)]
    pub const fn as_reduced_u32(&self) -> u32 {
        let mut c = self.0;
        if c >= Self::ORDER {
            c -= Self::ORDER;
        }
        c as u32
    }

    #[inline]
    pub fn mul_2exp_u64(&self, exp: u64) -> Self {
        // In a Mersenne field, multiplication by 2^k is just a left rotation by k bits.
        debug_assert!(exp < 31);
        let left = (self.0 << exp) & Self::ORDER;
        let right = self.0 >> (31 - exp);
        Self::new(left | right)
    }

    #[inline]
    pub fn div_2exp_u64(&self, exp: u64) -> Self {
        // In a Mersenne field, division by 2^k is just a right rotation by k bits.
        debug_assert!(exp < 31);
        let left = self.0 >> exp;
        let right = (self.0 << (31 - exp)) & Self::ORDER;
        Self::new(left | right)
    }

    #[inline(always)]
    pub fn from_negative_u64_with_reduction(x: u64) -> Self {
        let x_low = (x as u32) & Self::ORDER;
        let x_high = ((x >> 31) as u32) & Self::ORDER;
        let x_sign = (x >> 63) as u32;
        let res_wrapped = x_low.wrapping_add(x_high);
        let res_wrapped = res_wrapped - x_sign;
        let msb = res_wrapped & Self::MSBITMASK;
        let mut sum = res_wrapped;
        sum ^= msb;
        let mut res = sum + u32::from(msb != 0);
        if res >= Self::ORDER {
            res -= Self::ORDER;
        }
        M31(res)
    }

    #[inline(always)]
    pub fn mul_by_nonresidue(&mut self) {
        self.negate();
    }
}

impl Field for M31 {
    const ZERO: Self = Self(0);
    const ONE: Self = Self(1);
    const NUM_BYTES_IN_REPR: usize = 4;

    fn from_usize(value: usize) -> Self {
        let value = value as u64;
        Self::from_u64_with_reduction(value)
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    // Since the nonzero elements of GF(pn) form a finite group with respect to multiplication,
    // a^p^n−1 = 1 (for a ≠ 0), thus the inverse of a is a^p^n−2.
    // Code from https://github.com/Plonky3/Plonky3/blob/main/mersenne-31/src/mersenne_31.rs#L235
    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let mut p101 = *self;
        p101.exp_power_of_2(2);
        p101.mul_assign(&self);

        let mut p1111 = p101;
        p1111.square();
        p1111.mul_assign(&p101);

        let mut p11111111 = p1111;
        p11111111.exp_power_of_2(4);
        p11111111.mul_assign(&p1111);

        let mut p111111110000 = p11111111;
        p111111110000.exp_power_of_2(4);

        let mut p111111111111 = p111111110000;
        p111111111111.mul_assign(&p1111);

        let mut p1111111111111111 = p111111110000;
        p1111111111111111.exp_power_of_2(4);
        p1111111111111111.mul_assign(&p11111111);

        let mut p1111111111111111111111111111 = p1111111111111111;
        p1111111111111111111111111111.exp_power_of_2(12);
        p1111111111111111111111111111.mul_assign(&p111111111111);

        let mut p1111111111111111111111111111101 = p1111111111111111111111111111;
        p1111111111111111111111111111101.exp_power_of_2(3);
        p1111111111111111111111111111101.mul_assign(&p101);
        Some(p1111111111111111111111111111101)
    }

    fn add_assign(&mut self, other: &Self) {
        let sum = self.0.wrapping_add(other.0);
        // cond select of result based on overflow
        // avoids branching but idk if this is really that efficient
        let of = sum >= Self::ORDER;
        let reduced = sum.wrapping_sub(Self::ORDER);
        let mask = 0u32.wrapping_sub(of as u32);
        self.0 = sum ^ (mask & (sum ^ reduced));
    }

    fn sub_assign(&mut self, other: &Self) {
        let mut sum = self.0.wrapping_sub(other.0);
        let msb = sum & Self::MSBITMASK;
        sum ^= msb;
        self.0 = sum - u32::from(msb != 0);
    }

    fn mul_assign(&mut self, other: &Self) {
        let product = u64::from(self.0) * u64::from(other.0);
        let product_low = (product as u32) & ((1 << 31) - 1);
        let product_high = (product >> 31) as u32;
        *self = Self(product_low);
        self.add_assign(&Self(product_high));
    }

    fn negate(&mut self) {
        if !self.is_zero() {
            self.0 = Self::ORDER - self.0
        }
    }

    #[inline(always)]
    fn square(&mut self) {
        self.mul_assign(&self.clone());
    }

    fn double(&mut self) {
        *self = self.mul_2exp_u64(1);
    }

    fn div2(&mut self) {
        *self = self.div_2exp_u64(1);
    }

    fn to_le_bytes(self) -> [u8; Self::NUM_BYTES_IN_REPR] {
        self.as_reduced_u32().to_le_bytes()
    }
}

impl PrimeField for M31 {
    const TWO: Self = Self(2);
    const MINUS_ONE: Self = Self(Self::ORDER - 1);
    const CHAR_BITS: usize = 31;
    const CHARACTERISTICS: u64 = Self::ORDER as u64;

    #[inline(always)]
    fn as_u64(self) -> u64 {
        self.0 as u64
    }

    #[inline(always)]
    fn from_u64_unchecked(value: u64) -> Self {
        Self::new(value as u32)
    }
    #[inline(always)]
    fn from_u64(value: u64) -> Option<Self> {
        let value = value as u32;
        if value >= Self::ORDER {
            None
        } else {
            Some(Self(value))
        }
    }

    #[inline(always)]
    fn from_u64_with_reduction(value: u64) -> Self {
        Self((value % (Self::ORDER as u64)) as u32)
    }

    #[inline(always)]
    fn as_u64_reduced(&self) -> u64 {
        self.as_reduced_u32() as u64
    }

    fn as_boolean(&self) -> bool {
        assert!(self.0 == 0 || self.0 == 1);
        self.0 != 0
    }

    fn from_boolean(flag: bool) -> Self {
        if flag {
            Self::ONE
        } else {
            Self::ZERO
        }
    }

    fn increment_unchecked(&'_ mut self) {
        self.0 += 1;
    }
}

impl Hash for M31 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u32(self.as_reduced_u32())
    }
}

impl Display for M31 {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}

pub fn rand_fp_from_rng<R: rand::Rng>(rng: &mut R) -> M31 {
    M31(rng.gen_range(0..M31::ORDER))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_no_overflow() {
        let mut a = M31(3);
        let b = M31(5);
        a.add_assign(&b);
        a.add_assign(&b);
        assert_eq!(a, M31(13));
    }

    #[test]
    fn test_add_overflow() {
        let mut a = M31(2147483646);
        let b = M31(2147483646);
        a.add_assign(&b);
        a.add_assign(&b);
        assert_eq!(a, M31(2147483644));
    }

    #[test]
    fn test_sub_no_overflow() {
        let mut a = M31(8);
        let b = M31(5);
        a.sub_assign(&b);
        assert_eq!(a, M31(3));
    }

    #[test]
    fn test_sub_overflow() {
        let mut a = M31(1);
        let b = M31(2);
        a.sub_assign(&b);
        assert_eq!(a, M31(2147483646));
    }

    #[test]
    fn test_mul_no_overflow() {
        let mut a = M31(8);
        let b = M31(5);
        a.mul_assign(&b);
        assert_eq!(a, M31(40));
    }

    #[test]
    fn test_mul_overflow() {
        let mut a = M31(2147483646);
        let b = M31(2);
        a.mul_assign(&b);
        assert_eq!(a, M31(2147483645));
    }

    #[test]
    fn test_inverse() {
        let mut a = M31(2173);
        let inv = a.inverse().unwrap();
        a.mul_assign(&inv);
        assert_eq!(M31::ONE, a);
    }
}
