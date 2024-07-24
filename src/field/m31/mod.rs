// taken mostly from air compiler

pub mod complex;
pub mod quartic;

use crate::field::*;
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

/// The prime field `F_p` where `p = 2^31 - 1`.
// NOTE: using a 64 bit register shouldn't affect performance on 64-bit processors (which is
// where i assume most CPU provers will run) but should save us small extra costs incurred from
// register truncation and expansion in for example multiplications.
//
// XXX: since we use such a large register, lazy reduction could be really viable but only for
// additions/subtractions (does reduction still get to use the same hack in this case?)
#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
#[repr(transparent)]
pub struct M31(pub u64);

impl M31 {
    pub const ORDER: u64 = (1 << 31) - 1;
    pub const MSBITMASK: u64 = ((u32::MAX as u64) << 32) + (1 << 31);

    pub const fn new(value: u64) -> Self {
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
        let x_low = x & Self::ORDER;
        let x_high = (x >> 31) & Self::ORDER;
        let x_sign = x >> 63;
        let res_wrapped = x_low.wrapping_add(x_high);
        let res_wrapped = res_wrapped - x_sign;
        let msb = res_wrapped & Self::MSBITMASK;
        let mut sum = res_wrapped;
        sum ^= msb;
        let mut res = sum + (msb != 0) as u64;
        if res >= Self::ORDER {
            res -= Self::ORDER;
        }
        M31(res)
    }
}

impl Field for M31 {
    const ZERO: Self = Self(0);
    const ONE: Self = Self(1);

    fn from_usize(value: usize) -> Self {
        let value = value as u64;
        debug_assert!(value < Self::ORDER);
        Self(value)
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    fn inverse(&self) -> Option<Self> {
        //Since the nonzero elements of GF(pn) form a finite group with respect to multiplication,
        // a^p^n−1 = 1 (for a ≠ 0), thus the inverse of a is a^p^n−2.
        if self.is_zero() {
            return None;
        }

        let mut p101 = *self;
        p101.exp_power_of_2(2);
        p101 *= self;

        let mut p1111 = p101;
        p1111.square();
        p1111 *= p101;

        let mut p11111111 = p1111;
        p11111111.exp_power_of_2(4);
        p11111111 *= p1111;

        let mut p111111110000 = p11111111;
        p111111110000.exp_power_of_2(4);

        let mut p111111111111 = p111111110000;
        p111111111111 *= p1111;

        let mut p1111111111111111 = p111111110000;
        p1111111111111111.exp_power_of_2(4);
        p1111111111111111 *= p11111111;

        let mut p1111111111111111111111111111 = p1111111111111111;
        p1111111111111111111111111111.exp_power_of_2(12);
        p1111111111111111111111111111 *= p111111111111;

        let mut p1111111111111111111111111111101 = p1111111111111111111111111111;
        p1111111111111111111111111111101.exp_power_of_2(3);
        p1111111111111111111111111111101 *= p101;
        Some(p1111111111111111111111111111101)
    }

    #[inline(always)]
    fn square(&self) -> Self {
        *self * self
    }

    fn double(&self) -> Self {
        self.mul_2exp_u64(1)
    }

    fn div2(&self) -> Self {
        self.div_2exp_u64(1)
    }
}

impl Neg for M31 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if !self.is_zero() {
            Self(Self::ORDER - self.0)
        } else {
            self
        }
    }
}

impl Add<Self> for M31 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let sum = self.0.wrapping_add(rhs.0);
        // cond select of result based on overflow
        // avoids branching but idk if this is really that efficient
        let of = self.0 >= Self::ORDER;
        let reduced = self.0.wrapping_sub(Self::ORDER);
        let mask = 0u64.wrapping_sub(of as u64);
        Self(sum ^ (mask & (sum ^ reduced)))
    }
}

impl Add<&Self> for M31 {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        self + *rhs
    }
}

impl Sub<Self> for M31 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut res = self.0.wrapping_sub(rhs.0);
        let msb = res & Self::MSBITMASK;
        res ^= msb;
        Self(res - (msb != 0) as u64)
    }
}

impl Sub<&Self> for M31 {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        self - *rhs
    }
}

impl Mul<Self> for M31 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let product = self.0 * rhs.0; // since we're using u64 no need to care about overflow or
                                      // casting
        let product_low = product & Self::ORDER;
        let product_high = product >> 31;
        Self(product_low) + Self(product_high)
    }
}

impl Mul<&Self> for M31 {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        self * *rhs
    }
}

impl AddAssign<Self> for M31 {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl AddAssign<&Self> for M31 {
    fn add_assign(&mut self, other: &Self) {
        *self = *self + *other;
    }
}

impl SubAssign<Self> for M31 {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl SubAssign<&Self> for M31 {
    fn sub_assign(&mut self, other: &Self) {
        *self = *self - *other;
    }
}

impl MulAssign<Self> for M31 {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl MulAssign<&Self> for M31 {
    fn mul_assign(&mut self, other: &Self) {
        *self = *self * *other;
    }
}

impl PrimeField for M31 {
    const TWO: Self = Self(2);
    const MINUS_ONE: Self = Self(Self::ORDER - 1);
    const NUM_BYTES_IN_REPR: usize = 8;
    const CHAR_BITS: usize = 31;
    const CHARACTERISTICS: u64 = Self::ORDER;

    #[inline(always)]
    fn as_u64(self) -> u64 {
        self.0
    }

    #[inline(always)]
    fn from_u64_unchecked(value: u64) -> Self {
        Self::new(value)
    }
    #[inline(always)]
    fn from_u64(value: u64) -> Option<Self> {
        if value >= Self::ORDER {
            None
        } else {
            Some(Self(value))
        }
    }

    #[inline(always)]
    fn from_u64_with_reduction(value: u64) -> Self {
        Self(value % Self::ORDER)
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

    fn to_le_bytes(self) -> [u8; Self::NUM_BYTES_IN_REPR] {
        self.0.to_le_bytes()
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
    M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER))
}
