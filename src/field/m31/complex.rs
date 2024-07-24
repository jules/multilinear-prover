use super::M31;
use crate::field::{Field, PrimeField};
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct M31_2 {
    pub c0: M31,
    pub c1: M31,
}

impl M31_2 {
    pub fn mul_by_nonresidue(&mut self) {
        let (a, b) = (self.c0, self.c1);
        let mut real = a;
        real.double();
        real -= b;

        let mut imag = b;
        imag.double();
        imag += a;

        self.c0 = real;
        self.c1 = imag;
    }
}

impl Field for M31_2 {
    const ZERO: Self = Self {
        c0: M31::ZERO,
        c1: M31::ZERO,
    };
    const ONE: Self = Self {
        c0: M31::ONE,
        c1: M31::ZERO,
    };

    fn from_usize(v: usize) -> Self {
        Self {
            c0: M31::new(v as u64),
            c1: M31::ZERO,
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn inverse(&self) -> Option<Self> {
        None
    }

    fn square(&self) -> Self {
        let mut new = *self;
        let mut v0 = new.c0;
        v0 -= new.c1;
        let mut v3 = new.c0;
        let mut t0 = new.c1;
        t0.mul_by_nonresidue();
        v3 -= t0;
        let mut v2 = new.c0;
        v2 *= new.c1;
        v0 *= v3;
        v0 += v2;
        new.c1 = v2;
        new.c1.double();
        new.c0 = v0;
        v2.mul_by_nonresidue();
        new.c0 += v2;
        new
    }

    fn double(&self) -> Self {
        let mut new = *self;
        new.c0.double();
        new.c1.double();
        new
    }

    fn div2(&self) -> Self {
        let mut new = *self;
        new.c0.div2();
        new.c1.div2();
        new
    }
}

impl Neg for M31_2 {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.c0 = self.c0.neg();
        self.c1 = self.c1.neg();
        self
    }
}

impl Add<Self> for M31_2 {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.c0 += rhs.c0;
        self.c1 += rhs.c1;
        self
    }
}

impl Add<&Self> for M31_2 {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        self + *rhs
    }
}

impl Sub<Self> for M31_2 {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self.c0 -= rhs.c0;
        self.c1 -= rhs.c1;
        self
    }
}

impl Sub<&Self> for M31_2 {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        self - *rhs
    }
}

impl Mul<Self> for M31_2 {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        let mut v0 = self.c0;
        v0 *= rhs.c0;
        let mut v1 = self.c1;
        v1 *= rhs.c1;

        let t = self.c0;
        self.c1 += t;

        let mut t0 = rhs.c0;
        t0 += rhs.c1;
        self.c1 *= t0;
        self.c1 += v0;
        self.c1 += v1;
        self.c0 = v0;
        v1 = v1.neg();
        self.c0 += v1;
        self
    }
}

impl Mul<&Self> for M31_2 {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        self * *rhs
    }
}

impl AddAssign<Self> for M31_2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl AddAssign<&Self> for M31_2 {
    fn add_assign(&mut self, rhs: &Self) {
        *self = *self + *rhs
    }
}

impl SubAssign<Self> for M31_2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl SubAssign<&Self> for M31_2 {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = *self - *rhs
    }
}

impl MulAssign<Self> for M31_2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl MulAssign<&Self> for M31_2 {
    fn mul_assign(&mut self, rhs: &Self) {
        *self = *self * *rhs
    }
}

impl Hash for M31_2 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.c0.hash(state);
        self.c1.hash(state);
    }
}

impl Display for M31_2 {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "F2[{}, {}]",
            self.c0.as_reduced_u32(),
            self.c1.as_reduced_u32()
        )
    }
}

pub fn rand_fp2_from_rng<R: rand::Rng>(rng: &mut R) -> M31_2 {
    let a = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    let b = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    M31_2 { c0: a, c1: b }
}
