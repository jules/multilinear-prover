use super::{complex::M31_2, M31};
use crate::field::{ChallengeField, Field};
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

// TODO this can be derived or done in a macro as we just have two degree-2 extensions in sequence.

#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct M31_4 {
    pub c0: M31_2,
    pub c1: M31_2,
}

impl Field for M31_4 {
    const ZERO: Self = Self {
        c0: M31_2::ZERO,
        c1: M31_2::ZERO,
    };
    const ONE: Self = Self {
        c0: M31_2::ONE,
        c1: M31_2::ZERO,
    };

    fn from_usize(v: usize) -> Self {
        Self {
            c0: M31_2 {
                c0: M31::new(v as u64),
                c1: M31::ZERO,
            },
            c1: M31_2::ZERO,
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

impl ChallengeField<M31> for M31_4 {
    const DEGREE: usize = 4;

    fn new(values: [M31; Self::DEGREE]) -> Self {
        Self {
            c0: M31_2 {
                c0: values[0],
                c1: values[1],
            },
            c1: M31_2 {
                c0: values[2],
                c1: values[3],
            },
        }
    }
}

impl Into<Vec<M31>> for M31_4 {
    fn into(self) -> Vec<M31> {
        vec![self.c0.c0, self.c0.c1, self.c1.c0, self.c1.c1]
    }
}

impl Neg for M31_4 {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.c0 = self.c0.neg();
        self.c1 = self.c1.neg();
        self
    }
}

impl Add<Self> for M31_4 {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.c0 += rhs.c0;
        self.c1 += rhs.c1;
        self
    }
}

impl Add<&Self> for M31_4 {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        self + *rhs
    }
}

impl Sub<Self> for M31_4 {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self.c0 -= rhs.c0;
        self.c1 -= rhs.c1;
        self
    }
}

impl Sub<&Self> for M31_4 {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        self - *rhs
    }
}

impl Mul<Self> for M31_4 {
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
        v1.mul_by_nonresidue();
        self.c0 += v1;
        self
    }
}

impl Mul<&Self> for M31_4 {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        self * *rhs
    }
}

impl AddAssign<Self> for M31_4 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl AddAssign<&Self> for M31_4 {
    fn add_assign(&mut self, rhs: &Self) {
        *self = *self + *rhs
    }
}

impl SubAssign<Self> for M31_4 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl SubAssign<&Self> for M31_4 {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = *self - *rhs
    }
}

impl MulAssign<Self> for M31_4 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl MulAssign<&Self> for M31_4 {
    fn mul_assign(&mut self, rhs: &Self) {
        *self = *self * *rhs
    }
}

impl Add<M31> for M31_4 {
    type Output = Self;

    fn add(mut self, rhs: M31) -> Self::Output {
        self.c0.c0 += rhs;
        self
    }
}

impl Sub<M31> for M31_4 {
    type Output = Self;

    fn sub(mut self, rhs: M31) -> Self::Output {
        self.c0.c0 -= rhs;
        self
    }
}

impl Mul<M31> for M31_4 {
    type Output = Self;

    fn mul(mut self, rhs: M31) -> Self::Output {
        self.c0.c0 *= rhs;
        self.c0.c1 *= rhs;
        self.c1.c0 *= rhs;
        self.c1.c1 *= rhs;
        self
    }
}

impl From<M31> for M31_4 {
    fn from(v: M31) -> Self {
        Self {
            c0: M31_2 {
                c0: v,
                c1: M31::ZERO,
            },
            c1: M31_2::ZERO,
        }
    }
}

impl Hash for M31_4 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.c0.hash(state);
        self.c1.hash(state);
    }
}

impl Display for M31_4 {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "F4[{}, {}, {}, {}]",
            self.c0.c0.as_reduced_u32(),
            self.c0.c1.as_reduced_u32(),
            self.c1.c0.as_reduced_u32(),
            self.c1.c1.as_reduced_u32(),
        )
    }
}
