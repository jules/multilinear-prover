use super::{complex::M31_2, M31};
use crate::field::Field;
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct M31_4(pub [M31_2; 2]);

impl Field for M31_4 {
    const ZERO: Self = M31_4([M31_2::ZERO, M31_2::ZERO]);
    const ONE: Self = M31_4([M31_2::ONE, M31_2::ZERO]);

    fn from_usize(v: usize) -> Self {
        M31_4([M31_2([M31::new(v as u64), M31::ZERO]), M31_2::ZERO])
    }

    fn is_zero(&self) -> bool {
        self.0[0].is_zero() && self.0[1].is_zero()
    }

    fn inverse(&self) -> Option<Self> {
        None
    }

    fn square(&self) -> Self {
        let mut v0 = self.0[0];
        v0 -= self.0[1];
        let mut v3 = self.0[0];
        let mut t0 = self.0[1];
        M31::mul_by_nonresidue(&mut t0);
        v3 -= t0;
        let mut v2 = self.0[0];
        v2 *= self.0[1];
        v0 *= v3;
        v0 += v2;
        self.0[1] = v2;
        self.0[1].double();
        self.0[0] = v0;
        M31::mul_by_nonresidue(&mut v2);
        self.0[0] += v2;
        self
    }

    fn double(&self) -> Self {
        self.0[0].double();
        self.0[1].double();
    }

    fn div2(&self) -> Self {}
}

impl Neg for M31_4 {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.0[0] = self.0[0].neg();
        self.0[1] = self.0[1].neg();
        self
    }
}

impl Add<Self> for M31_4 {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.0[0] += rhs.0[0];
        self.0[1] += rhs.0[1];
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
        self.0[0] -= rhs.0[0];
        self.0[1] -= rhs.0[1];
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
        let mut v0 = self.0[0];
        v0 *= rhs.0[0];
        let mut v1 = self.0[1];
        v1 *= rhs.0[1];

        let t = self.0[0];
        self.0[1] += t;

        let mut t0 = rhs.0[0];
        t0 += rhs.0[1];
        self.0[1] *= t0;
        self.0[1] += v0;
        self.0[1] += v1;
        self.0[0] = v0;
        v1.mul_by_nonresidue();
        self.0[0] += v1;
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
        self.0[0].0[0] += rhs;
        self
    }
}

impl Sub<M31> for M31_4 {
    type Output = Self;

    fn sub(mut self, rhs: M31) -> Self::Output {
        self.0[0].0[0] -= rhs;
        self
    }
}

impl Mul<M31> for M31_4 {
    type Output = Self;

    fn mul(mut self, rhs: M31) -> Self::Output {
        self.0[0].0[0] *= rhs;
        self.0[0].0[1] *= rhs;
        self.0[1].0[0] *= rhs;
        self.0[1].0[1] *= rhs;
        self
    }
}

impl From<M31> for M31_4 {
    fn from(v: M31) -> Self {
        Self([M31_2([v, M31(0)]), M31_2([M31(0), M31(0)])])
    }
}

impl Hash for M31_4 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0[0].hash(state);
        self.0[1].hash(state);
    }
}

impl Display for M31_4 {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "F4[{}, {}, {}, {}]",
            self.0[0].0[0].as_reduced_u32(),
            self.0[0].0[1].as_reduced_u32(),
            self.0[1].0[0].as_reduced_u32(),
            self.0[1].0[1].as_reduced_u32(),
        )
    }
}
