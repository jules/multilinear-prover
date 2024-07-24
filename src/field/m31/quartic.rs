use super::{complex::M31_2, M31};
use core::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

#[derive(Copy, Clone, Default)]
pub struct M31_4(pub [M31_2; 2]);

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
        M31_2::mul_by_nonresidue(&mut v1);
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
