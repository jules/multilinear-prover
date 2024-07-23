use super::{complex::M31_2, M31};
use core::ops::{Add, Mul, Sub};

pub struct M31_4(pub [M31_2; 2]);

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
