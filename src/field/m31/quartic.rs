use crate::field::m31::{Mersenne31Complex, M31};
use core::ops::{Add, Mul, Sub};

pub struct Mersenne31Quartic([Mersenne31Complex; 2]);

impl Add<M31> for Mersenne31Quartic {
    type Output = Self;

    fn add(mut self, rhs: M31) -> Self::Output {
        self.0[0].0[0] += rhs;
        self
    }
}

impl Sub<M31> for Mersenne31Quartic {
    type Output = Self;

    fn sub(mut self, rhs: M31) -> Self::Output {
        self.0[0].0[0] -= rhs;
        self
    }
}

impl Mul<M31> for Mersenne31Quartic {
    type Output = Self;

    fn mul(mut self, rhs: M31) -> Self::Output {
        self.0[0].0[0] *= rhs;
        self.0[0].0[1] *= rhs;
        self.0[1].0[0] *= rhs;
        self.0[1].0[1] *= rhs;
        self
    }
}

impl From<M31> for Mersenne31Quartic {
    fn from(v: M31) -> Self {
        Self([
            Mersenne31Complex([v, M31(0)]),
            Mersenne31Complex([M31(0), M31(0)]),
        ])
    }
}
