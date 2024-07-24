use super::M31;
use crate::field::{Field, PrimeField};
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct M31_2(pub [M31; 2]);

impl M31_2 {
    pub const fn new(real: M31, imag: M31) -> Self {
        Self([real, imag])
    }

    pub fn real_part(&self) -> M31 {
        self.0[0]
    }

    pub fn imag_part(&self) -> M31 {
        self.0[1]
    }

    pub fn conjugate(&'_ mut self) -> &'_ mut Self {
        self.0[1] = self.0[1].neg();
        self
    }

    pub fn div_2exp_u64(&self, exp: u64) -> Self {
        Self::new(
            self.real_part().div_2exp_u64(exp),
            self.imag_part().div_2exp_u64(exp),
        )
    }

    pub fn mul_by_nonresidue(&mut self) {
        let (a, b) = (self.0[0], self.0[1]);
        let mut real = a;
        real.double();
        real -= b;

        let mut imag = b;
        imag.double();
        imag += a;

        self.0[0] = real;
        self.0[1] = imag;
    }
}

impl Field for M31_2 {
    const ZERO: Self = M31_2([M31::ZERO, M31::ZERO]);
    const ONE: Self = M31_2([M31::ONE, M31::ZERO]);

    fn from_usize(v: usize) -> Self {
        M31_2([M31::new(v as u64), M31::ZERO])
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

impl Neg for M31_2 {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.0[0] = self.0[0].neg();
        self.0[1] = self.0[1].neg();
        self
    }
}

impl Add<Self> for M31_2 {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.0[0] += rhs.0[0];
        self.0[1] += rhs.0[1];
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
        self.0[0] -= rhs.0[0];
        self.0[1] -= rhs.0[1];
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
        v1 = v1.neg();
        self.0[0] += v1;
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
        self.0[0].hash(state);
        self.0[1].hash(state);
    }
}

impl Display for M31_2 {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "F2[{}, {}]",
            self.0[0].as_reduced_u32(),
            self.0[1].as_reduced_u32()
        )
    }
}

pub fn rand_fp2_from_rng<R: rand::Rng>(rng: &mut R) -> M31_2 {
    let a = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    let b = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    M31_2::new(a, b)
}
