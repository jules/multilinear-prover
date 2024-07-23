use super::M31;
use crate::field::PrimeField;
use core::ops::Neg;

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
}

pub fn rand_fp2_from_rng<R: rand::Rng>(rng: &mut R) -> M31_2 {
    let a = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    let b = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    M31_2::new(a, b)
}
