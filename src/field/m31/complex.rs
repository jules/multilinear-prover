use super::M31;
use crate::field::{Field, PrimeField};
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
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
        real.sub_assign(&b);

        let mut imag = b;
        imag.double();
        imag.add_assign(&a);

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
    const NUM_BYTES_IN_REPR: usize = M31::NUM_BYTES_IN_REPR * 2;

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
        let mut v0 = self.c0;
        v0.square();
        let mut v1 = self.c1;
        v1.square();
        // v0 = v0 - beta * v1
        let mut v1_by_nonresidue = v1;
        v1_by_nonresidue.mul_by_nonresidue();
        v0.sub_assign(&v1_by_nonresidue);
        match v0.inverse() {
            Some(inversed) => {
                let mut c0 = self.c0;
                c0.mul_assign(&inversed);
                let mut c1 = self.c1;
                c1.mul_assign(&inversed);
                c1.negate();

                let new = Self { c0, c1 };
                Some(new)
            }
            None => None,
        }
    }

    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }

    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }

    fn mul_assign(&mut self, other: &Self) {
        let mut v0 = self.c0;
        v0.mul_assign(&other.c0);
        let mut v1 = self.c1;
        v1.mul_assign(&other.c1);

        let t = self.c0;
        self.c1.add_assign(&t);

        let mut t0 = other.c0;
        t0.add_assign(&other.c1);
        self.c1.mul_assign(&t0);
        self.c1.sub_assign(&v0);
        self.c1.sub_assign(&v1);
        self.c0 = v0;
        v1.mul_by_nonresidue();
        self.c0.add_assign(&v1);
    }

    fn negate(&mut self) {
        self.c0.negate();
        self.c1.negate();
    }

    fn square(&mut self) {
        let mut v0 = self.c0;
        v0.sub_assign(&self.c1);
        let mut v3 = self.c0;
        let mut t0 = self.c1;
        t0.mul_by_nonresidue();
        v3.sub_assign(&t0);
        let mut v2 = self.c0;
        v2.mul_assign(&self.c1);
        v0.mul_assign(&v3);
        v0.add_assign(&v2);

        self.c1 = v2;
        self.c1.double();
        self.c0 = v0;
        v2.mul_by_nonresidue();
        self.c0.add_assign(&v2);
    }

    fn double(&mut self) {
        self.c0.double();
        self.c1.double();
    }

    fn div2(&mut self) {
        self.c0.div2();
        self.c1.div2();
    }

    fn to_le_bytes(self) -> [u8; Self::NUM_BYTES_IN_REPR] {
        self.c0
            .as_reduced_u32()
            .to_le_bytes()
            .into_iter()
            .chain(self.c1.as_reduced_u32().to_le_bytes().into_iter())
            .collect::<Vec<u8>>()
            .try_into()
            .expect("should be able to create bytes array")
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
