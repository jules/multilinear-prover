// taken mostly from air compiler

use crate::field::*;
use core::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[derive(Clone, Copy, Default)]
#[repr(transparent)]
/// The prime field `F_p` where `p = 2^31 - 1`.
// NOTE: using a 64 bit register shouldn't affect performance on 64-bit processors (which is
// where i assume most CPU provers will run) but should save us small extra costs incurred from
// register truncation and expansion in for example multiplications.
//
// XXX: since we use such a large register, lazy reduction could be really viable but only for
// additions/subtractions (does reduction still get to use the same hack in this case?)
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

/*
impl BaseField for M31 {
    const QUADRATIC_NON_RESIDUE: M31 = M31::MINUS_ONE;

    fn mul_by_non_residue(elem: &mut Self) {
        elem.negate();
    }
}

pub type Mersenne31Complex = ExtensionField<M31, 2>;

impl Mersenne31Complex {
    pub const fn new(real: M31, imag: M31) -> Self {
        Self {
            coeffs: [real, imag],
        }
    }

    pub fn real_part(&self) -> M31 {
        self.coeffs[0]
    }

    pub fn imag_part(&self) -> M31 {
        self.coeffs[1]
    }

    pub fn conjugate(&'_ mut self) -> &'_ mut Self {
        self.coeffs[1].negate();
        self
    }

    pub fn div_2exp_u64(&self, exp: u64) -> Self {
        Self::new(
            self.real_part().div_2exp_u64(exp),
            self.imag_part().div_2exp_u64(exp),
        )
    }
}

impl BaseField for Mersenne31Complex {
    // 2 + i is non-residue
    const QUADRATIC_NON_RESIDUE: Mersenne31Complex = Mersenne31Complex {
        coeffs: [M31::TWO, M31::ONE],
    };

    fn mul_by_non_residue(elem: &mut Self) {
        // (a + b * i)(2 + i) = (2 * a - b) + (2 * b + a)i
        let [a, b] = elem.coeffs;
        let mut real = a;
        real.double();
        real.sub_assign(&b);

        let mut imag = b;
        imag.double();
        imag.add_assign(&a);

        elem.coeffs = [real, imag];
    }
}
pub type Mersenne31Quartic = ExtensionField<Mersenne31Complex, 2>;

impl TwoAdicField for Mersenne31Complex {
    const TWO_ADICITY: usize = 31;

    fn two_adic_generator() -> Self {
        // element of order p+1 - generator of cicrcle group
        Self::from_coeffs_in_base(&[
            M31::new(311014874),
            M31::new(1584694829),
        ])
    }

    fn two_adic_group_order() -> usize {
        1 << 31
    }
}

// TODO: for now it is a dirty hack and should definitely be derived automatically
impl FieldExtension<M31> for Mersenne31Quartic {
    const DEGREE: usize = 4;

    fn mul_assign_by_base(&mut self, elem: &M31) -> &mut Self {
        self.coeffs.iter_mut().for_each(|coef| {
            coef.mul_assign_by_base(elem);
        });
        self
    }

    fn into_coeffs_in_base(self) -> [M31; 4] {
        let Mersenne31Quartic { coeffs } = self;
        let [c0, c1] = coeffs[0].into_coeffs_in_base();
        let [c2, c3] = coeffs[1].into_coeffs_in_base();
        [c0, c1, c2, c3]
    }

    fn from_coeffs_in_base(coefs: &[M31]) -> Self {
        let c0 = Mersenne31Complex::from_coeffs_in_base(&coefs[0..2]);
        let c1 = Mersenne31Complex::from_coeffs_in_base(&coefs[2..4]);
        Self { coeffs: [c0, c1] }
    }

    #[inline(always)]
    fn from_coeffs_in_base_ref(coeffs: &[&M31]) -> Self {
        let c0 = Mersenne31Complex::from_coeffs_in_base_ref(&coeffs[0..2]);
        let c1 = Mersenne31Complex::from_coeffs_in_base_ref(&coeffs[2..4]);
        Self { coeffs: [c0, c1] }
    }

    fn coeffs_in_base(&self) -> &[M31] {
        unsafe { std::mem::transmute::<&[Mersenne31Complex], &[M31]>(&self.coeffs) }
    }

    fn add_assign_base(&mut self, elem: &M31) -> &mut Self {
        self.coeffs[0].add_assign_base(elem);
        self
    }

    fn sub_assign_base(&mut self, elem: &M31) -> &mut Self {
        self.coeffs[0].sub_assign_base(elem);
        self
    }

    fn from_base(elem: M31) -> Self {
        let c0 = Mersenne31Complex::from_base(elem);
        Self {
            coeffs: [c0, Mersenne31Complex::ZERO],
        }
    }

    fn get_coef_mut(&mut self, idx: usize) -> &mut M31 {
        let major_idx = if idx < 2 { 0 } else { 1 };
        self.coeffs[major_idx].get_coef_mut(idx % 2)
    }
}

pub fn rand_fp_from_rng<R: rand::Rng>(rng: &mut R) -> M31 {
    M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER))
}

pub fn rand_fp2_from_rng<R: rand::Rng>(rng: &mut R) -> Mersenne31Complex {
    let a = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    let b = M31::from_u64_unchecked(rng.gen_range(0..M31::ORDER));
    Mersenne31Complex::new(a, b)
}

/*
pub struct Mersenne31Engine {}

impl Engine for Mersenne31Engine {
    type BaseField = M31;
    type LdeOutputField = Mersenne31Complex;
    type ChallengeField = Mersenne31Quartic;

    // NB: FROBENIUS MULTIPLIERS were generated with the help of the following sage script:
    // p = 2^31 - 1
    // base_field = GF(p)
    // R.<x> = base_field[]
    // lde_field = base_field.extension(x^2+1,'i')

    // elem = lde_field(2 + x)
    // degree = (p-1)/2
    // alpha = elem ^ degree
    // print(alpha^(p+1))
    // print(alpha^(p^2 + p + 1))

    const FROBENIUS_MULTIPLIER_ONE: Self::LdeOutputField = Mersenne31Complex {
        coeffs: [
            M31::new(21189756),
            M31::new(42379512),
        ],
    };
    const FROBENIUS_MULTIPLIER_TWO: Self::LdeOutputField = Mersenne31Complex {
        coeffs: [M31::new(2147483646), M31::ZERO],
    };
    const FROBENIUS_MULTIPLIER_THREE: Self::LdeOutputField = Mersenne31Complex {
        coeffs: [
            M31::new(2126293891),
            M31::new(2105104135),
        ],
    };
}
*/
*/

impl PartialEq for Mersenne31Field {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl Eq for M31 {}

impl Hash for M31 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u32(self.as_reduced_u32())
    }
}

impl Ord for M31 {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for M31 {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for M31 {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}

impl Debug for M31 {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

/*
impl std::fmt::Debug for Mersenne31Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "F2[{}, {}]", self.coeffs[0], self.coeffs[1])
    }
}

impl std::fmt::Display for Mersenne31Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "F2[{}, {}]", self.coeffs[0], self.coeffs[1])
    }
}

impl std::fmt::Debug for Mersenne31Quartic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "F4[{}, {}, {}, {}]",
            self.coeffs[0].coeffs[0],
            self.coeffs[0].coeffs[1],
            self.coeffs[1].coeffs[0],
            self.coeffs[1].coeffs[1]
        )
    }
}

impl std::fmt::Display for Mersenne31Quartic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "F4[{}, {}, {}, {}]",
            self.coeffs[0].coeffs[0],
            self.coeffs[0].coeffs[1],
            self.coeffs[1].coeffs[0],
            self.coeffs[1].coeffs[1]
        )
    }
}
*/
