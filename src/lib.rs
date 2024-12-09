#![feature(generic_const_exprs)]
#![feature(associated_type_defaults)]
#![feature(isqrt)]
#![feature(allocator_api)]
#![feature(const_mut_refs)]
#![feature(slice_swap_unchecked)]
#![feature(const_swap)]
#![feature(raw_slice_split)]
#![cfg_attr(target_arch = "aarch64", feature(stdarch_aarch64_prefetch))]
#![allow(refining_impl_trait)]

mod arguments;
mod fft;
mod field;
mod iop;
mod linear_code;
mod merkle_tree;
mod pcs;
mod polynomial;
mod transcript;
mod univariate_utils;

#[cfg(test)]
mod test_utils;
