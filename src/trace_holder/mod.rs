use crate::field::m31::M31;
use std::alloc::Allocator;

#[derive(Clone, Copy, Debug)]
#[repr(align(16384))] // up to M1 family paging
struct Aligner {
    _inner: [u8; PAGE_SIZE],
}

pub const PAGE_SIZE: usize = 16384;

pub mod row_major;

pub use self::row_major::*;
