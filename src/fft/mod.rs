use crate::field::m31::M31;

pub mod column_major;
pub mod field_utils;
pub mod row_major;
pub mod utils;

pub use self::column_major::*;
pub use self::field_utils::*;
pub use self::row_major::*;
pub use self::utils::*;

pub trait GoodAllocator:
    std::alloc::Allocator + Clone + Default + Send + Sync + std::fmt::Debug
{
}
impl GoodAllocator for std::alloc::Global {}

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
pub const CACHE_LINE_WIDTH: usize = 128;

#[cfg(not(all(target_os = "macos", target_arch = "aarch64")))]
pub const CACHE_LINE_WIDTH: usize = 64;

pub const L1_CACHE_SIZE: usize = 1 << 17;

pub const CACHE_LINE_MULTIPLE: usize = const {
    assert!(core::mem::size_of::<M31>() >= core::mem::align_of::<M31>());

    CACHE_LINE_WIDTH / core::mem::size_of::<M31>()
};

#[cfg(target_arch = "aarch64")]
#[inline(always)]
unsafe fn prefetch_next_line(ptr: *const M31) {
    use core::arch::aarch64::{_PREFETCH_LOCALITY1, _PREFETCH_WRITE};
    core::arch::aarch64::_prefetch::<_PREFETCH_WRITE, _PREFETCH_LOCALITY1>(ptr.cast());
}

#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn prefetch_next_line(ptr: *const M31) {
    todo!()
}
