use super::*;

#[derive(Debug)]
pub struct RowMajorTrace<const N: usize, A: Allocator + Clone> {
    pub ptr: *mut M31,
    width: usize,
    length: usize,
    pub padded_width: usize,
    allocator: A,
}

unsafe impl<const N: usize, A: Allocator + Clone> Send for RowMajorTrace<N, A> {}
unsafe impl<const N: usize, A: Allocator + Clone> Sync for RowMajorTrace<N, A> {}

impl<const N: usize, A: Allocator + Clone> Drop for RowMajorTrace<N, A> {
    fn drop(&mut self) {
        unsafe {
            let required_size = self.length * self.padded_width;
            let reconstructed_capacity = Vec::<u8, A>::from_raw_parts_in(
                self.ptr.cast(),
                required_size,
                required_size,
                self.allocator.clone(),
            );
            drop(reconstructed_capacity);
        }
    }
}

impl<const N: usize, A: Allocator + Clone> Clone for RowMajorTrace<N, A> {
    fn clone(&self) -> Self {
        let new = Self::new_uninit_for_size(self.length, self.width, self.allocator());
        // memcopy
        let num_elements = self.padded_width * self.length;
        unsafe { core::ptr::copy_nonoverlapping(self.ptr.cast_const(), new.ptr, num_elements) }

        new
    }

    fn clone_from(&mut self, source: &Self) {
        assert_eq!(self.length, source.length);
        assert_eq!(self.width, source.width);
        let num_elements = self.padded_width * self.length;
        unsafe { core::ptr::copy_nonoverlapping(source.ptr.cast_const(), self.ptr, num_elements) }
    }
}

impl<const N: usize, A: Allocator + Clone> RowMajorTrace<N, A> {
    pub fn as_slice(&self) -> &[M31] {
        unsafe {
            core::slice::from_raw_parts(self.ptr.cast_const(), self.length * self.padded_width)
        }
    }
    pub fn len(&self) -> usize {
        self.length
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn allocator(&self) -> A {
        self.allocator.clone()
    }

    pub fn num_column_chunks(&self, width: usize) -> usize {
        assert!(width.is_power_of_two());
        assert!(self.padded_width % width == 0);

        self.padded_width / width
    }

    pub fn new_uninit_for_size(rows: usize, columns: usize, allocator: A) -> Self {
        assert!(N.is_power_of_two());
        let padded_columns = columns.next_multiple_of(N);
        let required_size = padded_columns * rows * std::mem::size_of::<M31>();
        let num_elements = required_size.next_multiple_of(PAGE_SIZE) / PAGE_SIZE;
        let capacity = Vec::<Aligner, A>::with_capacity_in(num_elements, allocator);
        let (ptr, _, _, alloc) = capacity.into_raw_parts_with_alloc();

        Self {
            ptr: ptr.cast(),
            width: columns,
            length: rows,
            padded_width: padded_columns,
            allocator: alloc,
        }
    }

    pub fn new_zeroed_for_size(rows: usize, columns: usize, allocator: A) -> Self {
        let new = Self::new_uninit_for_size(rows, columns, allocator);
        unsafe {
            let start = new.ptr.cast::<u8>();
            core::ptr::write_bytes(start, 0u8, new.length * new.padded_width);
        }

        new
    }

    #[track_caller]
    pub fn row_view(&self, range: std::ops::Range<usize>) -> RowMajorTraceView<N> {
        assert!(range.end <= self.length);
        let length = range.len();

        let ptr = unsafe { self.ptr.add(self.padded_width * range.start) };

        RowMajorTraceView {
            ptr,
            padded_width: self.padded_width,
            width: self.width,
            length,
            full_trace_starting_ptr: self.ptr,
            full_trace_ending_ptr: unsafe { self.ptr.add(self.padded_width * self.length) },
        }
    }

    #[track_caller]
    pub fn column_view(&self, offset: usize, width: usize) -> RowMajorTraceColumnsView<N> {
        assert!(width > 0);
        assert!(
            width + offset <= self.padded_width,
            "trace has padded width {}, but caller requested {} columns at offset {}",
            self.padded_width,
            width,
            offset
        );

        let ptr = unsafe { self.ptr.add(offset) };

        RowMajorTraceColumnsView {
            ptr,
            viewed_width: width,
            offset_per_row: self.padded_width,
            length: self.length,
        }
    }

    pub fn tile_view(
        &self,
        rows: std::ops::Range<usize>,
        width: usize,
        offset: usize,
    ) -> RowMajorTraceColumnsView<N> {
        assert!(rows.end <= self.length);
        let length = rows.len();
        assert!(width * (offset + 1) <= self.padded_width);

        let ptr = unsafe { self.ptr.add(self.padded_width * rows.start) };
        let ptr = unsafe { ptr.add(width * offset) };

        RowMajorTraceColumnsView {
            ptr,
            viewed_width: width,
            offset_per_row: self.padded_width,
            length: length,
        }
    }

    #[track_caller]
    pub fn column_view_fixed<const M: usize>(
        &self,
        offset: usize,
    ) -> RowMajorTraceFixedColumnsView<N, M> {
        assert!(M.is_power_of_two());
        assert!(
            M * (offset + 1) <= self.padded_width,
            "trying to get offset {} for column view of width {}, but our padded width is {}",
            offset,
            M,
            self.padded_width
        );

        let ptr = unsafe { self.ptr.add(M * offset) };

        RowMajorTraceFixedColumnsView {
            ptr,
            offset_per_row: self.padded_width,
            length: self.length,
        }
    }
}

pub struct RowMajorTraceView<const N: usize> {
    ptr: *mut M31,
    width: usize,
    length: usize,
    padded_width: usize,
    full_trace_starting_ptr: *const M31,
    full_trace_ending_ptr: *const M31,
}

unsafe impl<const N: usize> Send for RowMajorTraceView<N> {}
unsafe impl<const N: usize> Sync for RowMajorTraceView<N> {}

impl<const N: usize> RowMajorTraceView<N> {
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.length
    }

    #[inline(always)]
    pub fn padded_width(&self) -> usize {
        self.padded_width
    }

    #[inline(always)]
    pub fn current_row_ref(&self) -> &[M31] {
        unsafe { core::slice::from_raw_parts(self.ptr.cast_const(), self.width) }
    }

    #[inline(always)]
    pub fn current_and_next_row_ref(&self) -> (&[M31], &[M31]) {
        unsafe {
            let this_row_ptr = self.ptr.cast_const();
            let mut next_row_ptr = this_row_ptr.add(self.padded_width);
            if next_row_ptr >= self.full_trace_ending_ptr {
                next_row_ptr = self.full_trace_starting_ptr;
            }
            (
                core::slice::from_raw_parts(this_row_ptr, self.width),
                core::slice::from_raw_parts(next_row_ptr, self.width),
            )
        }
    }

    #[inline(always)]
    pub fn current_and_previous_row_ref(&self) -> (&[M31], &[M31]) {
        unsafe {
            let this_row_ptr = self.ptr.cast_const();
            let previous_row_ptr = if this_row_ptr == self.full_trace_starting_ptr {
                self.full_trace_ending_ptr.sub(self.padded_width)
            } else {
                this_row_ptr.sub(self.padded_width)
            };
            (
                core::slice::from_raw_parts(this_row_ptr, self.width),
                core::slice::from_raw_parts(previous_row_ptr, self.width),
            )
        }
    }

    #[inline(always)]
    pub fn current_row(&mut self) -> &mut [M31] {
        unsafe { core::slice::from_raw_parts_mut(self.ptr, self.width) }
    }

    #[inline(always)]
    pub fn current_row_padded(&mut self) -> &mut [M31] {
        unsafe { core::slice::from_raw_parts_mut(self.ptr, self.padded_width) }
    }

    #[inline(always)]
    pub fn advance_row(&mut self) {
        let offset = self.padded_width;
        self.length -= 1;
        unsafe {
            self.ptr = self.ptr.add(offset);
        }
    }

    pub fn as_slice(&self) -> &[M31] {
        unsafe {
            core::slice::from_raw_parts(self.ptr.cast_const(), self.length * self.padded_width)
        }
    }
}

pub struct RowMajorTraceColumnsView<const N: usize> {
    ptr: *mut M31,
    offset_per_row: usize,
    viewed_width: usize,
    length: usize,
}

unsafe impl<const N: usize> Send for RowMajorTraceColumnsView<N> {}
unsafe impl<const N: usize> Sync for RowMajorTraceColumnsView<N> {}

impl<const N: usize> RowMajorTraceColumnsView<N> {
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.length
    }

    #[inline(always)]
    pub fn width(&self) -> usize {
        self.viewed_width
    }

    #[inline(always)]
    pub fn current_row(&self) -> &[M31] {
        unsafe { core::slice::from_raw_parts(self.ptr.cast_const(), self.viewed_width) }
    }

    #[inline(always)]
    pub fn get_row(&self, idx: usize) -> &[M31] {
        debug_assert!(idx < self.length);
        unsafe {
            let ptr = self.ptr.add(self.offset_per_row * idx);
            core::slice::from_raw_parts(ptr.cast_const(), self.viewed_width)
        }
    }

    #[inline(always)]
    pub fn get_row_mut(&mut self, idx: usize) -> &mut [M31] {
        debug_assert!(idx < self.length);
        unsafe {
            let ptr = self.ptr.add(self.offset_per_row * idx);
            core::slice::from_raw_parts_mut(ptr, self.viewed_width)
        }
    }

    #[inline(always)]
    pub fn advance_row(&mut self) {
        unsafe {
            self.ptr = self.ptr.add(self.offset_per_row);
        }
    }

    #[inline(always)]
    pub fn advance_many(&mut self, offset: usize) {
        unsafe {
            self.ptr = self.ptr.add(self.offset_per_row * offset);
        }
    }
}

pub struct RowMajorTraceFixedColumnsView<const N: usize, const M: usize> {
    ptr: *mut M31,
    offset_per_row: usize,
    length: usize,
}

unsafe impl<const N: usize, const M: usize> Send for RowMajorTraceFixedColumnsView<N, M> {}
unsafe impl<const N: usize, const M: usize> Sync for RowMajorTraceFixedColumnsView<N, M> {}

impl<const N: usize, const M: usize> RowMajorTraceFixedColumnsView<N, M> {
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.length
    }

    #[inline(always)]
    pub fn current_row(&mut self) -> &mut [M31; M] {
        unsafe { &mut *self.ptr.cast::<[M31; M]>() }
    }

    #[inline(always)]
    pub fn get_row(&self, idx: usize) -> &[M31; M] {
        debug_assert!(idx < self.length);
        unsafe {
            let ptr = self.ptr.add(self.offset_per_row * idx);
            &*ptr.cast_const().cast::<[M31; M]>()
        }
    }

    #[inline(always)]
    pub fn get_row_mut(&mut self, idx: usize) -> &mut [M31; M] {
        debug_assert!(idx < self.length);
        unsafe {
            let ptr = self.ptr.add(self.offset_per_row * idx);
            &mut *ptr.cast::<[M31; M]>()
        }
    }

    #[inline(always)]
    pub fn advance_row(&mut self) {
        unsafe {
            self.ptr = self.ptr.add(self.offset_per_row);
        }
    }
}
