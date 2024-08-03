use crate::field::*;

// for: https://www.robinscheibler.org/2013/02/13/real-fft.html
pub fn pack_reals_as_complex<F: Field, E: ChallengeField<F>>(
    left: &[F],
    right: &[F],
    res_buf: &mut [E],
) {
    assert_eq!(left.len(), right.len());
    assert_eq!(left.len(), res_buf.len());
    assert_eq!(E::DEGREE, 2);

    left.iter()
        .zip(right.iter())
        .zip(res_buf.iter_mut())
        .for_each(|((left, right), out)| *out = E::new(vec![*left, *right]))
}

pub fn unpack_complex_into_reals<F: Field, E: ChallengeField<F>>(
    input: &[E],
    left_res: &mut [F],
    right_res: &mut [F],
) {
    assert_eq!(left_res.len(), right_res.len());
    assert_eq!(left_res.len(), input.len());
    assert_eq!(E::DEGREE, 2);

    left_res
        .iter_mut()
        .zip(right_res.iter_mut())
        .zip(input.iter())
        .for_each(|((left, right), input)| {
            let elems_in_base = Into::<Vec<F>>::into(*input);
            *left = elems_in_base[0];
            *right = elems_in_base[1];
        })
}

// swap via lookup table that itself fits into cache line
const SMALL_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE: usize = 6;
const SMALL_BITREVERSE_LOOKUP_TABLE: [u8; 1 << SMALL_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE] = const {
    let mut result = [0u8; 1 << SMALL_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE];
    let mut i = 0u64;
    let shift_right = 64 - SMALL_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE;
    while i < (1 << SMALL_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE) {
        let reversed = i.reverse_bits() >> shift_right;
        debug_assert!(reversed <= u8::MAX as u64);
        result[i as usize] = reversed as u8;
        i += 1;
    }

    result
};

// in this case we can easily swap bytes, and then swap bits in bytes
const MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE: usize = 8;
const MEDIUM_BITREVERSE_LOOKUP_TABLE: [u8; 1 << MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE] = const {
    let mut result = [0u8; 1 << MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE];
    let mut i = 0u64;
    let shift_right = 64 - MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE;
    while i < (1 << MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE) {
        let reversed = i.reverse_bits() >> shift_right;
        debug_assert!(reversed <= u8::MAX as u64);
        result[i as usize] = reversed as u8;
        i += 1;
    }

    result
};

// This operation is so cache-unfriendly, that parallelism is not used here
pub const fn bitreverse_enumeration_inplace<T>(input: &mut [T]) {
    if input.len() == 0 {
        return;
    }
    assert!(input.len().is_power_of_two());

    if input.len() <= SMALL_BITREVERSE_LOOKUP_TABLE.len() {
        bitreverse_enumeration_inplace_via_small_lookup(input);
    } else if input.len() <= MEDIUM_BITREVERSE_LOOKUP_TABLE.len() {
        bitreverse_enumeration_inplace_via_medium_lookup(input);
    } else {
        bitreverse_enumeration_inplace_hybrid(input);
    }
}

const fn bitreverse_enumeration_inplace_via_small_lookup<T>(input: &mut [T]) {
    assert!(input.len().is_power_of_two());
    assert!(input.len() <= SMALL_BITREVERSE_LOOKUP_TABLE.len());

    let shift_to_cleanup =
        (SMALL_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE as u32) - input.len().trailing_zeros();

    let mut i = 0;
    let work_size = input.len();
    while i < work_size {
        let mut j = SMALL_BITREVERSE_LOOKUP_TABLE[i] as usize;
        j >>= shift_to_cleanup; // if our table size is larger than work size
        if i < j {
            unsafe { input.swap_unchecked(i, j) };
        }

        i += 1;
    }
}

const fn bitreverse_enumeration_inplace_via_medium_lookup<T>(input: &mut [T]) {
    assert!(input.len().is_power_of_two());
    assert!(input.len() <= MEDIUM_BITREVERSE_LOOKUP_TABLE.len());

    let shift_to_cleanup =
        (MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE as u32) - input.len().trailing_zeros();

    let mut i = 0;
    let work_size = input.len();
    while i < work_size {
        let mut j = MEDIUM_BITREVERSE_LOOKUP_TABLE[i] as usize;
        j >>= shift_to_cleanup; // if our table size is larger than work size
        if i < j {
            unsafe { input.swap_unchecked(i, j) };
        }

        i += 1;
    }
}

const fn bitreverse_enumeration_inplace_hybrid<T>(input: &mut [T]) {
    assert!(input.len().is_power_of_two());
    assert!(input.len() > MEDIUM_BITREVERSE_LOOKUP_TABLE.len());
    assert!(input.len() <= 1usize << 31); // a reasonable upper bound to use u32 internally

    // there is a function usize::reverse_bits(), but if one looks into the compiler then
    // will see that it's something like (sorry for C code)
    // ```
    //     uint32_t bit_reverse32(uint32_t x)
    // {
    //     x = (x >> 16) | (x << 16);
    //     x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    //     x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    //     x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    //     return ((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1);
    // }
    // ```

    // since we bitreverse a continuous set of indexes, we can save a little by
    // doing two loops, such that one bitreverses (naively) some common bits,
    // and one that bitreversed uncommon via lookup

    let log_n = input.len().trailing_zeros();
    let common_part_log_n = log_n - (MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE as u32);

    // double loop. Note the swapping approach:
    // - lowest bits become highest bits and change every time
    // - highest bits change become lowest bits and change rarely
    // so our "i" counter is a counter over highest bits, and our source is in the form (i << 8) + j
    // and our dst is (reversed_j << common_part_log_n) + reversed_i
    // and since our source and destination are symmetrical we can formally swap them
    // and have our writes cache-friendly
    let mut i = 0;
    let work_size = 1u32 << common_part_log_n;
    while i < work_size {
        // bitreversing byte by byte
        let mut bytes = i.swap_bytes().to_le_bytes();
        bytes[0] = 0;
        bytes[1] = MEDIUM_BITREVERSE_LOOKUP_TABLE[bytes[1] as usize];
        bytes[2] = MEDIUM_BITREVERSE_LOOKUP_TABLE[bytes[2] as usize];
        bytes[3] = MEDIUM_BITREVERSE_LOOKUP_TABLE[bytes[3] as usize];
        let reversed_i = u32::from_le_bytes(bytes) >> (32 - common_part_log_n);

        debug_assert!(reversed_i == i.reverse_bits() >> (32 - common_part_log_n));

        let mut j = 0;
        while j < MEDIUM_BITREVERSE_LOOKUP_TABLE.len() {
            let reversed_j = MEDIUM_BITREVERSE_LOOKUP_TABLE[j];
            let dst = ((i as usize) << MEDIUM_BITREVERSE_LOOKUP_TABLE_LOG_2_SIZE) | j;
            let src = ((reversed_j as usize) << common_part_log_n) | (reversed_i as usize);
            if dst < src {
                unsafe { input.swap_unchecked(src, dst) };
            }

            j += 1;
        }

        i += 1;
    }
}
