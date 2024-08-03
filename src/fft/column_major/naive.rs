use super::*;
use crate::field::*;

// usually used to get bitreversed evaluation form from the monomial one
pub fn fft_natural_to_bitreversed<F: Field, E: Field + ChallengeField<F>>(
    input: &mut [E],
    coset_initial: F,
    coset_step: F,
    twiddles: &[F],
) {
    debug_assert!(input.len().is_power_of_two());
    debug_assert!(input.len() == twiddles.len() * 2);

    if coset_initial != F::ONE || coset_step != F::ONE {
        distribute_powers_serial(input, coset_initial, coset_step);
    }

    let log_n = input.len().trailing_zeros();
    cache_friendly_ntt_natural_to_bitreversed(input, log_n, twiddles);
}

pub fn fft_natural_to_natural<F: Field, E: Field + ChallengeField<F>>(
    input: &mut [E],
    coset_initial: F,
    coset_step: F,
    twiddles: &[F],
) {
    debug_assert!(input.len().is_power_of_two());
    debug_assert!(input.len() == twiddles.len() * 2);

    if coset_initial != F::ONE || coset_step != F::ONE {
        distribute_powers_serial(input, coset_initial, coset_step);
    }

    let log_n = input.len().trailing_zeros();
    cache_friendly_ntt_natural_to_bitreversed(input, log_n, twiddles);
    bitreverse_enumeration_inplace(input);
}

// usually used to get monomial from natural bitreversed evaluation
pub fn ifft_natural_to_natural<
    P: PrimeField,
    F: Field,
    E: Field + ChallengeField<P> + ChallengeField<F>,
>(
    input: &mut [E],
    coset: F,
    twiddles: &[F],
) {
    debug_assert!(input.len().is_power_of_two());
    debug_assert!(input.len() == twiddles.len() * 2);

    let log_n = input.len().trailing_zeros();
    serial_ct_ntt_natural_to_bitreversed(input, log_n, twiddles);
    bitreverse_enumeration_inplace(input);

    if coset != F::ONE {
        let coset = coset.inverse().expect("inverse of coset must exist");
        distribute_powers_serial(input, F::ONE, coset);
    }

    if input.len() > 1 {
        let n_inv = P::from_u64_with_reduction(input.len() as u64)
            .inverse()
            .unwrap();
        let mut i = 0;
        let work_size = input.len();
        while i < work_size {
            input[i].mul_base(&n_inv);
            i += 1;
        }
    }
}

pub fn serial_ct_ntt_natural_to_bitreversed<F: Field, E: Field + ChallengeField<F>>(
    a: &mut [E],
    log_n: u32,
    omegas_bit_reversed: &[F],
) {
    let n = a.len();
    if n == 1 {
        return;
    }

    if a.len() > 16 {
        debug_assert!(n == omegas_bit_reversed.len() * 2);
    }
    debug_assert!(n == (1 << log_n) as usize);

    let mut pairs_per_group = n / 2;
    let mut num_groups = 1;
    let mut distance = n / 2;

    {
        // special case for omega = 1
        debug_assert!(num_groups == 1);
        let idx_1 = 0;
        let idx_2 = pairs_per_group;

        let mut j = idx_1;

        while j < idx_2 {
            let u = a[j];
            let v = a[j + distance];

            let mut tmp = u;
            tmp.sub_assign(&v);

            a[j + distance] = tmp;
            a[j].add_assign(&v);

            j += 1;
        }

        pairs_per_group /= 2;
        num_groups *= 2;
        distance /= 2;
    }

    while num_groups < n {
        debug_assert!(num_groups > 1);
        let mut k = 0;
        while k < num_groups {
            let idx_1 = k * pairs_per_group * 2;
            let idx_2 = idx_1 + pairs_per_group;
            let s = omegas_bit_reversed[k];

            let mut j = idx_1;
            while j < idx_2 {
                let u = a[j];
                let mut v = a[j + distance];
                v.mul_base(&s);

                let mut tmp = u;
                tmp.sub_assign(&v);

                a[j + distance] = tmp;
                a[j].add_assign(&v);

                j += 1;
            }

            k += 1;
        }

        pairs_per_group /= 2;
        num_groups *= 2;
        distance /= 2;
    }
}

pub fn cache_friendly_ntt_natural_to_bitreversed<F: Field, E: Field + ChallengeField<F>>(
    a: &mut [E],
    log_n: u32,
    omegas_bit_reversed: &[F],
) {
    let n = a.len();
    if n == 1 {
        return;
    }

    // const CACHE_SIZE_LOG: u32 = 20; // 1 MB
    const CACHE_SIZE_LOG: u32 = 20; // 1 MB

    const WORD_SIZE_LOG: u32 = 3; // 8 B
    const CACHE_BUNCH_LOG: u32 = CACHE_SIZE_LOG - WORD_SIZE_LOG; // 2^17 B
    let cache_friendly_round = log_n.saturating_sub(CACHE_BUNCH_LOG); // 7 round

    if a.len() > 16 {
        debug_assert!(n == omegas_bit_reversed.len() * 2);
    }
    debug_assert!(n == (1 << log_n) as usize);

    let mut pairs_per_group = n / 2;
    let mut num_groups = 1;
    let mut distance = n / 2;
    let mut round = 0;

    {
        // special case for omega = 1
        debug_assert!(num_groups == 1);
        let idx_1 = 0;
        let idx_2 = pairs_per_group;

        let mut j = idx_1;
        while j < idx_2 {
            let mut u = a[j];
            let v = a[j + distance];
            u.sub_assign(&v);
            a[j + distance] = u;
            a[j].add_assign(&v);
            j += 1;
        }
        pairs_per_group /= 2;
        num_groups *= 2;
        distance /= 2;
        round += 1;
    }

    while round < cache_friendly_round {
        debug_assert!(num_groups > 1);
        let mut k = 0;
        while k < num_groups {
            let idx_1 = k * pairs_per_group * 2;
            let idx_2 = idx_1 + pairs_per_group;
            let s = omegas_bit_reversed[k];

            let mut j = idx_1;
            while j < idx_2 {
                let mut u = a[j];
                let mut v = a[j + distance];
                v.mul_base(&s);
                u.sub_assign(&v);
                a[j + distance] = u;
                a[j].add_assign(&v);
                j += 1;
            }
            k += 1;
        }

        pairs_per_group /= 2;
        num_groups *= 2;
        distance /= 2;
        round += 1;
    }
    let mut cache_bunch = 0;
    while cache_bunch < num_groups {
        // num_groups=128 // round loop
        let mut pairs_per_group_in_cache = pairs_per_group;
        let mut distance_in_cache = distance;
        let mut num_groups_in_cache = 1;
        let num_rounds_in_cache = log_n - round; // 17

        let mut round = 0;
        while round < num_rounds_in_cache {
            // experiment

            let mut k = 0;
            while k < num_groups_in_cache {
                // group loop
                let idx_1 = cache_bunch * pairs_per_group * 2 + k * pairs_per_group_in_cache * 2;
                let idx_2 = idx_1 + pairs_per_group_in_cache;
                let s = omegas_bit_reversed[cache_bunch * num_groups_in_cache + k];

                let mut j = idx_1;
                while j < idx_2 {
                    let mut u = a[j];
                    let mut v = a[j + distance_in_cache];
                    v.mul_base(&s);
                    u.sub_assign(&v);
                    a[j + distance_in_cache] = u;
                    a[j].add_assign(&v);

                    j += 1;
                }
                k += 1;
            }
            pairs_per_group_in_cache /= 2;
            num_groups_in_cache *= 2;
            distance_in_cache /= 2;
            round += 1;
        }
        cache_bunch += 1;
    }
}
