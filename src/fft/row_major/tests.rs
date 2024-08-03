use std::alloc::Global;

use super::*;

pub const TEST_WIDTH: usize = FFT_UNROLL_FACTOR;
// pub const TEST_WIDTH: usize = CACHE_LINE_MULTIPLE;

//#[test]
//fn test_partial_ifft() {
//    let log_n = 5;
//    let num_cores = 1;
//
//    let worker = Worker::new_with_num_threads(num_cores);
//    let trace_len = 1 << log_n;
//    let num_columns = 4;
//
//    let mut trace =
//        RowMajorTrace::<TEST_WIDTH, _>::new_zeroed_for_size(trace_len, num_columns, Global);
//
//    let values: Vec<_> = (0..trace_len)
//        .map(|i| M31(i as u32))
//        .map(|el| M31_2::from_base(el))
//        .collect();
//
//    let twiddles = Twiddles::<M31_2, Global>::new(trace_len, &worker);
//
//    let mut reference_result = values.clone();
//    serial_ct_ntt_natural_to_bitreversed(&mut reference_result, log_n, &twiddles.inverse_twiddles);
//
//    // copy to the trace
//    // let scale = M31(trace_len as u32).inverse().unwrap();
//
//    let mut row_view = trace.row_view(0..trace_len);
//    for i in 0..trace_len {
//        let row = row_view.current_row();
//        row[0] = values[i].c0;
//        row_view.advance_row()
//    }
//    // try to get some values
//
//    parallel_row_major_full_line_partial_ifft(&mut trace, &twiddles.inverse_twiddles, &worker);
//
//    let mut row_view = trace.row_view(0..trace_len);
//    for row_idx in 0..trace_len {
//        let reference_value_c0 = reference_result[row_idx].c0;
//        let reference_value_c1 = reference_result[row_idx].c1;
//        let row = row_view.current_row_ref();
//        assert_eq!(reference_value_c0, row[0], "c0 failed at row {}", row_idx);
//        assert_eq!(reference_value_c1, row[1], "c1 failed at row {}", row_idx);
//        row_view.advance_row();
//    }
//}

#[test]
fn test_fft_roundtrip() {
    let log_n = 20;
    let num_cores = 8;

    let worker = Worker::new_with_num_threads(num_cores);
    let trace_len = 1 << log_n;
    let num_columns = 32;

    let mut trace =
        RowMajorTrace::<TEST_WIDTH, _>::new_zeroed_for_size(trace_len, num_columns, Global);

    let values: Vec<_> = (0..trace_len)
        .map(|i| M31(i as u32))
        .map(|el| M31_2::from_base(el))
        .collect();

    let twiddles = Twiddles::<M31_2, Global>::new(trace_len, &worker);

    // copy to the trace
    // let scale = M31(trace_len as u32).inverse().unwrap();

    let mut row_view = trace.row_view(0..trace_len);
    for i in 0..trace_len {
        let row = row_view.current_row();
        row[0] = values[i].c0;
        row_view.advance_row();
    }
    // try to get some values
    let original_values = trace.clone();

    parallel_row_major_full_line_partial_ifft(&mut trace, &twiddles.inverse_twiddles, &worker);

    // let lde_precomputations = LdePrecomputations::new(trace_len, 2, &[0, 1], &worker);

    let lde_precomputations = LdePrecomputations::<Global>::new(trace_len, 2, &[0, 1], &worker);

    let scales = &lde_precomputations.domain_bound_precomputations[0]
        .as_ref()
        .unwrap()
        .bitreversed_powers[0];

    // let scale = M31(trace_len as u32).inverse().unwrap();
    // let scales = vec![M31_2::from_base(scale); trace_len];

    parallel_row_major_full_line_fft_dit(
        &mut trace,
        &twiddles.forward_twiddles_not_bitreversed,
        &scales,
        &worker,
    );

    let mut reference_row_view = original_values.row_view(0..trace_len);
    let mut row_view = trace.row_view(0..trace_len);
    for row_idx in 0..trace_len {
        let row = row_view.current_row_ref();
        let reference_row = reference_row_view.current_row_ref();
        assert_eq!(reference_row, row, "failed at row index {}", row_idx);
        reference_row_view.advance_row();
        row_view.advance_row();
    }
}

#[test]
fn test_fft_cross_domains_round_trip() {
    let log_n = 20;
    let num_cores = 8;

    let worker = Worker::new_with_num_threads(num_cores);
    let trace_len = 1 << log_n;
    let num_columns = 32;

    let mut trace =
        RowMajorTrace::<TEST_WIDTH, _>::new_zeroed_for_size(trace_len, num_columns, Global);

    let values: Vec<_> = (0..trace_len)
        .map(|i| M31(i as u32))
        .map(|el| M31_2::from_base(el))
        .collect();

    let twiddles = Twiddles::<M31_2, Global>::new(trace_len, &worker);

    let mut row_view = trace.row_view(0..trace_len);
    for i in 0..trace_len {
        let row = row_view.current_row();
        row[0] = values[i].c0;
        row_view.advance_row();
    }

    adjust_to_zero_c0_var_length(&mut trace, 0..num_columns, &worker);

    // try to get some values
    let original_values = trace.clone();

    // partial IFFT
    parallel_row_major_full_line_partial_ifft(&mut trace, &twiddles.inverse_twiddles, &worker);

    let lde_precomputations = LdePrecomputations::<Global>::new(trace_len, 2, &[0, 1], &worker);

    // now FFT into other domain 0 -> 1;
    let source_domain_idx = 0;
    let dest_domain_idx = 1;
    let scales = &lde_precomputations.domain_bound_precomputations[source_domain_idx]
        .as_ref()
        .unwrap()
        .bitreversed_powers[dest_domain_idx];
    parallel_row_major_full_line_fft_dit(
        &mut trace,
        &twiddles.forward_twiddles_not_bitreversed,
        scales,
        &worker,
    );

    // let domain_size_squared = M31_2::from_base(M31((trace_len as u32) * (trace_len as u32)));
    // dbg!(domain_size_squared);
    // for idx in 0..trace_len {
    //     let a = lde_precomputations.domain_bound_precomputations[0].as_ref().unwrap().bitreversed_powers[1][idx];
    //     let b = lde_precomputations.domain_bound_precomputations[1].as_ref().unwrap().bitreversed_powers[0][idx];
    //     let mut t = a;
    //     t.mul_assign(&b);

    //     assert_eq!(t.inverse().unwrap(), domain_size_squared);
    // }

    let mut row_view = trace.row_view(0..trace_len);
    for _ in 0..trace_len {
        let row = row_view.current_row_ref();
        // check that we get what we expect
        assert!(row[1].is_zero());
        for el in row[2..].iter() {
            assert!(el.is_zero());
        }
        row_view.advance_row();
    }

    // because our transformation 0 -> X is special, we must scale by tau^H/2 to get to "normal" values,
    // for all our future purposes

    let tau = &lde_precomputations.domain_bound_precomputations[dest_domain_idx]
        .as_ref()
        .unwrap()
        .coset_offset;
    let scale = tau.pow((trace_len as u32) / 2);
    let mut row_view = trace.row_view(0..trace_len);
    for _ in 0..trace_len {
        let row = row_view.current_row().as_mut_ptr().cast::<M31_2>();
        unsafe {
            let mut el = row.read();
            el.mul_assign(&scale);
            row.write(el);
        }
        row_view.advance_row();
    }

    // IFFT again
    parallel_row_major_full_line_partial_ifft(&mut trace, &twiddles.inverse_twiddles, &worker);

    // now FFT 1 -> 0 back
    let source_domain_idx = 1;
    let dest_domain_idx = 0;
    let scales = &lde_precomputations.domain_bound_precomputations[source_domain_idx]
        .as_ref()
        .unwrap()
        .bitreversed_powers[dest_domain_idx];
    parallel_row_major_full_line_fft_dit(
        &mut trace,
        &twiddles.forward_twiddles_not_bitreversed,
        scales,
        &worker,
    );

    let mut reference_row_view = original_values.row_view(0..trace_len);
    let mut row_view = trace.row_view(0..trace_len);
    for row_idx in 0..trace_len {
        let row = row_view.current_row_ref();
        let reference_row = reference_row_view.current_row_ref();
        assert_eq!(reference_row, row, "failed at row index {}", row_idx);
        reference_row_view.advance_row();
        row_view.advance_row();
    }
}

// #[test]
// fn test_bench_typical_fft() {
//     // let worker = Worker::new_with_num_threads(4);
//     let worker = Worker::new_with_num_threads(8);
//     // let worker = Worker::new_with_num_threads(12);
//     // let worker = Worker::new();
//     let trace_len = 1 << 20;
//     let num_columns = 225;

//     let mut trace = RowMajorTrace::<CACHE_LINE_MULTIPLE, _>::new_zeroed_for_size(
//         trace_len,
//         num_columns,
//         Global,
//     );
//     dbg!(trace.padded_width);
//     dbg!(worker.num_cores);

//     let precomputations = vec![M31_2::ONE; trace_len / 2];

//     let now = std::time::Instant::now();
//     parallel_row_major_full_line_ifft::<_, _, false>(
//         &mut trace,
//         &precomputations,
//         &[],
//         &worker,
//     );
//     dbg!(now.elapsed());
// }

// #[test]
// fn test_bench_typical_quasi_six_step_fft() {
//     let trace_len = 1 << 20;
//     let num_columns = 225;

//     for num_threads in [1, 2, 4, 8, 12, 16] {
//         dbg!(num_threads);
//         let worker = Worker::new_with_num_threads(num_threads);

//         let mut trace = RowMajorTrace::<CACHE_LINE_MULTIPLE, _>::new_zeroed_for_size(
//             trace_len,
//             num_columns,
//             Global,
//         );
//         // dbg!(trace.padded_width);

//         let precomputations = vec![M31_2::ONE; trace_len / 2];

//         let now = std::time::Instant::now();
//         parallel_row_major_full_line_quasi_six_step_ifft(&mut trace, &precomputations, &worker);
//         dbg!(now.elapsed());
//     }
// }

#[test]
fn test_adjust_c0() {
    let trace_len = 1 << 16;
    let num_columns = 8;

    let worker = Worker::new_with_num_threads(1);

    let mut trace = RowMajorTrace::<CACHE_LINE_MULTIPLE, _>::new_zeroed_for_size(
        trace_len,
        num_columns,
        Global,
    );
    // set first row
    let mut t = trace.row_view(0..trace_len);
    let mut i = 0u64;
    for _ in 0..trace_len {
        let row = t.current_row();
        for el in row.iter_mut() {
            *el = M31::from_u64_with_reduction(i);
            i += 1;
        }
        t.advance_row();
    }

    adjust_to_zero_c0_var_length(&mut trace, 0..num_columns, &worker);

    for column in 0..num_columns {
        let mut trace = trace.row_view(0..trace_len);
        let mut sum = M31::ZERO;
        for _ in 0..trace_len {
            let row = trace.current_row_ref();
            sum.add_assign(&row[column]);
            trace.advance_row();
        }
        assert_eq!(sum, M31::ZERO, "invalid for column {}", column);
    }
}
