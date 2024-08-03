use super::*;
use crate::fft::*;
use crate::field::m31::{complex::M31_2, quartic::M31_4, M31};
use crate::trace_holder::RowMajorTrace;
use crate::worker::Worker;
use core::marker::PhantomData;
use std::alloc::Global;

pub struct ReedSolomonCode<F: Field, E: ChallengeField<F>> {
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
}

impl LinearCode<M31, M31_4> for ReedSolomonCode<M31, M31_4> {
    const BLOWUP: usize = 2;

    fn encode(els: &[M31]) -> Vec<M31> {
        let worker = Worker::new_with_num_threads(8);
        let side = els.len().isqrt();
        println!("{}", els.len());
        println!("{side}");
        let mut trace = RowMajorTrace::<32, _>::new_zeroed_for_size(side, side, Global);
        let twiddles = Twiddles::<M31_2, Global>::new(side, &worker);
        // copy to trace
        let mut row_view = trace.row_view(0..side);
        for i in 0..side {
            let row = row_view.current_row();
            row[0] = els[i];
            row_view.advance_row();
        }

        // XXX needs to be pre-computed
        let lde_precomputations =
            LdePrecomputations::<Global>::new(side, Self::BLOWUP, &[0, 1], &worker);
        let scales = &lde_precomputations.domain_bound_precomputations[0]
            .as_ref()
            .unwrap()
            .bitreversed_powers[0];
        parallel_row_major_full_line_fft_dit(
            &mut trace,
            &twiddles.forward_twiddles_not_bitreversed,
            &scales,
            &worker,
        );

        trace.as_slice().to_vec()
    }

    fn encode_ext(els: &[M31_4]) -> Vec<M31_4> {
        let mut els = els.to_vec();
        els.extend(els.clone());
        els.extend(els.clone());
        els
    }
}
