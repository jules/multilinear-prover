//! Multilinear FRI in the extension.

use super::*;
use crate::{
    field::{ChallengeField, Field},
    linear_code::LinearCode,
    merkle_tree::MerkleTree,
    polynomial::mle::MultilinearExtension,
    transcript::Transcript,
};
use core::marker::PhantomData;

pub struct FRI<F: Field, E: ChallengeField<F>, T: Transcript<F>, LC: LinearCode<F>> {
    _f_marker: PhantomData<F>,
    _e_marker: PhantomData<E>,
    _t_marker: PhantomData<T>,
    _lc_marker: PhantomData<LC>,
}
