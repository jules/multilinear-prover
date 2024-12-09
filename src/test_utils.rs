use crate::{field::Field, polynomial::mle::MultilinearExtension};
use rand::Rng;

// #[derive(Default)]
// pub struct MockTranscript<F: Field> {
//     counter: usize,
//     _marker: PhantomData<F>,
// }
//
// impl<F: Field> Transcript<F> for MockTranscript<F> {
//     fn draw_challenge(&mut self) -> F {
//         self.counter += 1;
//         F::from_usize(self.counter)
//     }
//     fn draw_challenge_ext<E: ChallengeField<F>>(&mut self) -> E {
//         self.counter += 1;
//         E::from(F::from_usize(self.counter))
//     }
//     fn draw_bits(&mut self, _bits: usize) -> usize {
//         0
//     }
//     fn observe_witness(&mut self, _witness: F) {}
//     fn observe_witnesses(&mut self, _witness: &[F]) {}
// }

pub fn rand_poly<F: Field>(order_bits: u32) -> MultilinearExtension<F> {
    let mut evals = vec![F::default(); 2u32.pow(order_bits) as usize];
    evals
        .iter_mut()
        .for_each(|e| *e = F::from_usize(rand::thread_rng().gen::<usize>()));
    MultilinearExtension::new(evals)
}
