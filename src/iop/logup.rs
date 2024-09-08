//! Tooling for LogUp. https://eprint.iacr.org/2022/1530.pdf

use crate::{field::Field, polynomial::MultilinearExtension};
use rayon::prelude::*;

/// Compute multiplicities of table values over all columns, assuming the table is injective.
// XXX extremely inefficient we should just take note of multiplicities during witness generation
pub fn compute_multiplicities<F: Field>(
    columns: &[MultilinearExtension<F>],
    table: &MultilinearExtension<F>,
) -> MultilinearExtension<F> {
    let mut multiplicity_evals = Vec::with_capacity(table.len());

    table.evals.iter().for_each(|lookup_value| {
        multiplicity_evals.push(columns.iter().fold(F::ZERO, |acc, column| {
            column.evals.iter().fold(F::ZERO, |mut acc, value| {
                if value == lookup_value {
                    acc.add_assign(&F::ONE);
                }
                acc
            })
        }));
    });

    MultilinearExtension::new(multiplicity_evals)
}

/// Computes the LogUp helper columns. Currently fixed for a sum size of 1.
pub fn compute_helper_columns<F: Field>(
    x_plus_columns: &[MultilinearExtension<F>],
    x_plus_table: &MultilinearExtension<F>,
    multiplicities: &MultilinearExtension<F>,
) -> Vec<MultilinearExtension<F>> {
    let mut helper_columns = Vec::with_capacity(x_plus_columns.len() + 1);

    // First column: m / (x + t).
    let inverses = F::batch_inverse(&x_plus_table.evals).expect("table should be batch invertable");

    helper_columns.push(MultilinearExtension::new(
        multiplicities
            .evals
            .iter()
            .zip(inverses.iter())
            .map(|(m, i)| {
                let mut m = m.clone();
                m.mul_assign(i);
                m
            })
            .collect::<Vec<F>>(),
    ));

    // Rest is just the negative reciprocal of x + column.
    helper_columns.extend(
        x_plus_columns
            .par_iter()
            .map(|column| {
                let mut inverses =
                    F::batch_inverse(&column.evals).expect("column should be batch invertable");
                inverses.iter_mut().for_each(|value| {
                    value.negate();
                });
                MultilinearExtension::new(inverses)
            })
            .collect::<Vec<MultilinearExtension<F>>>(),
    );

    helper_columns
}
