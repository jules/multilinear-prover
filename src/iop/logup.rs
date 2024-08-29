//! Tooling for LogUp. https://eprint.iacr.org/2022/1530.pdf

use crate::{field::Field, polynomial::MultilinearExtension};

/// Compute multiplicities of table values over all columns, assuming the table is injective.
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
    columns: &[MultilinearExtension<F>],
    table: &MultilinearExtension<F>,
    multiplicities: &MultilinearExtension<F>,
    x: F,
) -> Vec<MultilinearExtension<F>> {
    let mut helper_columns = Vec::with_capacity(columns.len());
    // First column: m / (x + t).
    let inverses = table
        .evals
        .iter()
        .map(|value| {
            let mut value = value.clone();
            value.add_assign(&x);
            value.inverse().unwrap()
        })
        .collect::<Vec<F>>();

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

    // Rest is just the reciprocal of x + column.
    columns.iter().for_each(|column| {
        helper_columns.push(MultilinearExtension::new(
            column
                .evals
                .iter()
                .map(|value| {
                    let mut value = value.clone();
                    value.add_assign(&x);
                    value.inverse().unwrap()
                })
                .collect::<Vec<F>>(),
        ));
    });

    helper_columns
}
