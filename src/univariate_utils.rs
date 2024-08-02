use crate::field::{ChallengeField, Field};

// Standard lagrange interpolation, assuming indices for evals are 0, 1, 2, ...
// NOTE: can be sped up if we precompute lagrange coeffs for a given size
pub fn lagrange_interpolation<F: Field>(evals: &[F]) -> Vec<F> {
    let multiply_polys = |a: &[F], b: &[F]| -> Vec<F> {
        let mut result = vec![F::ZERO; a.len() + b.len() - 1];
        a.iter().enumerate().for_each(|(i, c1)| {
            b.iter().enumerate().for_each(|(j, c2)| {
                let mut m = c1.clone();
                m.mul_assign(c2);
                result[i + j].add_assign(&m);
            });
        });

        result
    };

    let mut polynomial = vec![F::ZERO; evals.len()];
    evals.iter().enumerate().for_each(|(i, eval)| {
        let mut coeffs = vec![F::ONE];
        evals.iter().enumerate().for_each(|(j, _)| {
            if i != j {
                let denom_inv = {
                    let mut i = F::from_usize(i);
                    let j = F::from_usize(j);
                    i.sub_assign(&j);
                    i.inverse().unwrap()
                };

                let mut new_coeffs = vec![F::ZERO; 2];
                new_coeffs[0] = F::from_usize(j);
                new_coeffs[0].negate();
                new_coeffs[0].mul_assign(&denom_inv);
                new_coeffs[1] = denom_inv;
                coeffs = multiply_polys(&coeffs, &new_coeffs);
            }
        });

        coeffs.iter().enumerate().for_each(|(k, c)| {
            let mut res = eval.clone();
            res.mul_assign(c);
            polynomial[k].add_assign(&res);
        });
    });

    polynomial
}

// Simple Horner's method evaluation.
pub fn univariate_eval<F: Field>(coeffs: &[F], point: F) -> F {
    coeffs
        .iter()
        .enumerate()
        .rev()
        .fold(F::ZERO, |mut acc, (i, coeff)| {
            acc.add_assign(coeff);
            if i != 0 {
                acc.mul_assign(&point);
            }
            acc
        })
}

// Horner's method evaluation that raises into an extension field.
pub fn univariate_eval_ext<F: Field, E: ChallengeField<F>>(coeffs: &[F], point: E) -> E {
    let mut init = point.clone();
    init.mul_base(&coeffs[coeffs.len() - 1]);

    coeffs
        .iter()
        .enumerate()
        .rev()
        .skip(1)
        .fold(init, |mut acc, (i, coeff)| {
            acc.add_base(coeff);
            if i != 0 {
                acc.mul_assign(&point);
            }
            acc
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::m31::{quartic::M31_4, M31};
    use rand::Rng;

    #[test]
    fn test_lagrange_2() {
        let mut evals = vec![M31::default(); 2];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = lagrange_interpolation(&evals);
        assert_eq!(evals[0], univariate_eval(&poly, M31::ZERO));
        assert_eq!(evals[1], univariate_eval(&poly, M31::ONE));
    }

    #[test]
    fn test_lagrange_3() {
        let mut evals = vec![M31::default(); 3];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = lagrange_interpolation(&evals);
        assert_eq!(evals[0], univariate_eval(&poly, M31::ZERO));
        assert_eq!(evals[1], univariate_eval(&poly, M31::ONE));
        assert_eq!(evals[2], univariate_eval(&poly, M31(2)));
    }

    #[test]
    fn test_lagrange_20() {
        let mut evals = vec![M31::default(); 20];
        evals
            .iter_mut()
            .for_each(|e| *e = M31(rand::thread_rng().gen_range(0..M31::ORDER)));
        let poly = lagrange_interpolation(&evals);
        for i in 0..20 {
            assert_eq!(evals[i], univariate_eval(&poly, M31::from_usize(i)));
        }
    }

    #[test]
    fn test_eval_2() {
        // f(x) = 5 + 2x
        let poly = vec![M31(5), M31(2)];

        // f(2) = 9
        assert_eq!(M31(9), univariate_eval(&poly, M31(2)));
    }

    #[test]
    fn test_eval_3() {
        // f(x) = 5 + 2x + 3x^2
        let poly = vec![M31(5), M31(2), M31(3)];

        // f(2) = 21
        assert_eq!(M31(21), univariate_eval(&poly, M31(2)));
    }

    #[test]
    fn test_eval_ext() {
        // f(x) = 5 + 2x
        let poly = vec![M31(5), M31(2)];

        // f([2, 0]) = 9
        assert_eq!(
            M31_4::from_usize(9),
            univariate_eval_ext(&poly, M31_4::from_usize(2))
        );
    }
}
