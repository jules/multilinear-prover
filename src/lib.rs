#![feature(generic_const_exprs)]
#![feature(associated_type_defaults)]

mod field;
mod m31;
mod mle;
mod pcs;
mod sumcheck;
mod transcript;

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
