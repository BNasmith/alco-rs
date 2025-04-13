use super::*;
use rayon::prelude::*;
use octavian::Octavian;
use std::collections::HashSet;

#[test]
/// Ensure that the definition of the identity Octavian is correct.
fn test_one() {
    let one = Octavian::<i8>::one();
    assert_eq!(one, Octavian::new([ -2, -3, -4, -6, -5, -4, -3, -2 ]));
}

#[test]
/// Ensure that addition works.
fn test_addition() {
    let one = Octavian::<i8>::one();
    assert_eq!(one.clone() + one, Octavian::new([ -4, -6, -8, -12, -10, -8, -6, -4 ]));
}

#[test]
/// Ensure that subtraction works.
fn test_subtraction() {
    let one = Octavian::<i8>::one();
    assert_eq!(one.clone() - one, Octavian::new([ 0i8; 8]));
}

#[test]
/// Ensure that the inner product works.
fn test_inner_product() {
    let b = Octavian::<isize>::BASIS;
    assert_eq!(-1, b[1].inner_product(b[3]));
}

#[test]
/// Ensure that the norm works.
fn test_norm() {
    let b = Octavian::<isize>::BASIS;
    assert_eq!(2, b[1].norm());
}

#[test]
/// Ensure that the 240 Octavian units form a closed set under multiplication.
fn closure_of_units() {
    let mut units: HashSet<Octavian<i8>> = HashSet::new();
    for u in Octavian::<i8>::OCTAVIAN_UNITS_COEFFICIENTS {
        let x = Octavian::new(u);
        units.insert(x);
    }
    assert_eq!(240, units.len());
    let mut result = HashSet::<Octavian<i8>>::new();
    for u in &units {
        for v in &units {
            result.insert(*u * *v);
        }
    }
    assert_eq!(240, result.len())
}    

#[test]
fn closure_of_units_parallel() {
    let units: HashSet<Octavian<i8>> = Octavian::<i8>::OCTAVIAN_UNITS_COEFFICIENTS
        .iter()
        .map(|&u| Octavian::new(u))
        .collect();

    assert_eq!(240, units.len());

    let result: HashSet<Octavian<i8>> = units
        .par_iter()
        .flat_map(|u| {
            units.par_iter().map(move |v| u.clone() * v.clone())
        })
        .collect();

    assert_eq!(240, result.len());
}
