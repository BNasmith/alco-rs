use num_traits::{ConstOne, ConstZero, Inv, Num, One, Pow, Zero};
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
use std::ops::{AddAssign, MulAssign};

/// The octavian integers are defined in Conway and Smith, and elsewhere. 
#[derive(PartialEq, Eq, Copy, Clone, Hash, Debug, Default)]
pub struct Octavian<T> {
/// The E8 lattice coordinates of the octavian (as defined in Nasmith 2023, pp. 89-90).  
pub coefficients: [T; 8],
// The norm of the octavian.
// pub norm: T,
// The trace of the octavian.
// pub trace: T,
}

impl<T> Octavian<T> {
/// Create a new `Octavian`.
pub const fn new(coefficients: [T; 8]) -> Self {
    Octavian {
        coefficients: coefficients,
        // norm: T::from(1),
        // trace: 1,
    }
}
}

impl<T: Clone + Copy + Num + std::iter::Sum + From<i8> + ConstZero + ConstOne> Octavian<T> {
/// Multiplies `self` by the scalar `t`.
#[inline]
pub fn scale(&self, t: T) -> Self {
    Self::new(self.coefficients.clone().map(|x| x * t.clone()))
}

/// Divides `self` by the scalar `t`.
pub fn unscale(&self, t: T) -> Self {
    Self::new(self.coefficients.clone().map(|x| x / t.clone()))
}

/// Defines the inner product between the basis vectors.
pub const GRAM_MATRIX: [[i8; 8]; 8] = [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2]
];

/// Returns the Gram matrix with elements converted to type `T`.
pub fn gram_matrix_typed() -> [[T; 8]; 8] {
    Self::GRAM_MATRIX.map(|row: [i8; 8] | row.map(|x: i8| T::from(x)))
}

/// Returns the inner product of two `Octavian` elements.
pub fn inner_product(&self, rhs: Self) -> T {
    let g: [[i8; 8]; 8] = Self::GRAM_MATRIX;
    let temp = g.iter().map(|row| {
        row.iter()
            .zip(&rhs.coefficients)
            .map(|(&x, &y)| T::from(x) * y)
            .sum()
    });
    self.coefficients
        .iter()
        .zip(temp)
        .map(|(&x, y)| x * y)
        .sum()
}

/// Returns the norm of an octavian.
pub fn norm(&self) -> T {
    self.inner_product(self.clone())
}

/// Returns the trace of an octavian.
pub fn trace(&self) -> T {
    self.inner_product(Octavian::<T>::one())
}

/// Returns the conjugate of an `Octavian` element, which is the trace of the element (multiplied by the identity) minus the element.
pub fn conj(self) -> Self {
    Octavian::<T>::one().scale(self.trace()) + self.scale((-1).into())
}


pub const OCTAVIAN_ADJOINT_MATRICES: [[[i8; 8]; 8]; 8] = 
[ [ [ 2, -1, -1, 0, 1, 0, -1, 0 ], [ 3, -1, -1, 0, 1, 0, -2, 1 ],
    [ 4, -2, -2, 0, 2, 0, -2, 1 ], [ 6, -2, -3, 0, 3, -1, -3, 2 ],
    [ 5, -1, -3, 0, 2, 0, -3, 2 ], [ 4, -1, -3, 1, 1, 0, -2, 1 ],
    [ 3, 0, -2, 0, 1, 0, -1, 0 ], [ 2, 0, -1, 0, 0, 0, 0, 0 ] ],
[ [ 1, 2, -2, 0, 0, 0, 0, 0 ], [ 1, 3, -2, 0, -1, 1, -1, 0 ],
    [ 2, 4, -3, 0, -1, 1, -1, 0 ], [ 2, 6, -4, 0, -2, 2, -2, 1 ],
    [ 1, 5, -3, 0, -2, 2, -1, 0 ], [ 1, 4, -2, 0, -2, 1, 0, 0 ],
    [ 0, 3, -1, 0, -1, 0, 0, 0 ], [ 0, 2, 0, -1, 0, 0, 0, 0 ] ],
[ [ -1, 2, 2, -2, 0, 0, 0, 0 ], [ -2, 2, 3, -3, 0, 1, 0, 0 ],
    [ -2, 3, 4, -4, 0, 1, 0, -1 ], [ -3, 4, 6, -6, 0, 2, 0, -1 ],
    [ -2, 3, 5, -5, 0, 1, 1, -1 ], [ -1, 2, 4, -4, 0, 1, 0, 0 ],
    [ -1, 1, 3, -2, -1, 1, 0, 0 ], [ -1, 0, 2, -1, 0, 0, 0, 0 ] ],
[ [ 0, -2, 0, 2, -2, 1, 0, 0 ], [ 0, -3, 0, 3, -2, 0, 1, -1 ],
    [ 0, -4, 0, 4, -3, 0, 1, 0 ], [ 0, -6, 0, 6, -4, 0, 1, -1 ],
    [ 0, -5, 0, 5, -3, 0, 0, 0 ], [ -1, -4, 0, 4, -2, 0, 0, 0 ],
    [ 0, -3, -1, 3, -1, 0, 0, 0 ], [ 0, -1, -1, 2, -1, 0, 0, 0 ] ],
[ [ -1, 0, 0, 0, 2, -2, 0, 0 ], [ -1, 1, 0, -1, 3, -3, 0, 1 ],
    [ -2, 1, 0, -1, 4, -3, -1, 1 ], [ -3, 2, 0, -2, 6, -5, 0, 1 ],
    [ -2, 2, 0, -2, 5, -4, 0, 0 ], [ -1, 2, 0, -2, 4, -3, 0, 0 ],
    [ -1, 1, 1, -2, 3, -2, 0, 0 ], [ 0, 0, 0, -1, 2, -1, 0, 0 ] ],
[ [ 0, 0, 0, -1, 0, 2, 0, -1 ], [ 0, -1, -1, 0, 0, 3, -1, -1 ],
    [ 0, -1, -1, 0, -1, 4, 0, -2 ], [ 1, -2, -2, 0, -1, 6, -1, -2 ],
    [ 0, -2, -1, 0, -1, 5, -1, -1 ], [ 0, -1, -1, 0, -1, 4, -1, -1 ],
    [ 0, 0, -1, 0, -1, 3, -1, 0 ], [ 0, 0, 0, 0, -1, 2, -1, 0 ] ],
[ [ 1, 0, 0, 0, 0, -2, 2, 0 ], [ 2, 1, 0, -1, 0, -2, 3, -1 ],
    [ 2, 1, 0, -1, 1, -4, 4, -1 ], [ 3, 2, 0, -1, 0, -5, 6, -2 ],
    [ 3, 1, -1, 0, 0, -4, 5, -2 ], [ 2, 0, 0, 0, 0, -3, 4, -2 ],
    [ 1, 0, 0, 0, 0, -2, 3, -2 ], [ 0, 0, 0, 0, 0, -1, 2, -1 ] ],
[ [ -1, 0, 0, 0, 0, 1, -2, 2 ], [ -1, -1, 0, 1, -1, 1, -2, 3 ],
    [ -1, 0, 0, 0, -1, 2, -3, 4 ], [ -2, -1, 1, 0, -1, 2, -4, 6 ],
    [ -2, 0, 1, 0, -1, 1, -3, 5 ], [ -1, 0, 0, 0, 0, 0, -2, 4 ],
    [ 0, 0, 0, 0, 0, 0, -2, 3 ], [ 0, 0, 0, 0, 0, 0, -1, 1 ] ] ];

pub fn octavian_adjoint_matrices_typed() -> [[[T; 8]; 8]; 8] {
    Self::OCTAVIAN_ADJOINT_MATRICES.map(|matrix| {
        matrix.map(|row| row.map(|x| T::from(x)))
    })
}

pub const OCTAVIAN_UNITS_COEFFICIENTS: [[i8; 8]; 240] = 
  [ [  -2,  -3,  -4,  -6,  -5,  -4,  -3,  -2 ],
[  -2,  -3,  -4,  -6,  -5,  -4,  -3,  -1 ],
[  -2,  -3,  -4,  -6,  -5,  -4,  -2,  -1 ],
[  -2,  -3,  -4,  -6,  -5,  -3,  -2,  -1 ],
[  -2,  -3,  -4,  -6,  -4,  -3,  -2,  -1 ],
[  -2,  -3,  -4,  -5,  -4,  -3,  -2,  -1 ],
[  -2,  -3,  -3,  -5,  -4,  -3,  -2,  -1 ],
[  -2,  -2,  -4,  -5,  -4,  -3,  -2,  -1 ],
[  -2,  -2,  -3,  -5,  -4,  -3,  -2,  -1 ],
[  -2,  -2,  -3,  -4,  -4,  -3,  -2,  -1 ],
[  -2,  -2,  -3,  -4,  -3,  -3,  -2,  -1 ],
[  -2,  -2,  -3,  -4,  -3,  -2,  -2,  -1 ],
[  -2,  -2,  -3,  -4,  -3,  -2,  -1,  -1 ],
[  -2,  -2,  -3,  -4,  -3,  -2,  -1,   0 ],
[  -1,  -3,  -3,  -5,  -4,  -3,  -2,  -1 ],
[  -1,  -2,  -3,  -5,  -4,  -3,  -2,  -1 ],
[  -1,  -2,  -3,  -4,  -4,  -3,  -2,  -1 ],
[  -1,  -2,  -3,  -4,  -3,  -3,  -2,  -1 ],
[  -1,  -2,  -3,  -4,  -3,  -2,  -2,  -1 ],
[  -1,  -2,  -3,  -4,  -3,  -2,  -1,  -1 ],
[  -1,  -2,  -3,  -4,  -3,  -2,  -1,   0 ],
[  -1,  -2,  -2,  -4,  -4,  -3,  -2,  -1 ],
[  -1,  -2,  -2,  -4,  -3,  -3,  -2,  -1 ],
[  -1,  -2,  -2,  -4,  -3,  -2,  -2,  -1 ],
[  -1,  -2,  -2,  -4,  -3,  -2,  -1,  -1 ],
[  -1,  -2,  -2,  -4,  -3,  -2,  -1,   0 ],
[  -1,  -2,  -2,  -3,  -3,  -3,  -2,  -1 ],
[  -1,  -2,  -2,  -3,  -3,  -2,  -2,  -1 ],
[  -1,  -2,  -2,  -3,  -3,  -2,  -1,  -1 ],
[  -1,  -2,  -2,  -3,  -3,  -2,  -1,   0 ],
[  -1,  -2,  -2,  -3,  -2,  -2,  -2,  -1 ],
[  -1,  -2,  -2,  -3,  -2,  -2,  -1,  -1 ],
[  -1,  -2,  -2,  -3,  -2,  -2,  -1,   0 ],
[  -1,  -2,  -2,  -3,  -2,  -1,  -1,  -1 ],
[  -1,  -2,  -2,  -3,  -2,  -1,  -1,   0 ],
[  -1,  -2,  -2,  -3,  -2,  -1,   0,   0 ],
[  -1,  -1,  -2,  -3,  -3,  -3,  -2,  -1 ],
[  -1,  -1,  -2,  -3,  -3,  -2,  -2,  -1 ],
[  -1,  -1,  -2,  -3,  -3,  -2,  -1,  -1 ],
[  -1,  -1,  -2,  -3,  -3,  -2,  -1,   0 ],
[  -1,  -1,  -2,  -3,  -2,  -2,  -2,  -1 ],
[  -1,  -1,  -2,  -3,  -2,  -2,  -1,  -1 ],
[  -1,  -1,  -2,  -3,  -2,  -2,  -1,   0 ],
[  -1,  -1,  -2,  -3,  -2,  -1,  -1,  -1 ],
[  -1,  -1,  -2,  -3,  -2,  -1,  -1,   0 ],
[  -1,  -1,  -2,  -3,  -2,  -1,   0,   0 ],
[  -1,  -1,  -2,  -2,  -2,  -2,  -2,  -1 ],
[  -1,  -1,  -2,  -2,  -2,  -2,  -1,  -1 ],
[  -1,  -1,  -2,  -2,  -2,  -2,  -1,   0 ],
[  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1 ],
[  -1,  -1,  -2,  -2,  -2,  -1,  -1,   0 ],
[  -1,  -1,  -2,  -2,  -2,  -1,   0,   0 ],
[  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1 ],
[  -1,  -1,  -2,  -2,  -1,  -1,  -1,   0 ],
[  -1,  -1,  -2,  -2,  -1,  -1,   0,   0 ],
[  -1,  -1,  -2,  -2,  -1,   0,   0,   0 ],
[  -1,  -1,  -1,  -2,  -2,  -2,  -2,  -1 ],
[  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1 ],
[  -1,  -1,  -1,  -2,  -2,  -2,  -1,   0 ],
[  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1 ],
[  -1,  -1,  -1,  -2,  -2,  -1,  -1,   0 ],
[  -1,  -1,  -1,  -2,  -2,  -1,   0,   0 ],
[  -1,  -1,  -1,  -2,  -1,  -1,  -1,  -1 ],
[  -1,  -1,  -1,  -2,  -1,  -1,  -1,   0 ],
[  -1,  -1,  -1,  -2,  -1,  -1,   0,   0 ],
[  -1,  -1,  -1,  -2,  -1,   0,   0,   0 ],
[  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1 ],
[  -1,  -1,  -1,  -1,  -1,  -1,  -1,   0 ],
[  -1,  -1,  -1,  -1,  -1,  -1,   0,   0 ],
[  -1,  -1,  -1,  -1,  -1,   0,   0,   0 ],
[  -1,  -1,  -1,  -1,   0,   0,   0,   0 ],
[  -1,   0,  -1,  -1,  -1,  -1,  -1,  -1 ],
[  -1,   0,  -1,  -1,  -1,  -1,  -1,   0 ],
[  -1,   0,  -1,  -1,  -1,  -1,   0,   0 ],
[  -1,   0,  -1,  -1,  -1,   0,   0,   0 ],
[  -1,   0,  -1,  -1,   0,   0,   0,   0 ],
[  -1,   0,  -1,   0,   0,   0,   0,   0 ],
[  -1,   0,   0,   0,   0,   0,   0,   0 ],
[   0,  -1,  -1,  -2,  -2,  -2,  -2,  -1 ],
[   0,  -1,  -1,  -2,  -2,  -2,  -1,  -1 ],
[   0,  -1,  -1,  -2,  -2,  -2,  -1,   0 ],
[   0,  -1,  -1,  -2,  -2,  -1,  -1,  -1 ],
[   0,  -1,  -1,  -2,  -2,  -1,  -1,   0 ],
[   0,  -1,  -1,  -2,  -2,  -1,   0,   0 ],
[   0,  -1,  -1,  -2,  -1,  -1,  -1,  -1 ],
[   0,  -1,  -1,  -2,  -1,  -1,  -1,   0 ],
[   0,  -1,  -1,  -2,  -1,  -1,   0,   0 ],
[   0,  -1,  -1,  -2,  -1,   0,   0,   0 ],
[   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1 ],
[   0,  -1,  -1,  -1,  -1,  -1,  -1,   0 ],
[   0,  -1,  -1,  -1,  -1,  -1,   0,   0 ],
[   0,  -1,  -1,  -1,  -1,   0,   0,   0 ],
[   0,  -1,  -1,  -1,   0,   0,   0,   0 ],
[   0,  -1,   0,  -1,  -1,  -1,  -1,  -1 ],
[   0,  -1,   0,  -1,  -1,  -1,  -1,   0 ],
[   0,  -1,   0,  -1,  -1,  -1,   0,   0 ],
[   0,  -1,   0,  -1,  -1,   0,   0,   0 ],
[   0,  -1,   0,  -1,   0,   0,   0,   0 ],
[   0,  -1,   0,   0,   0,   0,   0,   0 ],
[   0,   0,  -1,  -1,  -1,  -1,  -1,  -1 ],
[   0,   0,  -1,  -1,  -1,  -1,  -1,   0 ],
[   0,   0,  -1,  -1,  -1,  -1,   0,   0 ],
[   0,   0,  -1,  -1,  -1,   0,   0,   0 ],
[   0,   0,  -1,  -1,   0,   0,   0,   0 ],
[   0,   0,  -1,   0,   0,   0,   0,   0 ],
[   0,   0,   0,  -1,  -1,  -1,  -1,  -1 ],
[   0,   0,   0,  -1,  -1,  -1,  -1,   0 ],
[   0,   0,   0,  -1,  -1,  -1,   0,   0 ],
[   0,   0,   0,  -1,  -1,   0,   0,   0 ],
[   0,   0,   0,  -1,   0,   0,   0,   0 ],
[   0,   0,   0,   0,  -1,  -1,  -1,  -1 ],
[   0,   0,   0,   0,  -1,  -1,  -1,   0 ],
[   0,   0,   0,   0,  -1,  -1,   0,   0 ],
[   0,   0,   0,   0,  -1,   0,   0,   0 ],
[   0,   0,   0,   0,   0,  -1,  -1,  -1 ],
[   0,   0,   0,   0,   0,  -1,  -1,   0 ],
[   0,   0,   0,   0,   0,  -1,   0,   0 ],
[   0,   0,   0,   0,   0,   0,  -1,  -1 ],
[   0,   0,   0,   0,   0,   0,  -1,   0 ],
[   0,   0,   0,   0,   0,   0,   0,  -1 ],
[   0,   0,   0,   0,   0,   0,   0,   1 ],
[   0,   0,   0,   0,   0,   0,   1,   0 ],
[   0,   0,   0,   0,   0,   0,   1,   1 ],
[   0,   0,   0,   0,   0,   1,   0,   0 ],
[   0,   0,   0,   0,   0,   1,   1,   0 ],
[   0,   0,   0,   0,   0,   1,   1,   1 ],
[   0,   0,   0,   0,   1,   0,   0,   0 ],
[   0,   0,   0,   0,   1,   1,   0,   0 ],
[   0,   0,   0,   0,   1,   1,   1,   0 ],
[   0,   0,   0,   0,   1,   1,   1,   1 ],
[   0,   0,   0,   1,   0,   0,   0,   0 ],
[   0,   0,   0,   1,   1,   0,   0,   0 ],
[   0,   0,   0,   1,   1,   1,   0,   0 ],
[   0,   0,   0,   1,   1,   1,   1,   0 ],
[   0,   0,   0,   1,   1,   1,   1,   1 ],
[   0,   0,   1,   0,   0,   0,   0,   0 ],
[   0,   0,   1,   1,   0,   0,   0,   0 ],
[   0,   0,   1,   1,   1,   0,   0,   0 ],
[   0,   0,   1,   1,   1,   1,   0,   0 ],
[   0,   0,   1,   1,   1,   1,   1,   0 ],
[   0,   0,   1,   1,   1,   1,   1,   1 ],
[   0,   1,   0,   0,   0,   0,   0,   0 ],
[   0,   1,   0,   1,   0,   0,   0,   0 ],
[   0,   1,   0,   1,   1,   0,   0,   0 ],
[   0,   1,   0,   1,   1,   1,   0,   0 ],
[   0,   1,   0,   1,   1,   1,   1,   0 ],
[   0,   1,   0,   1,   1,   1,   1,   1 ],
[   0,   1,   1,   1,   0,   0,   0,   0 ],
[   0,   1,   1,   1,   1,   0,   0,   0 ],
[   0,   1,   1,   1,   1,   1,   0,   0 ],
[   0,   1,   1,   1,   1,   1,   1,   0 ],
[   0,   1,   1,   1,   1,   1,   1,   1 ],
[   0,   1,   1,   2,   1,   0,   0,   0 ],
[   0,   1,   1,   2,   1,   1,   0,   0 ],
[   0,   1,   1,   2,   1,   1,   1,   0 ],
[   0,   1,   1,   2,   1,   1,   1,   1 ],
[   0,   1,   1,   2,   2,   1,   0,   0 ],
[   0,   1,   1,   2,   2,   1,   1,   0 ],
[   0,   1,   1,   2,   2,   1,   1,   1 ],
[   0,   1,   1,   2,   2,   2,   1,   0 ],
[   0,   1,   1,   2,   2,   2,   1,   1 ],
[   0,   1,   1,   2,   2,   2,   2,   1 ],
[   1,   0,   0,   0,   0,   0,   0,   0 ],
[   1,   0,   1,   0,   0,   0,   0,   0 ],
[   1,   0,   1,   1,   0,   0,   0,   0 ],
[   1,   0,   1,   1,   1,   0,   0,   0 ],
[   1,   0,   1,   1,   1,   1,   0,   0 ],
[   1,   0,   1,   1,   1,   1,   1,   0 ],
[   1,   0,   1,   1,   1,   1,   1,   1 ],
[   1,   1,   1,   1,   0,   0,   0,   0 ],
[   1,   1,   1,   1,   1,   0,   0,   0 ],
[   1,   1,   1,   1,   1,   1,   0,   0 ],
[   1,   1,   1,   1,   1,   1,   1,   0 ],
[   1,   1,   1,   1,   1,   1,   1,   1 ],
[   1,   1,   1,   2,   1,   0,   0,   0 ],
[   1,   1,   1,   2,   1,   1,   0,   0 ],
[   1,   1,   1,   2,   1,   1,   1,   0 ],
[   1,   1,   1,   2,   1,   1,   1,   1 ],
[   1,   1,   1,   2,   2,   1,   0,   0 ],
[   1,   1,   1,   2,   2,   1,   1,   0 ],
[   1,   1,   1,   2,   2,   1,   1,   1 ],
[   1,   1,   1,   2,   2,   2,   1,   0 ],
[   1,   1,   1,   2,   2,   2,   1,   1 ],
[   1,   1,   1,   2,   2,   2,   2,   1 ],
[   1,   1,   2,   2,   1,   0,   0,   0 ],
[   1,   1,   2,   2,   1,   1,   0,   0 ],
[   1,   1,   2,   2,   1,   1,   1,   0 ],
[   1,   1,   2,   2,   1,   1,   1,   1 ],
[   1,   1,   2,   2,   2,   1,   0,   0 ],
[   1,   1,   2,   2,   2,   1,   1,   0 ],
[   1,   1,   2,   2,   2,   1,   1,   1 ],
[   1,   1,   2,   2,   2,   2,   1,   0 ],
[   1,   1,   2,   2,   2,   2,   1,   1 ],
[   1,   1,   2,   2,   2,   2,   2,   1 ],
[   1,   1,   2,   3,   2,   1,   0,   0 ],
[   1,   1,   2,   3,   2,   1,   1,   0 ],
[   1,   1,   2,   3,   2,   1,   1,   1 ],
[   1,   1,   2,   3,   2,   2,   1,   0 ],
[   1,   1,   2,   3,   2,   2,   1,   1 ],
[   1,   1,   2,   3,   2,   2,   2,   1 ],
[   1,   1,   2,   3,   3,   2,   1,   0 ],
[   1,   1,   2,   3,   3,   2,   1,   1 ],
[   1,   1,   2,   3,   3,   2,   2,   1 ],
[   1,   1,   2,   3,   3,   3,   2,   1 ],
[   1,   2,   2,   3,   2,   1,   0,   0 ],
[   1,   2,   2,   3,   2,   1,   1,   0 ],
[   1,   2,   2,   3,   2,   1,   1,   1 ],
[   1,   2,   2,   3,   2,   2,   1,   0 ],
[   1,   2,   2,   3,   2,   2,   1,   1 ],
[   1,   2,   2,   3,   2,   2,   2,   1 ],
[   1,   2,   2,   3,   3,   2,   1,   0 ],
[   1,   2,   2,   3,   3,   2,   1,   1 ],
[   1,   2,   2,   3,   3,   2,   2,   1 ],
[   1,   2,   2,   3,   3,   3,   2,   1 ],
[   1,   2,   2,   4,   3,   2,   1,   0 ],
[   1,   2,   2,   4,   3,   2,   1,   1 ],
[   1,   2,   2,   4,   3,   2,   2,   1 ],
[   1,   2,   2,   4,   3,   3,   2,   1 ],
[   1,   2,   2,   4,   4,   3,   2,   1 ],
[   1,   2,   3,   4,   3,   2,   1,   0 ],
[   1,   2,   3,   4,   3,   2,   1,   1 ],
[   1,   2,   3,   4,   3,   2,   2,   1 ],
[   1,   2,   3,   4,   3,   3,   2,   1 ],
[   1,   2,   3,   4,   4,   3,   2,   1 ],
[   1,   2,   3,   5,   4,   3,   2,   1 ],
[   1,   3,   3,   5,   4,   3,   2,   1 ],
[   2,   2,   3,   4,   3,   2,   1,   0 ],
[   2,   2,   3,   4,   3,   2,   1,   1 ],
[   2,   2,   3,   4,   3,   2,   2,   1 ],
[   2,   2,   3,   4,   3,   3,   2,   1 ],
[   2,   2,   3,   4,   4,   3,   2,   1 ],
[   2,   2,   3,   5,   4,   3,   2,   1 ],
[   2,   2,   4,   5,   4,   3,   2,   1 ],
[   2,   3,   3,   5,   4,   3,   2,   1 ],
[   2,   3,   4,   5,   4,   3,   2,   1 ],
[   2,   3,   4,   6,   4,   3,   2,   1 ],
[   2,   3,   4,   6,   5,   3,   2,   1 ],
[   2,   3,   4,   6,   5,   4,   2,   1 ],
[   2,   3,   4,   6,   5,   4,   3,   1 ],
[   2,   3,   4,   6,   5,   4,   3,   2 ] ];

pub const BASIS: [Octavian<T>; 8] = [
    Octavian::new([T::ONE, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO]), 
    Octavian::new([T::ZERO, T::ONE, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO]), 
    Octavian::new([T::ZERO, T::ZERO, T::ONE, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO]), 
    Octavian::new([T::ZERO, T::ZERO, T::ZERO, T::ONE, T::ZERO, T::ZERO, T::ZERO, T::ZERO]), 
    Octavian::new([T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ONE, T::ZERO, T::ZERO, T::ZERO]), 
    Octavian::new([T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ONE, T::ZERO, T::ZERO]), 
    Octavian::new([T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ONE, T::ZERO]), 
    Octavian::new([T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ZERO, T::ONE])
];

}



/// The negative of an `Octavian` element is the octavian with the opposite coefficients.
impl<T: Clone + Copy + Num + std::iter::Sum + From<i8> + Neg<Output = T> + ConstZero + ConstOne> Neg for Octavian<T> {
type Output = Self;
fn neg(self) -> Self::Output {
    self.scale((-1).into())
}
}

/// The zero `Octavian` is the octavian with all zero coefficients. 
impl<T: ConstZero + From<i8>> Octavian<T> {
/// A constant `Octavian` 0.
pub const ZERO: Self = Self::new([T::ZERO; 8]);

/// A method to select zero from any `Octavian`.
pub fn zero(self) -> Self {
    Octavian::<T>::ZERO
}

/// The constant multiplicative identity `Octavian`.
pub fn one() -> Self {
    Self::new([
        T::from(-2), T::from(-3), T::from(-4), T::from(-6),
        T::from(-5), T::from(-4), T::from(-3), T::from(-2)
    ])
}
}

/// Implements addition for `Octavian` elements, which is just the sum of the coefficients.
impl<T: Clone + Copy + Num> Add for Octavian<T>
{
type Output = Self;
fn add(self, other: Self) -> Self {
    Self::Output::new(
        self.coefficients
        .iter()
        .zip(other.coefficients.iter())
        .map(|(x, y)| x.clone() + y.clone())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap_or_else(|v: Vec<T>| panic!("Expected a Vec of length {} but it was {}", 8, v.len()))
    )
}
}

/// Implement right scalar multiplication on an Octavian<T> where T is the scalar. 
impl<T: Clone + Copy + Num + std::iter::Sum + From<i8> + ConstZero + ConstOne> Mul<T> for Octavian<T>
{
type Output = Octavian<T>;
fn mul(self, rhs: T) -> Octavian<T> {
    self.scale(rhs)
}
}

/// Implement right scalar division on an Octavian<T> where T is the scalar. 
impl<T: Clone + Copy + Num + std::iter::Sum + From<i8> + ConstZero + ConstOne> Div<T> for Octavian<T>
{
type Output = Octavian<T>;
fn div(self, rhs: T) -> Octavian<T> {
    self.unscale(rhs)
}
}

/// Implements subtraction for `Octavian` elements, which is just the difference of the coefficients.
impl<T: Clone + Copy + Num> Sub for Octavian<T>
{
type Output = Self;
fn sub(self, other: Self) -> Self {
    Self::Output::new(
        self.coefficients
        .iter()
        .zip(other.coefficients.iter())
        .map(|(x, y)| x.clone() - y.clone())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap_or_else(|v: Vec<T>| panic!("Expected a Vec of length {} but it was {}", 8, v.len()))
    )
}
}

// /// Multiply a matrix by a scalar. 
// fn scale_matrix<T>(matrix: [[T; 8]; 8], t: T) -> [[T; 8]; 8]
// where
//     T: Clone + Mul<Output = T>,
// {
//     matrix.map(|row| row.map(|x| x * t.clone()))
// }

impl<T: Clone + Copy + Num + std::iter::Sum + From<i8> + ConstZero + ConstOne> Octavian<T> {
/// Computes the left adjoint matrix of an `Octavian` element in the basis given by the coefficients.
pub fn left_adjoint_matrix(&self) -> [[T; 8]; 8] {
    // Get the typed adjoint matrices.
    let adj_matrices = Self::OCTAVIAN_ADJOINT_MATRICES;

    // Initialize a zero matrix.
    let mut result = [[T::zero(); 8]; 8];

    // Iterate over the adjoint matrices and coefficients.
    for (matrix, &coeff) in adj_matrices.iter().zip(&self.coefficients) {
        for (i, row) in matrix.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                result[i][j] = result[i][j] + T::from(value) * coeff;
            }
        }
    }

    result
}
}

/// Implements multiplication for `Octavian` elements.
impl<T: Clone + Copy + Num + std::iter::Sum + From<i8> + ConstZero + ConstOne> Mul for Octavian<T>
{
type Output = Self;
fn mul(self, other: Self) -> Self {
    // Compute the left adjoint matrix of `self`.
    let left_matrix = self.left_adjoint_matrix();

    // Multiply the matrix with the coefficients of `other`.
    let mut result_coefficients = [T::zero(); 8];
    for (i, row) in left_matrix.iter().enumerate() {
        result_coefficients[i] = row
            .iter()
            .zip(other.coefficients.iter())
            .map(|(a, b)| *a * *b)
            .sum();
    }

    // Convert the resulting coefficients back to an `Octavian`.
    Octavian::new(result_coefficients)
}
}

/// Implements AddAssign for `Octavian` elements.
impl<T: Clone + Copy + Num> AddAssign for Octavian<T> {
fn add_assign(&mut self, other: Self) {
    for (a, b) in self.coefficients.iter_mut().zip(other.coefficients.iter()) {
        *a = *a + *b;
    }
}
}

/// Implements MulAssign for `Octavian` elements.
impl<T: Clone + Copy + Num + std::iter::Sum + From<i8> + ConstZero + ConstOne> MulAssign for Octavian<T> {
fn mul_assign(&mut self, other: Self) {
    *self = *self * other;
}
}
