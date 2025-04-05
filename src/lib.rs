use std::ops::{Add, Mul, Div};
use std::iter::Sum;

/// A struct representing an octavian integer, provided that the coefficients are integer valued. Generalized to general type `T` to allow for general octonion calculations.
#[derive(Clone, Debug, PartialEq)]
pub struct Octavian<T> 
{
    /// Coefficients correspond to the basis vectors given in my thesis.
    pub coefficients: Vec<T>,
    pub norm: T,
    pub trace: T,
}

impl<T> Octavian<T>
where 
    T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    /// Creates a new Octavian with the given coefficients.
    pub fn new(coefficients: Vec<T>) -> Self {
        assert_eq!(coefficients.len(), 8);
        
        Self {
            // Specify the coefficients 
            coefficients: coefficients.clone(), 
            // Calculate the norm
            norm: Octavian::<T>::inner_product(
                coefficients.clone(), 
                coefficients.clone()
            )/2,
            // Calculate the trace
            trace: Octavian::<T>::inner_product(
                coefficients.clone(), 
                // Octavian::one().coefficients,
                [ -2, -3, -4, -6, -5, -4, -3, -2 ].map(|a| a.into()).to_vec(),
            ),
            
        }
    }

    /// Creates a new Octavian from an array of coefficients.
    pub fn from_array(coefficients: [T; 8]) -> Self {
        let vec: Vec<T> = coefficients.to_vec();
        Self {
            // Specify the coefficients 
            coefficients: vec.clone(), 
            // Calculate the norm
            norm: Octavian::<T>::inner_product(
                vec.clone(), 
                vec.clone()
            )/2,
            // Calculate the trace
            trace: Octavian::<T>::inner_product(
                vec.clone(), 
                // Octavian::one().coefficients,
                [ -2, -3, -4, -6, -5, -4, -3, -2 ].map(|a| a.into()).to_vec(),
            ),
            
        }
    }

    /// Creates the identity Octavian.
    pub fn one() -> Self {
        Octavian {
            coefficients: [ -2, -3, -4, -6, -5, -4, -3, -2 ]
                .map(|a| a.into())
                .to_vec(),
            norm: 1.into(),
            trace: 2.into(),
        }
    }

    /// Creates the zero Octavian.
    pub fn zero() -> Self {
        Octavian {
            coefficients: [ 0, 0, 0, 0, 0, 0, 0, 0 ]
                .map(|a| a.into())
                .to_vec(),
            norm: 0.into(),
            trace: 0.into(),
        }
    }

    /// Defines the inner product between the basis vectors.
    pub const GRAM_MATRIX: [[i8; 8]; 8] = 
  [ [   2,   0,  -1,   0,   0,   0,   0,   0 ],
    [   0,   2,   0,  -1,   0,   0,   0,   0 ],
    [  -1,   0,   2,  -1,   0,   0,   0,   0 ],
    [   0,  -1,  -1,   2,  -1,   0,   0,   0 ],
    [   0,   0,   0,  -1,   2,  -1,   0,   0 ],
    [   0,   0,   0,   0,  -1,   2,  -1,   0 ],
    [   0,   0,   0,   0,   0,  -1,   2,  -1 ],
    [   0,   0,   0,   0,   0,   0,  -1,   2 ] ];

    /// Takes the inner product given a pair of Octavian coefficients.
    pub fn inner_product(x: Vec<T>, y: Vec<T>) -> T {
        // Compute the norm using the Gram Matrix
        let mut result: T = 0.into();
        // Multiply x by the Gram matrix
        let gram = Octavian::<T>::GRAM_MATRIX;
        let a = x.clone();
        let b = y.clone();
        let mut intermediate: Vec<T> = vec![0.into(); 8];
        for i in 0..8 {
            for j in 0..8 {
                intermediate[i] = intermediate[i] + (a[j] * gram[j][i].into());
            }
        }
        // Take the dot product of the intermediate vector with y
        for i in 0..8 {
            result = result + (intermediate[i] * b[i]);
        }
        result
    }

    // pub const BASIS: [Octavian<i8>; 8] = [
    //     Octavian::from_array([1, 0, 0, 0, 0, 0, 0, 0]), 
    //     Octavian::from_array([0, 1, 0, 0, 0, 0, 0, 0]), 
    //     Octavian::from_array([0, 0, 1, 0, 0, 0, 0, 0]), 
    //     Octavian::from_array([0, 0, 0, 1, 0, 0, 0, 0]), 
    //     Octavian::from_array([0, 0, 0, 0, 1, 0, 0, 0]), 
    //     Octavian::from_array([0, 0, 0, 0, 0, 1, 0, 0]), 
    //     Octavian::from_array([0, 0, 0, 0, 0, 0, 1, 0]), 
    //     Octavian::from_array([0, 0, 0, 0, 0, 0, 0, 1])
    // ];

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

}

    ///     Implement right scalar multiplication for Octavian<T> where T is the scalar. 
    impl<T> Mul<T> for Octavian<T>
    where 
        T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
    {
        type Output = Octavian<T>;
        fn mul(self, rhs: T) -> Octavian<T> {
            let mut z: Vec<T> = Vec::new();
            let scalar = &rhs; 
            for a in self.clone().coefficients.into_iter() {
                z.push(a.mul(rhs.clone()));
            }
            Octavian { 
                coefficients: z, 
                norm: self.norm.clone()*scalar.clone()*scalar.clone(), 
                trace: self.trace.clone() * scalar.clone()
            }
        }
    }

/// Implements addition for Octavians.
impl<T> Add for Octavian<T>
where 
    T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut z: Vec<T> = Vec::new();
        for (a, b) in self.coefficients.clone().into_iter().zip(rhs.coefficients.clone()) {
            z.push(a.add(b));
        }
        Octavian { 
            coefficients: z.clone(),
            norm: Octavian::<T>::inner_product(z.clone(), z.clone())/2,
            trace: self.trace + rhs.trace 
        }
    }
}


/// Left adjoint matrix for an Octavian. 
pub fn left_adjoint_matrix<T>(oct: Octavian<T>) -> [[T; 8]; 8]
where 
    T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    let adj = Octavian::<T>::OCTAVIAN_ADJOINT_MATRICES.clone();
    let coeffs = oct.clone().coefficients;
    let mut result: [[T; 8]; 8] = [[0.into(); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            result[i][j] = coeffs
                .iter()
                .zip(adj.iter())
                .map(|(n, matrix)| *n * T::from(matrix[i][j]))
                .fold(T::default(), |acc, x| acc + x);
        }
    }
    result
}

// /// Transposed matrix. 
// fn transpose<T: Default + Clone>(matrix: [[T;8];8]) -> Vec<Vec<T>> {
//     let rows = matrix.len();
//     let cols = matrix[0].len();
//     let mut transposed = vec![vec![matrix[0][0].clone(); rows]; cols];

//     for i in 0..rows {
//         for j in 0..cols {
//             transposed[j][i] = matrix[i][j].clone();
//         }
//     }

//     transposed
// }

// /// Implement multiplication on Octavians using the OCTAVIAN_ADJOINT_MATRICES
impl<T> Mul for Octavian<T>
where 
    T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Octavian<T> {
        // Compute the new coefficients as the linear combination of the columns of the adjoint matrix. 
        let mut matrix: [[T; 8]; 8] = left_adjoint_matrix(self.clone());
        let coeffs = rhs.coefficients.clone();
        for i in 0..8 {
            for j in 0..8 {
                matrix[i][j] = matrix[i][j]*coeffs[j];
            }
        }
        // Next the result is a column vector taken by summing over the row entries.
        let mut result = Vec::new();
        for k in 0..8 {
            let mut entry: T = 0.into();
            for l in 0..8 {
                entry = entry + matrix[k][l];
            }
            result.push(entry);
        }
        // Construct the product. 
        Octavian {
            // Specify the coefficients
            coefficients: result.clone(),
            // Calculate the norm (using the composition rule)
            norm: self.norm * rhs.clone().norm,
            // Calculate the trace
            trace: Octavian::<T>::inner_product(
                result.clone(), 
                // Octavian::one().coefficients,
                [ -2, -3, -4, -6, -5, -4, -3, -2 ].map(|a| a.into()).to_vec(),
            ),
        }
    }
}



