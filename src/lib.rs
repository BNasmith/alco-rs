use std::ops::{Add, Sub, Mul, Div};
use std::iter::{Sum, FromIterator};
use std::hash::{Hash, Hasher};
// use num::Integer;


/// A struct representing an octavian integer, provided that the coefficients are integer valued. Generalized to general type `T` to allow for general octonion calculations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Octavian<T> 
{
    /// Coefficients correspond to the basis vectors given in my thesis.
    pub coefficients: Vec<T>,
    // Records the norm of the octavian.
    pub norm: T,
    // Records the trace of the octavian, which is twice the real component. 
    pub trace: T,
}

impl<T> Octavian<T>
where 
    T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    /// Creates a new Octavian with the given coefficients. The norm and trace are also computed and stored. 
    pub fn new(coefficients: Vec<T>) -> Self {
        assert_eq!(coefficients.len(), 8);
        
        Self {
            // Specify the coefficients. 
            coefficients: coefficients.clone(), 
            // Calculate the norm.
            norm: Octavian::<T>::inner_product(
                coefficients.clone(), 
                coefficients.clone()
            )/2,
            // Calculate the trace.
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

    const OCTAVIAN_ADJOINT_MATRICES: [[[i8; 8]; 8]; 8] = 
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

impl<T> Octavian<T>
where 
    T: Copy + Default + Clone + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>, 
{
/// Implement octonion conjugation, which has the formula x.conj() = Trace(x)*x.one() - x.
/// ```
/// use alco::Octavian;
/// let x = Octavian::<i8>::from_array([1,2,0,0,0,0,0,1]);
/// let one = Octavian::<i8>::one();
/// assert_eq!(one * x.norm, x.clone()*x.conj());
/// ```
    pub fn conj(self) -> Self {
        let result = (Octavian::<T>::one() * self.trace) - self;
        result
    }

/// Compute the norm of an Octavian directly.
/// ```
/// use alco::Octavian;
/// let x = Octavian::<i8>::from_array([0,0,0,3,0,0,0,1]);
/// let one = Octavian::<i8>::one();
/// assert_eq!(one * x.norm, x.clone()*x.conj());
/// ```
    pub fn check_norm(self) -> T {
        let result = self.clone()*self.conj();
        result.trace/2
    }

/// Compute the trace of an Octavian directly.
/// ```
/// use alco::Octavian;
/// let x = Octavian::<i8>::from_array([0,0,0,3,0,0,0,1]);
/// let one = Octavian::<i8>::one();
/// assert_eq!(one * x.trace, x.clone() + x.conj());
/// ```
    pub fn check_trace(self) -> T {
        let result = self.clone() + self.conj();
        result.trace
    }
}



/// Implement right scalar multiplication on an Octavian<T> where T is the scalar. 
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

/// Implements subtraction for Octavians.
impl<T> Sub for Octavian<T>
where 
    T: Copy + Default + Clone + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut z: Vec<T> = Vec::new();
        for (a, b) in self.coefficients.clone().into_iter().zip(rhs.coefficients.clone()) {
            z.push(a.sub(b));
        }
        Octavian { 
            coefficients: z.clone(),
            norm: Octavian::<T>::inner_product(z.clone(), z.clone())/2,
            trace: self.trace - rhs.trace 
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

/// Implement multiplication on Octavians using the `OCTAVIAN_ADJOINT_MATRICES` constant.
/// ```
/// use::alco::Octavian;
/// let x = Octavian::from_array([0,0,0,0,0,0,0,1]);
/// let y = Octavian::from_array([0,0,0,0,0,0,1,0]);
/// let z = x * y;
/// assert_eq!( z.clone(), Octavian::from_array([-2, -2, -3, -4, -3, -2, -2, -1]));
/// assert_eq!(z.clone().trace, 1);
/// assert_eq!(z.clone().norm, 1);
/// ``` 
impl<T> Mul for Octavian<T>
where 
    T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Octavian<T> {
        let mut result = Vec::new(); 
        for row in left_adjoint_matrix(self.clone()) {
            result.push(
                row.into_iter()
                .zip(rhs.coefficients.clone())
                .map(|(a,b)| a * b)
                .sum()
            )
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

impl<T> FromIterator<T> for Octavian<T>
where
    T: Copy + Default + Clone + Add<Output = T> + Mul<Output = T> + From<i8> + Div<i8, Output = T> + Sum<T>,
{
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let coefficients: Vec<T> = iter.into_iter().collect();
        assert_eq!(coefficients.len(), 8, "Octavian requires exactly 8 coefficients.");
        Octavian::new(coefficients)
    }
}

impl<T> Hash for Octavian<T>
where
    T: Hash + Eq,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.coefficients.hash(state);
        self.norm.hash(state);
        self.trace.hash(state);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    
    use std::collections::HashSet;

    #[test]
    /// Ensure that the definition of the identity Octavian is correct.
    fn test_one() {
        let one = Octavian::<i8>::one();
        assert_eq!(one, Octavian::new(vec![ -2, -3, -4, -6, -5, -4, -3, -2 ]));
    }

    #[test]
    /// Ensure that addition works.
    fn test_addition() {
        let one = Octavian::<i8>::one();
        assert_eq!(one.clone() + one, Octavian::new(vec![ -4, -6, -8, -12, -10, -8, -6, -4 ]));
    }

    #[test]
    /// Ensure that subtraction works.
    fn test_subtraction() {
        let one = Octavian::<i8>::one();
        assert_eq!(one.clone() - one, Octavian::new(vec![ 0i8; 8]));
    }

    #[test]
    /// Ensure that the 240 Octavian units form a closed set under multiplication.
    fn closure_of_units() {
        let mut units: HashSet<Octavian<i8>> = HashSet::new();
        for u in Octavian::<i8>::OCTAVIAN_UNITS_COEFFICIENTS {
            let x = Octavian::from_array(u);
            units.insert(x);
        }
        assert_eq!(240, units.len());
        let mut result = HashSet::<Octavian<i8>>::new();
        for u in &units {
            for v in &units {
                result.insert(u.clone() * v.clone());    
            }
        }
        assert_eq!(240, result.len())
    }    
}


// #[test]
