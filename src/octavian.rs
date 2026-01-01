// use num::integer;
use core::ops::{Add, Mul, Neg, Sub};
use num_traits::{FromPrimitive, Num};
use std::fmt::Debug;

/// The octavian integers are defined in Conway and Smith's book, [On Quaternions and Octonions](https://www.routledge.com/On-Quaternions-and-Octonions/Conway-Smith/p/book/9781568811345), and elsewhere.
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T>,
{
    /// The E8 lattice coordinates of the octavian (as defined in [Nasmith 2023](https://espace.rmc.ca/jspui/handle/11264/1423), pp. 89-90).  
    pub coefficients: [T; 8],
}

impl<T> Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T> + From<i8>,
{
    /// Create a new `Octavian`.
    pub const fn new(coefficients: [T; 8]) -> Self {
        Octavian {
            coefficients,
        }
    }

    /// Returns the trace of an octavian.
    /// In the coordinates chosen, each component is trace-free except for the last one.
    pub fn trace(&self) -> T {
        -self.coefficients[7]
    }

    /// Returns the inner product of two octavians.
    pub fn inner_product(&self, rhs: &Octavian<T>) -> T {
        let x = self.coefficients;
        let y = rhs.coefficients;
        // From the Gram Matrix and Dynkin diagram for E8, addition and subtraction alternate to reduce chance of overflow.
        x[0] * y[0] + x[0] * y[0] - x[0] * y[2] - x[2] * y[0] + x[1] * y[1] + x[1] * y[1]
            - x[1] * y[3]
            - x[3] * y[1]
            + x[2] * y[2]
            + x[2] * y[2]
            - x[2] * y[3]
            - x[3] * y[2]
            + x[3] * y[3]
            + x[3] * y[3]
            - x[3] * y[4]
            - x[4] * y[3]
            + x[4] * y[4]
            + x[4] * y[4]
            - x[4] * y[5]
            - x[5] * y[4]
            + x[5] * y[5]
            + x[5] * y[5]
            - x[5] * y[6]
            - x[6] * y[5]
            + x[6] * y[6]
            + x[6] * y[6]
            - x[6] * y[7]
            - x[7] * y[6]
            + x[7] * y[7]
            + x[7] * y[7]
    }

    /// Returns the norm of an octavian scaled to the E8 lattice.
    /// Accordingly the norm is always an even number.  
    pub fn norm(&self) -> T {
        self.inner_product(self) / 2.into()
    }

    /// Multiplies `self` by the scalar `t`.
    pub fn scale(&self, t: T) -> Self {
        Self::new(self.coefficients.map(|x| x * t))
    }

    /// The constant multiplicative identity `Octavian`.
    pub fn one() -> Self {
        Self::new([2i8, 3, 4, 6, 5, 4, 3, 2].map(|x| -T::from(x)))
    }

    /// The constant multiplicative identity `Octavian`.
    pub fn zero() -> Self {
        Self::new([0i8, 0, 0, 0, 0, 0, 0, 0].map(|x| T::from(x)))
    }

    /// Conjugation of an octavian.
    /// Reverses the sign of the imaginary component.
    pub fn conjugate(&self) -> Self {
        Octavian::one().scale(self.trace()) - *self
    }

    /// Bimultiplication of octavians.
    /// B(x) = - Norm(x)*ref(1)*ref(x), where Norm(x) = 1 on the roots of E8.
    /// This simplifes to B(x)y = 2 <x,y> conjugate(x) - Norm(x) conjugate(y)
    /// Need to check this more closely.
    pub fn bimultiplication(self, other: Octavian<T>) -> Octavian<T> {
        other.conjugate() * self.inner_product(&other) - self.conjugate() * (self.norm())
    }

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

impl<T> Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T> + From<i8>,
{
    /// Defines the inner product between the basis vectors.
    pub const GRAM_MATRIX: [[i8; 8]; 8] = [
        [2, 0, -1, 0, 0, 0, 0, 0],
        [0, 2, 0, -1, 0, 0, 0, 0],
        [-1, 0, 2, -1, 0, 0, 0, 0],
        [0, -1, -1, 2, -1, 0, 0, 0],
        [0, 0, 0, -1, 2, -1, 0, 0],
        [0, 0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, 0, -1, 2, -1],
        [0, 0, 0, 0, 0, 0, -1, 2],
    ];

    pub const OCTAVIAN_ADJOINT_MATRICES: [[[i8; 8]; 8]; 8] = [
        [
            [2, -1, -1, 0, 1, 0, -1, 0],
            [3, -1, -1, 0, 1, 0, -2, 1],
            [4, -2, -2, 0, 2, 0, -2, 1],
            [6, -2, -3, 0, 3, -1, -3, 2],
            [5, -1, -3, 0, 2, 0, -3, 2],
            [4, -1, -3, 1, 1, 0, -2, 1],
            [3, 0, -2, 0, 1, 0, -1, 0],
            [2, 0, -1, 0, 0, 0, 0, 0],
        ],
        [
            [1, 2, -2, 0, 0, 0, 0, 0],
            [1, 3, -2, 0, -1, 1, -1, 0],
            [2, 4, -3, 0, -1, 1, -1, 0],
            [2, 6, -4, 0, -2, 2, -2, 1],
            [1, 5, -3, 0, -2, 2, -1, 0],
            [1, 4, -2, 0, -2, 1, 0, 0],
            [0, 3, -1, 0, -1, 0, 0, 0],
            [0, 2, 0, -1, 0, 0, 0, 0],
        ],
        [
            [-1, 2, 2, -2, 0, 0, 0, 0],
            [-2, 2, 3, -3, 0, 1, 0, 0],
            [-2, 3, 4, -4, 0, 1, 0, -1],
            [-3, 4, 6, -6, 0, 2, 0, -1],
            [-2, 3, 5, -5, 0, 1, 1, -1],
            [-1, 2, 4, -4, 0, 1, 0, 0],
            [-1, 1, 3, -2, -1, 1, 0, 0],
            [-1, 0, 2, -1, 0, 0, 0, 0],
        ],
        [
            [0, -2, 0, 2, -2, 1, 0, 0],
            [0, -3, 0, 3, -2, 0, 1, -1],
            [0, -4, 0, 4, -3, 0, 1, 0],
            [0, -6, 0, 6, -4, 0, 1, -1],
            [0, -5, 0, 5, -3, 0, 0, 0],
            [-1, -4, 0, 4, -2, 0, 0, 0],
            [0, -3, -1, 3, -1, 0, 0, 0],
            [0, -1, -1, 2, -1, 0, 0, 0],
        ],
        [
            [-1, 0, 0, 0, 2, -2, 0, 0],
            [-1, 1, 0, -1, 3, -3, 0, 1],
            [-2, 1, 0, -1, 4, -3, -1, 1],
            [-3, 2, 0, -2, 6, -5, 0, 1],
            [-2, 2, 0, -2, 5, -4, 0, 0],
            [-1, 2, 0, -2, 4, -3, 0, 0],
            [-1, 1, 1, -2, 3, -2, 0, 0],
            [0, 0, 0, -1, 2, -1, 0, 0],
        ],
        [
            [0, 0, 0, -1, 0, 2, 0, -1],
            [0, -1, -1, 0, 0, 3, -1, -1],
            [0, -1, -1, 0, -1, 4, 0, -2],
            [1, -2, -2, 0, -1, 6, -1, -2],
            [0, -2, -1, 0, -1, 5, -1, -1],
            [0, -1, -1, 0, -1, 4, -1, -1],
            [0, 0, -1, 0, -1, 3, -1, 0],
            [0, 0, 0, 0, -1, 2, -1, 0],
        ],
        [
            [1, 0, 0, 0, 0, -2, 2, 0],
            [2, 1, 0, -1, 0, -2, 3, -1],
            [2, 1, 0, -1, 1, -4, 4, -1],
            [3, 2, 0, -1, 0, -5, 6, -2],
            [3, 1, -1, 0, 0, -4, 5, -2],
            [2, 0, 0, 0, 0, -3, 4, -2],
            [1, 0, 0, 0, 0, -2, 3, -2],
            [0, 0, 0, 0, 0, -1, 2, -1],
        ],
        [
            [-1, 0, 0, 0, 0, 1, -2, 2],
            [-1, -1, 0, 1, -1, 1, -2, 3],
            [-1, 0, 0, 0, -1, 2, -3, 4],
            [-2, -1, 1, 0, -1, 2, -4, 6],
            [-2, 0, 1, 0, -1, 1, -3, 5],
            [-1, 0, 0, 0, 0, 0, -2, 4],
            [0, 0, 0, 0, 0, 0, -2, 3],
            [0, 0, 0, 0, 0, 0, -1, 1],
        ],
    ];

    pub const OCTAVIAN_UNITS_COEFFICIENTS: [[i8; 8]; 240] = [
        [-2, -3, -4, -6, -5, -4, -3, -2],
        [-2, -3, -4, -6, -5, -4, -3, -1],
        [-2, -3, -4, -6, -5, -4, -2, -1],
        [-2, -3, -4, -6, -5, -3, -2, -1],
        [-2, -3, -4, -6, -4, -3, -2, -1],
        [-2, -3, -4, -5, -4, -3, -2, -1],
        [-2, -3, -3, -5, -4, -3, -2, -1],
        [-2, -2, -4, -5, -4, -3, -2, -1],
        [-2, -2, -3, -5, -4, -3, -2, -1],
        [-2, -2, -3, -4, -4, -3, -2, -1],
        [-2, -2, -3, -4, -3, -3, -2, -1],
        [-2, -2, -3, -4, -3, -2, -2, -1],
        [-2, -2, -3, -4, -3, -2, -1, -1],
        [-2, -2, -3, -4, -3, -2, -1, 0],
        [-1, -3, -3, -5, -4, -3, -2, -1],
        [-1, -2, -3, -5, -4, -3, -2, -1],
        [-1, -2, -3, -4, -4, -3, -2, -1],
        [-1, -2, -3, -4, -3, -3, -2, -1],
        [-1, -2, -3, -4, -3, -2, -2, -1],
        [-1, -2, -3, -4, -3, -2, -1, -1],
        [-1, -2, -3, -4, -3, -2, -1, 0],
        [-1, -2, -2, -4, -4, -3, -2, -1],
        [-1, -2, -2, -4, -3, -3, -2, -1],
        [-1, -2, -2, -4, -3, -2, -2, -1],
        [-1, -2, -2, -4, -3, -2, -1, -1],
        [-1, -2, -2, -4, -3, -2, -1, 0],
        [-1, -2, -2, -3, -3, -3, -2, -1],
        [-1, -2, -2, -3, -3, -2, -2, -1],
        [-1, -2, -2, -3, -3, -2, -1, -1],
        [-1, -2, -2, -3, -3, -2, -1, 0],
        [-1, -2, -2, -3, -2, -2, -2, -1],
        [-1, -2, -2, -3, -2, -2, -1, -1],
        [-1, -2, -2, -3, -2, -2, -1, 0],
        [-1, -2, -2, -3, -2, -1, -1, -1],
        [-1, -2, -2, -3, -2, -1, -1, 0],
        [-1, -2, -2, -3, -2, -1, 0, 0],
        [-1, -1, -2, -3, -3, -3, -2, -1],
        [-1, -1, -2, -3, -3, -2, -2, -1],
        [-1, -1, -2, -3, -3, -2, -1, -1],
        [-1, -1, -2, -3, -3, -2, -1, 0],
        [-1, -1, -2, -3, -2, -2, -2, -1],
        [-1, -1, -2, -3, -2, -2, -1, -1],
        [-1, -1, -2, -3, -2, -2, -1, 0],
        [-1, -1, -2, -3, -2, -1, -1, -1],
        [-1, -1, -2, -3, -2, -1, -1, 0],
        [-1, -1, -2, -3, -2, -1, 0, 0],
        [-1, -1, -2, -2, -2, -2, -2, -1],
        [-1, -1, -2, -2, -2, -2, -1, -1],
        [-1, -1, -2, -2, -2, -2, -1, 0],
        [-1, -1, -2, -2, -2, -1, -1, -1],
        [-1, -1, -2, -2, -2, -1, -1, 0],
        [-1, -1, -2, -2, -2, -1, 0, 0],
        [-1, -1, -2, -2, -1, -1, -1, -1],
        [-1, -1, -2, -2, -1, -1, -1, 0],
        [-1, -1, -2, -2, -1, -1, 0, 0],
        [-1, -1, -2, -2, -1, 0, 0, 0],
        [-1, -1, -1, -2, -2, -2, -2, -1],
        [-1, -1, -1, -2, -2, -2, -1, -1],
        [-1, -1, -1, -2, -2, -2, -1, 0],
        [-1, -1, -1, -2, -2, -1, -1, -1],
        [-1, -1, -1, -2, -2, -1, -1, 0],
        [-1, -1, -1, -2, -2, -1, 0, 0],
        [-1, -1, -1, -2, -1, -1, -1, -1],
        [-1, -1, -1, -2, -1, -1, -1, 0],
        [-1, -1, -1, -2, -1, -1, 0, 0],
        [-1, -1, -1, -2, -1, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1],
        [-1, -1, -1, -1, -1, -1, -1, 0],
        [-1, -1, -1, -1, -1, -1, 0, 0],
        [-1, -1, -1, -1, -1, 0, 0, 0],
        [-1, -1, -1, -1, 0, 0, 0, 0],
        [-1, 0, -1, -1, -1, -1, -1, -1],
        [-1, 0, -1, -1, -1, -1, -1, 0],
        [-1, 0, -1, -1, -1, -1, 0, 0],
        [-1, 0, -1, -1, -1, 0, 0, 0],
        [-1, 0, -1, -1, 0, 0, 0, 0],
        [-1, 0, -1, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, -1, -2, -2, -2, -2, -1],
        [0, -1, -1, -2, -2, -2, -1, -1],
        [0, -1, -1, -2, -2, -2, -1, 0],
        [0, -1, -1, -2, -2, -1, -1, -1],
        [0, -1, -1, -2, -2, -1, -1, 0],
        [0, -1, -1, -2, -2, -1, 0, 0],
        [0, -1, -1, -2, -1, -1, -1, -1],
        [0, -1, -1, -2, -1, -1, -1, 0],
        [0, -1, -1, -2, -1, -1, 0, 0],
        [0, -1, -1, -2, -1, 0, 0, 0],
        [0, -1, -1, -1, -1, -1, -1, -1],
        [0, -1, -1, -1, -1, -1, -1, 0],
        [0, -1, -1, -1, -1, -1, 0, 0],
        [0, -1, -1, -1, -1, 0, 0, 0],
        [0, -1, -1, -1, 0, 0, 0, 0],
        [0, -1, 0, -1, -1, -1, -1, -1],
        [0, -1, 0, -1, -1, -1, -1, 0],
        [0, -1, 0, -1, -1, -1, 0, 0],
        [0, -1, 0, -1, -1, 0, 0, 0],
        [0, -1, 0, -1, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, -1, -1, -1, -1, -1],
        [0, 0, -1, -1, -1, -1, -1, 0],
        [0, 0, -1, -1, -1, -1, 0, 0],
        [0, 0, -1, -1, -1, 0, 0, 0],
        [0, 0, -1, -1, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0, 0, 0],
        [0, 0, 0, -1, -1, -1, -1, -1],
        [0, 0, 0, -1, -1, -1, -1, 0],
        [0, 0, 0, -1, -1, -1, 0, 0],
        [0, 0, 0, -1, -1, 0, 0, 0],
        [0, 0, 0, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, -1, -1, -1],
        [0, 0, 0, 0, -1, -1, -1, 0],
        [0, 0, 0, 0, -1, -1, 0, 0],
        [0, 0, 0, 0, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, -1, -1],
        [0, 0, 0, 0, 0, -1, -1, 0],
        [0, 0, 0, 0, 0, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, -1],
        [0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, -1],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 1, 1, 1],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 1, 1, 0],
        [0, 0, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, 1, 1, 0, 0],
        [0, 0, 0, 1, 1, 1, 1, 0],
        [0, 0, 0, 1, 1, 1, 1, 1],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 0, 0],
        [0, 0, 1, 1, 1, 1, 1, 0],
        [0, 0, 1, 1, 1, 1, 1, 1],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 0, 0, 0, 0],
        [0, 1, 0, 1, 1, 0, 0, 0],
        [0, 1, 0, 1, 1, 1, 0, 0],
        [0, 1, 0, 1, 1, 1, 1, 0],
        [0, 1, 0, 1, 1, 1, 1, 1],
        [0, 1, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 0],
        [0, 1, 1, 1, 1, 1, 1, 1],
        [0, 1, 1, 2, 1, 0, 0, 0],
        [0, 1, 1, 2, 1, 1, 0, 0],
        [0, 1, 1, 2, 1, 1, 1, 0],
        [0, 1, 1, 2, 1, 1, 1, 1],
        [0, 1, 1, 2, 2, 1, 0, 0],
        [0, 1, 1, 2, 2, 1, 1, 0],
        [0, 1, 1, 2, 2, 1, 1, 1],
        [0, 1, 1, 2, 2, 2, 1, 0],
        [0, 1, 1, 2, 2, 2, 1, 1],
        [0, 1, 1, 2, 2, 2, 2, 1],
        [1, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 1, 1, 0, 0, 0, 0],
        [1, 0, 1, 1, 1, 0, 0, 0],
        [1, 0, 1, 1, 1, 1, 0, 0],
        [1, 0, 1, 1, 1, 1, 1, 0],
        [1, 0, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 2, 1, 0, 0, 0],
        [1, 1, 1, 2, 1, 1, 0, 0],
        [1, 1, 1, 2, 1, 1, 1, 0],
        [1, 1, 1, 2, 1, 1, 1, 1],
        [1, 1, 1, 2, 2, 1, 0, 0],
        [1, 1, 1, 2, 2, 1, 1, 0],
        [1, 1, 1, 2, 2, 1, 1, 1],
        [1, 1, 1, 2, 2, 2, 1, 0],
        [1, 1, 1, 2, 2, 2, 1, 1],
        [1, 1, 1, 2, 2, 2, 2, 1],
        [1, 1, 2, 2, 1, 0, 0, 0],
        [1, 1, 2, 2, 1, 1, 0, 0],
        [1, 1, 2, 2, 1, 1, 1, 0],
        [1, 1, 2, 2, 1, 1, 1, 1],
        [1, 1, 2, 2, 2, 1, 0, 0],
        [1, 1, 2, 2, 2, 1, 1, 0],
        [1, 1, 2, 2, 2, 1, 1, 1],
        [1, 1, 2, 2, 2, 2, 1, 0],
        [1, 1, 2, 2, 2, 2, 1, 1],
        [1, 1, 2, 2, 2, 2, 2, 1],
        [1, 1, 2, 3, 2, 1, 0, 0],
        [1, 1, 2, 3, 2, 1, 1, 0],
        [1, 1, 2, 3, 2, 1, 1, 1],
        [1, 1, 2, 3, 2, 2, 1, 0],
        [1, 1, 2, 3, 2, 2, 1, 1],
        [1, 1, 2, 3, 2, 2, 2, 1],
        [1, 1, 2, 3, 3, 2, 1, 0],
        [1, 1, 2, 3, 3, 2, 1, 1],
        [1, 1, 2, 3, 3, 2, 2, 1],
        [1, 1, 2, 3, 3, 3, 2, 1],
        [1, 2, 2, 3, 2, 1, 0, 0],
        [1, 2, 2, 3, 2, 1, 1, 0],
        [1, 2, 2, 3, 2, 1, 1, 1],
        [1, 2, 2, 3, 2, 2, 1, 0],
        [1, 2, 2, 3, 2, 2, 1, 1],
        [1, 2, 2, 3, 2, 2, 2, 1],
        [1, 2, 2, 3, 3, 2, 1, 0],
        [1, 2, 2, 3, 3, 2, 1, 1],
        [1, 2, 2, 3, 3, 2, 2, 1],
        [1, 2, 2, 3, 3, 3, 2, 1],
        [1, 2, 2, 4, 3, 2, 1, 0],
        [1, 2, 2, 4, 3, 2, 1, 1],
        [1, 2, 2, 4, 3, 2, 2, 1],
        [1, 2, 2, 4, 3, 3, 2, 1],
        [1, 2, 2, 4, 4, 3, 2, 1],
        [1, 2, 3, 4, 3, 2, 1, 0],
        [1, 2, 3, 4, 3, 2, 1, 1],
        [1, 2, 3, 4, 3, 2, 2, 1],
        [1, 2, 3, 4, 3, 3, 2, 1],
        [1, 2, 3, 4, 4, 3, 2, 1],
        [1, 2, 3, 5, 4, 3, 2, 1],
        [1, 3, 3, 5, 4, 3, 2, 1],
        [2, 2, 3, 4, 3, 2, 1, 0],
        [2, 2, 3, 4, 3, 2, 1, 1],
        [2, 2, 3, 4, 3, 2, 2, 1],
        [2, 2, 3, 4, 3, 3, 2, 1],
        [2, 2, 3, 4, 4, 3, 2, 1],
        [2, 2, 3, 5, 4, 3, 2, 1],
        [2, 2, 4, 5, 4, 3, 2, 1],
        [2, 3, 3, 5, 4, 3, 2, 1],
        [2, 3, 4, 5, 4, 3, 2, 1],
        [2, 3, 4, 6, 4, 3, 2, 1],
        [2, 3, 4, 6, 5, 3, 2, 1],
        [2, 3, 4, 6, 5, 4, 2, 1],
        [2, 3, 4, 6, 5, 4, 3, 1],
        [2, 3, 4, 6, 5, 4, 3, 2],
    ];

    /// The unit octavians as an array in a canonical order.
    pub fn unit_vectors() -> [Self; 240] {
        Octavian::<T>::OCTAVIAN_UNITS_COEFFICIENTS
            .map(|coeffs| Octavian::new(coeffs.map(|x| x.into())))
    }

    /// The standard basis vectors for the octavian integers.
    pub fn basis_vectors() -> [Self; 8] {
        [
            Octavian::new([1i8, 0, 0, 0, 0, 0, 0, 0].map(|x| x.into())),
            Octavian::new([0i8, 1, 0, 0, 0, 0, 0, 0].map(|x| x.into())),
            Octavian::new([0i8, 0, 1, 0, 0, 0, 0, 0].map(|x| x.into())),
            Octavian::new([0i8, 0, 0, 1, 0, 0, 0, 0].map(|x| x.into())),
            Octavian::new([0i8, 0, 0, 0, 1, 0, 0, 0].map(|x| x.into())),
            Octavian::new([0i8, 0, 0, 0, 0, 1, 0, 0].map(|x| x.into())),
            Octavian::new([0i8, 0, 0, 0, 0, 0, 1, 0].map(|x| x.into())),
            Octavian::new([0i8, 0, 0, 0, 0, 0, 0, 1].map(|x| x.into())),
        ]
    }
}

/// Implements addition for `Octavian` elements, which is just the sum of the coefficients.
impl<T: Add<Output = T>> Add for Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T> + From<i8>,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let mut x = self.coefficients;
        let y = other.coefficients;
        for i in 0..8 {
            x[i] = x[i] + y[i];
        }
        Octavian::new(x)
    }
}

/// Implements subtraction for `Octavian` elements, which is just the difference of the coefficients.
impl<T: Sub<Output = T>> Sub for Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T> + From<i8>,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        let mut x = self.coefficients;
        let y = other.coefficients;
        for i in 0..8 {
            x[i] = x[i] - y[i];
        }
        Octavian::new(x)
    }
}

/// Implements negation for `Octavian` elements, which is just the negative of the coefficients.
impl<T: Neg<Output = T>> Neg for Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T> + From<i8>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut x = self.coefficients;
        // for i in 0..8 {
        //     x[i] = -x[i]
        // }
        for c in &mut x {
            *c = -*c;
        }
        Octavian::new(x)
    }
}

/// Implement right scalar multiplication on an `Octavian<T>` where `T` is the scalar.
impl<T: Mul<Output = T>> Mul<T> for Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T> + From<i8>,
{
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        self.scale(rhs)
    }
}

/// Implements multiplication for `Octavian` elements.
impl<T: Mul<Output = T>> Mul for Octavian<T>
where
    T: FromPrimitive + Num + Copy + Neg<Output = T> + From<i8>,
{
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        // Compute the left adjoint matrix of `self`.
        let left_matrix = self.left_adjoint_matrix();
        let mut coefficients = [T::zero(); 8];
        for i in 0..8 {
            for j in 0..8 {
                coefficients[i] = coefficients[i] + left_matrix[i][j] * other.coefficients[j];
            }
        }
        Self::new(coefficients)
    }
}
