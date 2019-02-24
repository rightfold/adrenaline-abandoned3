//! Complex numbers and operations on complex numbers.

use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

/// A 128-bit complex number consists of a 64-bit real part and a 64-bit
/// imaginary part.
#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct c128(pub f64, pub f64);

impl c128 {
    /// The complex number with the given real part and a zero imaginary part.
    pub const fn from_real(real: f64) -> c128 {
        c128(real, 0.0)
    }

    /// The complex number with the given imaginary part and a zero real part.
    pub const fn from_imag(real: f64) -> c128 {
        c128(0.0, real)
    }

    /// The complex number at the given polar coordinates.
    #[inline(always)]
    pub fn from_polar(r: f64, th: f64) -> c128 {
        let (s, c) = th.sin_cos();
        c128(r * c, r * s)
    }

    /// The real part of the complex number.
    pub const fn real(self) -> f64 {
        self.0
    }

    /// The imaginary part of the complex number.
    pub const fn imag(self) -> f64 {
        self.1
    }

    /// The complex conjugate of the complex number.
    #[inline(always)]
    pub fn conj(self) -> c128 {
        c128(self.0, -self.1)
    }
}

impl Add for c128 {
    type Output = c128;

    #[inline(always)]
    fn add(self, rhs: c128) -> c128 {
        c128(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl Sub for c128 {
    type Output = c128;

    #[inline(always)]
    fn sub(self, rhs: c128) -> c128 {
        c128(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl Mul for c128 {
    type Output = c128;

    #[inline(always)]
    fn mul(self, rhs: c128) -> c128 {
        c128(self.0 * rhs.0 - self.1 * rhs.1,
             self.0 * rhs.1 + self.1 * rhs.0)
    }
}

impl Div for c128 {
    type Output = c128;

    #[inline(always)]
    fn div(self, rhs: c128) -> c128 {
        let num = self * rhs.conj();
        let den = rhs  * rhs.conj();
        c128(num.0 / den.0, num.1 / den.0)
    }
}
