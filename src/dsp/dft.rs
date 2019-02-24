//! Discrete Fourier transform subroutines using the Cooley&ndash;Tukey fast
//! Fourier transform algorithm.

use std::f64::consts::PI;

use dsp::complex::c128;

/// Compute the forward discrete Fourier transform of the input.
///
/// When calling this subroutine, you must beware of certain restrictions and
/// liberties:
///
///  - The input slice must have a length _n_ &ge; 1.
///  - The output slice must have at least _n_ elements.
///  - The first _n_ elements of the output slice will be overwritten.
///  - It is safe to call this subroutine with the elements of the output slice
///    uninitialized.
///  - This subroutine does not have any side-effects other than overwriting
///    the elements of the output slice.
#[inline(always)]
pub fn fdft(input: &[c128], output: &mut [c128]) {
    assert!( !input.is_empty()           , "The input slice is empty"      );
    assert!( output.len() >= input.len() , "The output slice is too small" );
    unsafe { fft(input, output, input.len(), 1, |c| c); }
}

/// Compute the inverse discrete Fourier transform of the input.
///
/// The same restrictions and liberties apply as those to the [fdft]
/// subroutine.
///
/// [fdft]: fn.fdft.html
#[inline(always)]
pub fn idft(input: &[c128], output: &mut [c128]) {
    assert!( !input.is_empty()           , "The input slice is empty"      );
    assert!( output.len() >= input.len() , "The output slice is too small" );
    unsafe { fft(input, output, input.len(), 1, |c| c.conj()); }
    for r in &mut output[.. input.len()] {
        *r = r.conj() / c128::from_real(input.len() as f64);
    }
}

unsafe fn fft<F>(i: &[c128], o: &mut [c128], n: usize, s: usize, f: F)
    where F: Copy + Fn(c128) -> c128 {
    macro_rules! i { [$offset:expr] => { *i.get_unchecked    ($offset) }; }
    macro_rules! o { [$offset:expr] => { *o.get_unchecked_mut($offset) }; }

    if n == 1 {
        o![0] = f(i![0]);
        return;
    }

    fft(&i![ ..], &mut o![   ..], n/2, 2*s, f);
    fft(&i![s..], &mut o![n/2..], n/2, 2*s, f);

    for k in 0..n/2 {
        let (kf, nf) = (k as f64, n as f64);
        let tf = c128::from_polar(1.0, -2.0*PI*kf/nf) * o![k+n/2];
        let ok = o![k];
        o![k    ] = ok+tf;
        o![k+n/2] = ok-tf;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_aq {
        ($a:expr, $b:expr) => {{
            let (a, b) = ($a, $b);
            assert!(f64::abs(a.real() - b.real()) <= 0.01, "{:?} ≉ {:?}", a, b);
            assert!(f64::abs(a.imag() - b.imag()) <= 0.01, "{:?} ≉ {:?}", a, b);
        }};
    }

    #[test]
    fn test_fdft() {
        let     input  = [c128(1.00,  0.00), c128(1.00,  0.00),
                          c128(1.00,  0.00), c128(1.00,  0.00),
                          c128(0.00,  0.00), c128(0.00,  0.00),
                          c128(0.00,  0.00), c128(0.00,  0.00)];
        let mut output = [c128(0.0, 0.0); 8];
        fdft(&input, &mut output);
        assert_aq!(output[0], c128( 4.00 ,  0.00 ));
        assert_aq!(output[1], c128( 1.00 , -2.41 ));
        assert_aq!(output[2], c128( 0.00 ,  0.00 ));
        assert_aq!(output[3], c128( 1.00 , -0.41 ));
        assert_aq!(output[4], c128( 0.00 ,  0.00 ));
        assert_aq!(output[5], c128( 0.99 ,  0.41 ));
        assert_aq!(output[6], c128( 0.00 ,  0.00 ));
        assert_aq!(output[7], c128( 0.99 ,  2.41 ));
    }

    #[test]
    fn test_idft() {
        let     input  = [c128(4.00,  0.00), c128(1.00, -2.41),
                          c128(0.00,  0.00), c128(1.00, -0.41),
                          c128(0.00,  0.00), c128(0.99,  0.41),
                          c128(0.00,  0.00), c128(0.99,  2.41)];
        let mut output = [c128(0.0, 0.0); 8];
        idft(&input, &mut output);
        assert_aq!(output[0], c128( 1.00 ,  0.00 ));
        assert_aq!(output[1], c128( 1.00 ,  0.00 ));
        assert_aq!(output[2], c128( 1.00 ,  0.00 ));
        assert_aq!(output[3], c128( 1.00 ,  0.00 ));
        assert_aq!(output[4], c128( 0.00 ,  0.00 ));
        assert_aq!(output[5], c128( 0.00 ,  0.00 ));
        assert_aq!(output[6], c128( 0.00 ,  0.00 ));
        assert_aq!(output[7], c128( 0.00 ,  0.00 ));
    }
}
