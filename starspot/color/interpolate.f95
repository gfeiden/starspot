! The MIT License (MIT)
!
! Copyright (c) 2015 Gregory Feiden
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module interpolate
    implicit none
    use utils

    public  :: akima_setup, akima_spline
    private :: akima_spline_edge

contains

    subroutine akima_setup()
    end subroutine akima_setup

    subroutine akima_spline_edge(x1, x2, x3, f1, f2, f3, f00, f01)
        real(dp) :: df21, df31, dx21, dx31, dx32, den, g1, g2
        real(dp), intent(in)  :: x1, x2, x2, f1, f2, f3
        real(dp), intent(out) :: f00, f01

        ! compute differences
        df21 = f2 - f1
        df31 = f3 - f1
        dx21 = x2 - x1
        dx31 = x3 - x1
        dx32 = x3 - x2

        ! compute coefficients
        den = dx21*dx32*dx31
        g1 = (df21*dx31*dx31 - df31*dx21*dx21)/den
        g2 = (df31*dx21 - df21*dx31)/den

        ! compute extended end points
        f00 = f1 - g1*dx32 + g2*dx32**2
        f01 = f1 - g1*dx31 + g2*dx31**2
    end subroutine akima_spline_edge

    subroutine akima_spline()
    end subroutine akima_spline

    subroutine lagrange(x, coeffs, x_new, degree)
        ! N-point Lagrangian interpolation
        integer  :: i, j, degree
        real(dp) :: x_new, denom
        real(dp), dimension(degree), intent(in)  :: x
        real(dp), dimension(degree), intent(out) :: coeffs

        if (degree <= 0) then
            call log_error('invalid degree for Lagrange interpolation, cannot continue')
            stop
        else if (degree < 3) then
            call log_warn('Lagrange interpolation attempted with degree < 3, beware of results')
        end if

        do i = 1, degree
            denom = 1.0_dp
            coeffs(i) = 1.0_dp
            do j = 1, degree
                if (j /= i) then
                    denom = denom*(x(i) - x(j))
                    coeffs(i) = coeffs(i)*(x_new - x(j))
                end if
            end do
            coeffs(i) = coeffs(i)/denom
        end do

    end subroutine lagrange


end module interpolate
