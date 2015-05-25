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
    use utils
    implicit none

    public  :: lagrange

contains

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
        else
            continue
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
