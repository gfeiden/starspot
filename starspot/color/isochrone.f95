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

module isochrone
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)

    real(dp), dimension(:,:), allocatable :: isochrone, mags

contains

    subroutine init(feh, afe, isochrone_length, numbers_of_mags, filters)
        ! initialize arrays for input isochrone and magnitudes and 
        ! generate bolometric correction and color tables at the given
        ! metallicity.

        ! allocate input and output arrays
        if (allocated(isochrone) .eqv. .false.) then
            allocate(isochrone(isochrone_length, 12))
        end if

        if (allocated(mags) .eqv. .false.) then
            allocate(mags(isochrone_length, 12))
        end if

        ! interpolate tables for the given composition
        call some_subroutine(feh, afe)
    end subroutine init

    subroutine add_color(feh, afe, isochrone, isochrone_length, mags)
        ! color in the isochrone, including contributions from spots.
    end subroutine add_color

    subroutine clean()
        ! free up memory from allocated arrays
        if (allocated(isochrone) .eqv. .true.) deallocate(isochrone)
        if (allocated(mags) .eqv. .true.) deallocate(mags)

    end subroutine

end module isochrone
