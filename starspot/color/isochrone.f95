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
    use utils
    implicit none

    real(dp), dimension(:,:), allocatable :: isochrone, mags

contains

    subroutine init(feh, afe, filters, isochrone_length, isochrone_width)
        ! initialize arrays for input isochrone and magnitudes and 
        ! generate bolometric correction and color tables at the given
        ! metallicity.
        use bolcorrection, only: bc_init

        integer :: isochrone_length, isochrone_width, numbers_of_mags

        ! allocate input and output arrays
        if (allocated(isochrone) .eqv. .false.) then
            allocate(isochrone(isochrone_length, isochrone_width))
        end if

        if (allocated(mags) .eqv. .false.) then
            allocate(mags(isochrone_length, numbers_of_mags))
        end if

        ! interpolate tables for the given composition
        call bc_init(feh, afe, brand, filters)

    end subroutine init

    subroutine add_color(feh, isochrone)
        ! color in the isochrone, including contributions from spots.
        real(dp), intent(in) :: feh
        real(dp) :: logg, logt_phot, logt_spot, logl_phot, logl_spot, surface_coverge

        do i = 1, size(isochrone, 1)
            ! define parameters for star
            if (surface_coverge > 0.0) then
                logg = isochrone(i, 2)
                logt_phot = isochrone(i, 4)
                logt_spot = isochrone(i, 5)
                logl_phot = isochrone(i, 6)
                logl_spot = isochrone(i, 7)
            else
                logg = isochrone(i, 3)
                logl_phot = isochrone(i, 4)
                logt_phot = isochrone(i, 2)
                logt_spot = logt_phot
                logl_spot = 0.0
            end if

            ! unspotted region
            call bc_eval(10.0**logt_phot, logg, logl_phot, mags)

            ! spotted region
            if (surface_coverge > 0.0)
                call bc_eval(10.0**logt_spot, logg, logl_spot, mags)
            end if
        end do
    end subroutine add_color

    subroutine clean()
        ! free up memory from allocated arrays
        if (allocated(isochrone) .eqv. .true.) deallocate(isochrone)
        if (allocated(mags) .eqv. .true.) deallocate(mags)

    end subroutine

end module isochrone
