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

module bolcorrection
    implicit none
    use utils

    integer :: n_passbands

    real(dp), dimension(:,:,:), allocatable :: bc_tables

    character(len=15) :: bc_type
    character(len=5), dimension(:), allocatable :: passbands

    public  :: bc_init, bc_eval, bc_clean
    private :: marcs08, phoenix_amescond, vc03

contains

    subroutine bc_init(feh, afe, brand, filters)
        character(len=15), intent(in) :: brand
        character(len=5), dimension(:), intent(in) :: filters

        n_passbands = size(filters, 1)
        if (allocated(passbands) .eqv. .false.) allocate(passbands(n_passbands))

        bc_type = trim(brand)
        do i = 1, n_passbands
            passbands(i) = trim(filters(i))
        end do

        ! initialize bolometric correction tables
        select case (bc_type)
            case ('marcs08')
                call marcs08(feh, afe)
            case ('phx_amescond')
                call phoenix_amescond(feh, afe)
            case ('vc03')
                call vc03(feh, afe)
            case default
                call log_warn('invalid bc_type in bc_init: default to marcs08')
                call marcs08(feh, afe)
        end select

        call log_note('bolometric correction tables successfully initialized')
    end subroutine bc_init

    subroutine bc_eval()
    end subroutine bc_eval

    subroutine bc_clean()
        if (allocated(passbands) .eqv. .true.) deallocate(passbands)
        if (allocated(bc_table) .eqv. .true.) deallocate(bc_table)
    end subroutine bc_clean

    subroutine marcs08()
    end subroutine marcs08

    subroutine phoenix_amescond()
    end subroutine phoenix_amescond

    subroutine vc03()
    end subroutine vc03


end module bolcorrection