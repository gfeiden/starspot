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

    integer :: n_passbands, n_teffs, n_loggs, n_fehs
    integer, dimension(:), allocatable :: end_teff

    real(dp), dimension(:), allocatable :: fehs, teffs, loggs
    real(dp), dimension(:,:,:,:), allocatable :: bc_tables

    character(len=15) :: bc_type
    character(len=5), dimension(:), allocatable :: passbands

    public  :: bc_init, bc_eval, bc_clean
    private :: marcs, phoenix_amescond, vc_semiemp

contains

    subroutine bc_init(brand, filters)
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
            case ('marcs')
                call marcs()
            case ('phx_amescond')
                call phoenix_amescond()
            case ('vc_semiemp')
                call vc_semiemp()
            case default
                call log_warn('invalid bc_type in bc_init: default to marcs08')
                call marcs()
        end select

        call log_note('bolometric correction tables successfully initialized')
    end subroutine bc_init

    subroutine bc_eval(feh, afe)
    end subroutine bc_eval

    subroutine bc_clean()
        if (allocated(passbands) .eqv. .true.) deallocate(passbands)
        if (allocated(bc_tables) .eqv. .true.) deallocate(bc_tables)
        if (allocated(teffs) .eqv. .true.) deallocate(teffs)
        if (allocated(loggs) .eqv. .true.) deallocate(loggs)
        if (allocated(end_teff) .eqv. .true.) deallocate(end_teff)
    end subroutine bc_clean

    subroutine marcs()
        integer :: i, j, k, m  ! loop iterators
        integer :: ioerr       ! error flag

        character :: directory, filename

        real(dp) :: ebv

        if (afe < -10.0) then
            directory = 'tab/std/'
        else if (afe < -0.2) then
            directory = 'tab/m04/'
        else if (afe < 0.2) then
            directory = 'tab/p00/'
        else if (afe < 0.4) then
            directory = 'tab/p04/'
        else
            call log_warn('strange [a/Fe] requested. default to standard')
            directory = 'tab/std/'

        ! open data header, allocate memory, read header, and close
        open(unit=90, file=directory // 'header.data', status='old', iostat=ioerr)
        if (ioerr /= 0) then
            call log_error('header file IO error in marcs08, cannot continue')
            stop
        end if
        read(90, '(1x, i3, 14x, i2, 14x, i2, 14x, i2, 21x, f6.3)') n_teffs, n_loggs, n_fehs, ebv

        if (allocated(fehs) .eqv. .false.) allocate(fehs(n_fehs))
        if (allocated(teffs) .eqv. .false.) allocate(teffs(n_teffs))
        if (allocated(loggs) .eqv. .false.) allocate(loggs(n_loggs))
        if (allocated(end_teff) .eqv. .false.) allocate(end_teff(n_loggs))
        read(90, '(' // n_teffs // 'f6.0)') (teffs(i), i = 1, n_teffs)
        read(90, '(' // n_loggs // 'f4.1)') (loggs(i), i = 1, n_loggs)
        read(90, '(' // n_loggs // 'i4') (end_teff(i), i = 1, n_loggs)
        close(90)

        ! allocate memory for tabulated data
        if (allocated(bc_tables) .eqv. .false.) then
            allocate(bc_tables(n_passbands, n_fehs, n_loggs, n_teffs))
        end if

        ! pass each passband table to memory
        do i = 1, n_passbands

            ! format directory and filename
            select case(trim(passbands(i)))
                case('U', 'B', 'V', 'R', 'I', 'Rc', 'Ic')
                    directory = directory // 'ubvri12/'
                    filename = 'jc_' // trim(passbands(i)) // '.data'
                case('Ux', 'Bx', 'U_90', 'B_90', 'V_90', 'R_90', 'I_90', 'Rc_90', 'Ic_90')
                    directory = directory // 'ubvri90/'
                    filename = 'jc_' // trim(passbands(i)) // '.data'
                case('J', 'H', 'K')
                    directory = directory // '2mass/'
                    filename = '2mass_' // trim(passbands(i)) // '.data'
                case('u', 'g', 'r', 'i', 'z')
                    directory = directory // 'sdss/'
                    filename = 'sdss_' // trim(passbands(i)) // '.data'
                case default
                    call log_error('invalid passband requested')
                    stop
            end select

            ! read bolometric corrections into memory
            open(90, file=filename, status='old', iostat=ioerr)
            if (ioerr /= 0) then
                call log_error('bc table file IO error in marcs08, cannot continue')
                stop
            end if
            do j = 1, n_fehs
                read(90, '(11x, f5.2)') fehs(j)
                do k = 1, n_loggs
                    read(90, '(' // end_teff(k) // 'f8.4)') (bc_tables(i, j, k, m), m = 1, end_teff(k))
                end do ! logg loop
            end do ! [Fe/H] loop
            close(90)

        end do ! passband loop

    end subroutine marcs

    subroutine phoenix_amescond()
    end subroutine phoenix_amescond

    subroutine vc_semiemp()
    end subroutine vc_semiemp


end module bolcorrection