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
    use utils
    implicit none

    integer :: n_passbands, n_teffs, n_loggs, n_fehs
    integer, dimension(:), allocatable :: end_teff

    real(dp), dimension(:), allocatable :: fehs, teffs, loggs
    real(dp), dimension(:,:,:), allocatable :: bc_table

    character(len=15) :: bc_type
    character(len=5), dimension(:), allocatable :: passbands

    ! no private routines allowed with f2py
    !private :: marcs, phoenix_amescond, vc_semiemp
    !public  :: bc_init, bc_eval, bc_clean

contains

    subroutine bc_init(feh, afe, brand, filters)
        integer :: i

        real(dp), intent(in) :: feh, afe

        character(len=15), intent(in) :: brand
        character(len=1), dimension(:), intent(in) :: filters

        !f2py intent(in) feh, afe, brand, filters

        call log_note('*** starting run ***')
        n_passbands = size(filters, 1)
        if (allocated(passbands) .eqv. .false.) allocate(passbands(n_passbands))

        bc_type = trim(brand)
        do i = 1, n_passbands
            passbands(i) = trim(filters(i))
        end do

        call log_note('initializing bolometric correction tables')
        ! initialize bolometric correction tables
        select case (trim(bc_type))
            case ('marcs')
                call marcs(feh, afe)
            case ('phx_amescond')
                call phoenix_amescond()
            case ('vc_semiemp')
                call vc_semiemp()
            case ('mamjek')
                call log_warn('Pecaut and Mamajek only valid for MS, solar metallicity')
                call bc_mamjek(0.0, 0.0, 1)
            case default
                call log_warn('invalid bc_type in bc_init: default to marcs08')
                call marcs(feh, afe)
        end select

        call log_note('bolometric correction tables successfully initialized')
    end subroutine bc_init

    subroutine bc_eval(teff, logg, logl, mag_length, magnitudes)
        use interpolate, only: lagrange

        integer :: i, j, k, teff_index, logg_index, mag_length
        real(dp) :: m_bol
        real(dp), intent(in) :: teff, logg, logl
        real(dp), parameter  :: m_bol_sun = 4.74
        real(dp), dimension(4) :: coeffs
        real(dp), dimension(:,:), allocatable :: bc_loggs
        real(dp), dimension(mag_length), intent(out) :: magnitudes

        !f2py intent(in) teff, logg, logl, mag_length
        !f2py intent(out) magnitudes
        !f2py depend(mag_length) magnitudes

        ! allocatable output arrays don't play well with f2py
        !
        !if (allocated(magnitudes) .eqv. .false.) then
        !    call log_note('allocating memory for magnitudes in bc_eval')
        !    allocate(magnitudes(n_passbands))
        !else
        !    call log_warn('magnitudes in bc_eval already allocated')
        !    deallocate(magnitudes)
        !    allocate(magnitudes(n_passbands))
        !end if

        ! define bolometric magnitude
        m_bol = m_bol_sun - 2.5*logl

        if (trim(brand) == 'mamjek') then
            call bc_mamjek(teff, logl, mag_length, 0, magnitudes)
            return
        end if

        ! hunt for logg and teff indicies
        do i = 1, n_loggs
            if (loggs(i) > logg) then
                j = i
                exit
            else
                cycle
            end if
        end do

        do i = 1, n_teffs
            if (teffs(i) > teff) then
                k = i
                exit
            else
                cycle
            end if
        end do

        ! protect against edge effects
        if (n_loggs - j < 1) then
            logg_index = n_loggs - 1
        else if (j < 3) then
            logg_index = 3
        else
            logg_index = j
        end if

        if (n_teffs - k < 1) then
            teff_index = n_teffs - 1
        else if (k < 3) then
            teff_index = 3
        else
            teff_index = k
        end if

        ! get interpolation coefficients for teff
        call lagrange(teffs(teff_index - 2:teff_index + 1), coeffs, teff, 4)

        ! intepolate in teff
        if (allocated(bc_loggs) .eqv. .false.) allocate(bc_loggs(n_passbands, n_loggs))
        do i = 1, n_passbands
            do j = 1, n_loggs
                bc_loggs(i, j) = coeffs(1)*bc_table(i, j, teff_index - 2) + &
                                 coeffs(2)*bc_table(i, j, teff_index - 1) + &
                                 coeffs(3)*bc_table(i, j, teff_index    ) + &
                                 coeffs(4)*bc_table(i, j, teff_index + 1)
            end do
        end do

        ! get interpolation coefficients for logg
        call lagrange(loggs(logg_index - 2:logg_index + 1), coeffs, logg, 4)

        ! interpolate at requested logg
        do i = 1, n_passbands
            magnitudes(i) = coeffs(1)*bc_loggs(i, logg_index - 2) + &
                            coeffs(2)*bc_loggs(i, logg_index - 1) + &
                            coeffs(3)*bc_loggs(i, logg_index    ) + &
                            coeffs(4)*bc_loggs(i, logg_index + 1)

            magnitudes(i) = m_bol - magnitudes(i)
        end do

        if (allocated(bc_loggs) .eqv. .true.) deallocate(bc_loggs)

    end subroutine bc_eval

    subroutine bc_clean()
        if (allocated(passbands) .eqv. .true.) deallocate(passbands)
        if (allocated(bc_table) .eqv. .true.) deallocate(bc_table)
        if (allocated(teffs) .eqv. .true.) deallocate(teffs)
        if (allocated(loggs) .eqv. .true.) deallocate(loggs)
        if (allocated(end_teff) .eqv. .true.) deallocate(end_teff)
    end subroutine bc_clean

    subroutine marcs(feh, afe)
        use interpolate, only: lagrange

        integer :: i, j, k, m  ! loop iterators
        integer :: ioerr       ! error flag
        integer :: feh_index, bcv

        character(len=132) :: directory, filename

        real(dp) :: ebv
        real(dp), intent(in) :: feh, afe
        real(dp), dimension(4) :: coeffs
        real(dp), dimension(:,:,:,:), allocatable :: bc_tables

        if (afe < -10.0) then
            directory = 'color/tab/std/'
        else if (afe < -0.2) then
            directory = 'color/tab/m04/'
        else if (afe < 0.2) then
            directory = 'color/tab/p00/'
        else if (afe < 0.4) then
            directory = 'color/tab/p04/'
        else
            call log_warn('strange [a/Fe] requested. default to standard')
            directory = 'tab/std/'
        end if

        ! open data header, allocate memory, read header, and close
        open(unit=90, file=trim(directory) // 'header.data', status='old', iostat=ioerr)
        if (ioerr /= 0) then
            call log_error('header file IO error in marcs08, cannot continue')
            stop
        end if
        read(90, '(i3, 14x, i2, 14x, i2, 14x, i2, 21x, f6.3)') n_teffs, n_loggs, n_fehs, bcv, ebv

        if (allocated(fehs) .eqv. .false.) allocate(fehs(n_fehs))
        if (allocated(teffs) .eqv. .false.) allocate(teffs(n_teffs))
        if (allocated(loggs) .eqv. .false.) allocate(loggs(n_loggs))
        if (allocated(end_teff) .eqv. .false.) allocate(end_teff(n_loggs))
        read(90, *) (teffs(i), i = 1, n_teffs)
        read(90, *) (loggs(i), i = 1, n_loggs)
        read(90, *) (end_teff(i), i = 1, n_loggs)
        !read(90, '(' // n_teffs // 'f6.0)') (teffs(i), i = 1, n_teffs)
        !read(90, '(' // n_loggs // 'f4.1)') (loggs(i), i = 1, n_loggs)
        !read(90, '(' // n_loggs // 'i4') (end_teff(i), i = 1, n_loggs)
        close(90)

        ! allocate memory for tabulated data
        if (allocated(bc_tables) .eqv. .false.) then
            call log_note('allocating memory for bc tables')
            allocate(bc_tables(n_passbands, n_fehs, n_loggs, n_teffs))
        else
            call log_warn('bc tables already allocated to memory, reallocating')
            deallocate(bc_tables)
            allocate(bc_tables(n_passbands, n_fehs, n_loggs, n_teffs))
        end if

        ! pass each passband table to memory
        do i = 1, n_passbands

            ! format directory and filename
            select case(trim(passbands(i)))
                case('U', 'B', 'V', 'R', 'I', 'Rc', 'Ic')
                    filename = trim(directory) // 'ubvri12/' // 'jc_' // trim(passbands(i)) // '.data'
                case('Ux', 'Bx', 'U_90', 'B_90', 'V_90', 'R_90', 'I_90', 'Rc_90', 'Ic_90')
                    filename = trim(directory) // 'ubvri90/' // 'jc_' // trim(passbands(i)) // '.data'
                case('J', 'H', 'K')
                    filename = trim(directory) // '2mass/' // '2mass_' // trim(passbands(i)) // '.data'
                case('u', 'g', 'r', 'i', 'z')
                    filename = trim(directory) // 'sdss/' // 'sdss_' // trim(passbands(i)) // '.data'
                case default
                    call log_error('invalid passband requested')
                    stop
            end select

            ! read bolometric corrections into memory
            open(90, file=trim(filename), status='old', iostat=ioerr)
            if (ioerr /= 0) then
                call log_error('bc table file IO error in marcs, cannot continue: ' // trim(filename))
                stop
            end if
            call log_note('reading bc table ' // trim(filename) // ' into memory')
            do j = 1, n_fehs
                read(90, '(10x, f6.2)') fehs(j)
                do k = 1, n_loggs
                    read(90, *) (bc_tables(i, j, k, m), m = 1, end_teff(k))
                    !read(90, '(' // end_teff(k) // 'f8.4)') (bc_tables(i, j, k, m), m = 1, end_teff(k))
                end do ! logg loop
            end do ! [Fe/H] loop
            close(90)

        end do ! passband loop

        if (allocated(bc_table) .eqv. .false.) then
            call log_note('allocating memory for bolometric correction table')
            allocate(bc_table(n_passbands, n_loggs, n_teffs))
        else
            call log_warn('bolometric correction table already defined: redefining')
            deallocate(bc_table)
            allocate(bc_table(n_passbands, n_loggs, n_teffs))
        end if

        ! find indices of nearest 4 [Fe/H] values
        do i = 1, n_fehs
            if (fehs(i) > feh) then
                j = i
                exit
            else
                cycle
            end if
        end do

        ! protect against edge errors
        if (n_fehs - j < 1) then
            feh_index = n_fehs - 1
        else if (j < 3) then
            feh_index = 3
        else
            feh_index = j
        end if

        ! get interpolation coefficients with 4th order lagrange method
        call lagrange(fehs(feh_index - 2:feh_index + 1), coeffs, feh, 4)

        ! interpolate at the requested [Fe/H]
        do i = 1, n_passbands
            do j = 1, n_loggs
                do k = 1, end_teff(j)
                    bc_table(i, j, k) = coeffs(1)*bc_tables(i, feh_index - 2, j, k) + &
                                        coeffs(2)*bc_tables(i, feh_index - 1, j, k) + &
                                        coeffs(3)*bc_tables(i, feh_index    , j, k) + &
                                        coeffs(4)*bc_tables(i, feh_index + 1, j, k)
                end do
            end do
        end do

        if (allocated(bc_tables) .eqv. .true.) deallocate(bc_tables)
    end subroutine marcs

    subroutine phoenix_amescond()
        return
    end subroutine phoenix_amescond

    subroutine vc_semiemp()
        return
    end subroutine vc_semiemp

    subroutine bc_mamjek(teff, logl, mag_length, init, magnitudes)
        use interpolate, only: lagrange

        integer,  intent(in) :: init, mag_length
        real(dp), intent(in) :: teff, logl
        real(dp), dimension(mag_length), intent(out) :: magnitudes

        if (init == 1) then
            ! read in the bolometric correction table
            n_teffs = 102
            n_loggs = 1
            n_fehs  = 0

            open(90, file="color/tab/pm13/")
            ! read file header

        end if

    end subroutine bc_mamjek


end module bolcorrection
