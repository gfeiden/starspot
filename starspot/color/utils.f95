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

module utils
    implicit none

    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: qp = selected_real_kind(33, 4931)
    integer, parameter :: log_unit = 99

    !private :: log_error, log_note, log_warn
    !public  :: log_init, log_close

contains
    subroutine log_init(filename)
        character(len=*) :: filename
        open(log_unit, file=trim(filename), status='unknown')
    end subroutine log_init

    subroutine log_warn(message)
        character(len=*) :: message
        write(log_unit, '(" ---- WARNING: ", A)') message
    end subroutine log_warn

    subroutine log_error(message)
        character(len=*) :: message
        write(log_unit, '(" **** ERROR: ", A)') message
    end subroutine log_error

    subroutine log_note(message)
        character(len=*) :: message
        write(log_unit, '(" NOTE: ", A)') message
    end subroutine log_note

    subroutine log_close()
        close(log_unit)
    end subroutine log_close

end module utils
