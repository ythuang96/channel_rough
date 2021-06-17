module boundary_planes
    ! provide variables and functions for velocity boundary conditions
    use types, only: sp
    implicit none
    private
    include 'ctes3D'
    complex(kind=sp), dimension(0:mx1, 0:mz1) :: &
        uWallBottom, vWallBottom, wWallBottom,   &
        uWallTop, vWallTop, wWallTop
    ! boundary planes must be visible to other program units
    ! (Jimenez' code design)
    public :: uWallBottom, vWallBottom, wWallBottom, &
        uWallTop, vWallTop, wWallTop, set_zero_bc, &
        set_bc_from_restart_file
contains
    subroutine set_zero_bc()
        ! initialize all boundary planes to zero (this corresponds
        ! to no-slip and no-through boundary conditions)
        uWallBottom = (0.0_sp, 0.0_sp)
        uWallTop = (0.0_sp, 0.0_sp)
        vWallBottom = (0.0_sp, 0.0_sp)
        vWallTop = (0.0_sp, 0.0_sp)
        wWallBottom = (0.0_sp, 0.0_sp)
        wWallTop = (0.0_sp, 0.0_sp)
    end subroutine

    subroutine set_bc_from_restart_file(uWallBottomRestart, uWallTopRestart, &
                                        vWallBottomRestart, vWallTopRestart, &
                                        wWallBottomRestart, wWallTopRestart)
        ! set velocity boundary conditions to values read from restart file
        real(kind=sp), intent(in) :: &
            uWallBottomRestart(:,:), uWallTopRestart(:,:), &
            vWallBottomRestart(:,:), vWallTopRestart(:,:), &
            wWallBottomRestart(:,:), wWallTopRestart(:,:)
        call broadcast_bc_restart_file(uWallBottomRestart, uWallTopRestart, &
                                       vWallBottomRestart, vWallTopRestart, &
                                       wWallBottomRestart, wWallTopRestart)
        call convert_real_to_complex_bc(uWallBottomRestart, uWallTopRestart, &
                                        vWallBottomRestart, vWallTopRestart, &
                                        wWallBottomRestart, wWallTopRestart)
    end subroutine

    subroutine broadcast_bc_restart_file(uWallBottomRestart, uWallTopRestart, &
                                         vWallBottomRestart, vWallTopRestart, &
                                         wWallBottomRestart, wWallTopRestart)
        ! MPI master broadcasts data read from restart file
        include "mpif.h"
        real(kind=sp), intent(in) :: &
            uWallBottomRestart(:,:), uWallTopRestart(:,:), &
            vWallBottomRestart(:,:), vWallTopRestart(:,:), &
            wWallBottomRestart(:,:), wWallTopRestart(:,:)
        integer, parameter :: idMaster = 0  ! MPI master has id zero
        integer errorStatus
        call MPI_Bcast(uWallBottomRestart, size(uWallBottomRestart), &
            MPI_REAL, idMaster, MPI_COMM_WORLD, errorStatus)
        call MPI_Bcast(uWallTopRestart, size(uWallTopRestart), &
            MPI_REAL, idMaster, MPI_COMM_WORLD, errorStatus)
        call MPI_Bcast(vWallBottomRestart, size(vWallBottomRestart), &
            MPI_REAL, idMaster, MPI_COMM_WORLD, errorStatus)
        call MPI_Bcast(vWallTopRestart, size(vWallTopRestart), &
            MPI_REAL, idMaster, MPI_COMM_WORLD, errorStatus)
        call MPI_Bcast(wWallBottomRestart, size(wWallBottomRestart), &
            MPI_REAL, idMaster, MPI_COMM_WORLD, errorStatus)
        call MPI_Bcast(wWallTopRestart, size(wWallTopRestart), &
            MPI_REAL, idMaster, MPI_COMM_WORLD, errorStatus)
    end subroutine

    subroutine convert_real_to_complex_bc(uWallBottomRestart, uWallTopRestart, &
                                          vWallBottomRestart, vWallTopRestart, &
                                          wWallBottomRestart, wWallTopRestart)
        ! convert the real values read from the restart file to complex values
        ! used in the DNS
        real(kind=sp), intent(in) :: &
            uWallBottomRestart(:,:), uWallTopRestart(:,:), &
            vWallBottomRestart(:,:), vWallTopRestart(:,:), &
            wWallBottomRestart(:,:), wWallTopRestart(:,:)
        integer :: i, k, ii, kk
        integer :: mxRestart, mzRestart, mx1Restart, mz1Restart, &
            mx1Min, mz1Min
        ! initialize everything to zero (important if number of
        ! Fourier modes in current run is higher than number
        ! of Fourier modes in previous run)
        call set_zero_bc()
        ! obtain array size of restart file (all velocity boundary
        ! planes have the same size, we arbitrarily choose
        ! uWallBottomRestart here)
        mxRestart = size(uWallBottomRestart, 1)
        mzRestart = size(uWallBottomRestart, 2)
        mx1Restart = mxRestart / 2 - 1  ! identical to ctes3D
        mz1Restart = mzRestart - 1  ! identical to ctes3D
        ! only assign the smaller of the two number of modes
        mx1Min = min(mx1, mx1Restart)
        mz1Min = min(mz1, mz1Restart)
        ! real --> complex values and take care of size
        do k = 0, mz1Min
            kk = k + 1  ! recall: mz = mz1 + 1
            do i = 0, mx1Min
                ii = 2 * i + 2  ! recall: mx = 2*mx1 + 2
                uWallBottom(i,k) = &
                    cmplx(uWallBottomRestart(ii-1,kk), uWallBottomRestart(ii,kk))
                uWallTop(i,k) = &
                    cmplx(uWallTopRestart(ii-1,kk), uWallTopRestart(ii,kk))
                vWallBottom(i,k) = &
                    cmplx(vWallBottomRestart(ii-1,kk), vWallBottomRestart(ii,kk))
                vWallTop(i,k) = &
                    cmplx(vWallTopRestart(ii-1,kk), vWallTopRestart(ii,kk))
                wWallBottom(i,k) = &
                    cmplx(wWallBottomRestart(ii-1,kk), wWallBottomRestart(ii,kk))
                wWallTop(i,k) = &
                    cmplx(wWallTopRestart(ii-1,kk), wWallTopRestart(ii,kk))
            enddo
        enddo
    end subroutine

end module
