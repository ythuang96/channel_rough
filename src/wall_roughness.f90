module wall_roughness
    use types, only: sp
    use mpi
    implicit none
    private
    include 'ctes3D'
    real(kind=sp), parameter :: disturbanceU = 0.05_sp
    real(kind=sp), parameter :: disturbanceV = -0.05_sp
    real(kind=sp), parameter :: targetWavenumberX = 0.0_sp
    real(kind=sp), parameter :: absoluteValueTargetWavenumberZ = 3.0_sp
    integer :: idxWavenumberX, idxPosWavenumberZ, idxNegWavenumberZ
    public :: initialize_wall_roughness, set_wall_roughness

contains
    subroutine initialize_wall_roughness(complexExponentVectorX, complexExponentVectorZ)
        complex(kind=sp), intent(in) :: complexExponentVectorX(0:), complexExponentVectorZ(0:)
        idxWavenumberX = find_index_target_wavenumber(targetWavenumberX, complexExponentVectorX)
        idxPosWavenumberZ = find_index_target_wavenumber(absoluteValueTargetWavenumberZ, complexExponentVectorZ)
        idxNegWavenumberZ = find_index_target_wavenumber(-absoluteValueTargetWavenumberZ, complexExponentVectorZ)
        call verify_wavenumber_index(idxWavenumberX, complexExponentVectorX, targetWavenumberX)
        call verify_wavenumber_index(idxPosWavenumberZ, complexExponentVectorZ, absoluteValueTargetWavenumberZ)
        call verify_wavenumber_index(idxNegWavenumberZ, complexExponentVectorZ, -absoluteValueTargetWavenumberZ)
        call print_information_to_console()
    end subroutine

    subroutine set_wall_roughness(uWallBottom, uWallTop, vWallBottom, vWallTop)
        complex(kind=sp), intent(out) :: uWallBottom(0:,0:), uWallTop(0:,0:), vWallBottom(0:,0:), vWallTop(0:,0:)
        uWallBottom = (0.0_sp, 0.0_sp)
        uWallTop = (0.0_sp, 0.0_sp)
        vWallBottom = (0.0_sp, 0.0_sp)
        vWallTop = (0.0_sp, 0.0_sp)
        ! perturb u
        uWallBottom(idxWavenumberX, idxPosWavenumberZ) = (disturbanceU, 0.0_sp)
        uWallBottom(idxWavenumberX, idxNegWavenumberZ) = (disturbanceU, 0.0_sp)
        uWallTop(idxWavenumberX, idxPosWavenumberZ) = (disturbanceU, 0.0_sp)
        uWallTop(idxWavenumberX, idxNegWavenumberZ) = (disturbanceU, 0.0_sp)
        ! perturb v
        vWallBottom(idxWavenumberX, idxPosWavenumberZ) = (disturbanceV, 0.0_sp)
        vWallBottom(idxWavenumberX, idxNegWavenumberZ) = (disturbanceV, 0.0_sp)
        ! invert sign at upper wall to have a symmetric configuration
        vWallTop(idxWavenumberX, idxPosWavenumberZ) = -(disturbanceV, 0.0_sp)
        vWallTop(idxWavenumberX, idxNegWavenumberZ) = -(disturbanceV, 0.0_sp)
    end subroutine

    pure function find_index_target_wavenumber(targetWavenumber, complexExponentVector) result(idxTargetWavenumber)
        real(kind=sp), intent(in) :: targetWavenumber
        complex(kind=sp), intent(in) :: complexExponentVector(0:)
        integer :: idxTargetWavenumber
        real(kind=sp), dimension(size(complexExponentVector)) :: deltaWavenumberVector
        integer, parameter :: offsetArrayStart = 1
        deltaWavenumberVector = aimag(complexExponentVector) - targetWavenumber
        ! idxTargetWavenumber has to be in (0, 1, ..., mx1) or (0, 1, ..., mz1)
        ! minloc returns values in (1, ... mx1+1) and (1, ..., mz1+1) -> subtract offsetArrayStart
        idxTargetWavenumber = minloc(abs(deltaWavenumberVector), dim = 1) - offsetArrayStart  ! dim required for scalar output
    end function

    subroutine verify_wavenumber_index(idx, complexExponentVector, targetWavenumber)
        integer, intent(in) :: idx
        complex(kind=sp), intent(in) :: complexExponentVector(0:)
        real(kind=sp), intent(in) :: targetWavenumber
        real(kind=sp) :: effectiveWavenumber
        real(kind=sp), parameter :: tolerance = 1e-4
        effectiveWavenumber = aimag(complexExponentVector(idx))
        if(abs(targetWavenumber - effectiveWavenumber) > tolerance) then
            stop 'Error in wall_roughness. Target wavenumber not found.'
        endif
    end subroutine

    subroutine print_information_to_console()
        integer :: mpiRank, mpiError
        integer, parameter :: mpiMaster = 0
        call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)
        if(mpiRank == mpiMaster) then
            write(*,*) 'Initialize wall roughness. Parameters'
            write(*,'(a16, f6.2)') 'x wavenumber = ', targetWavenumberX
            write(*,'(a22, f6.2)') 'z wavenumber (pos) = ', absoluteValueTargetWavenumberZ
            write(*,'(a22, f6.2)') 'z wavenumber (neg) = ', -absoluteValueTargetWavenumberZ
            write(*,'(a12, f6.2)') '|u_wall| = ', disturbanceU
            write(*,'(a12, f6.2)') '|v_wall| = ', disturbanceV
        endif
    end subroutine

end module

