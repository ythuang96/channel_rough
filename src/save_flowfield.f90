module save_flowfield
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'ctes3D'
    ! Buffer for Fourier space velocity and vorticity
    complex(kind=sp) :: bufferStreamwiseVelocity(0:mx1,0:mz1,myp), bufferWallnormalVelocity(0:mx1,0:mz1,myp), &
        bufferSpanwiseVelocity(0:mx1,0:mz1,myp)
    complex(kind=sp) :: bufferStreamwiseVorticity(0:mx1,0:mz1,myp), bufferWallnormalVorticity(0:mx1,0:mz1,myp), &
        bufferSpanwiseVorticity(0:mx1,0:mz1,myp)
    ! Buffer for Fourier space phi = nabla^2 v
    ! complex(kind=sp) :: bufferPhi(0:mx1,0:mz1,myp)

    ! Buffer for Fourier space dvdy and do2dy
    ! complex(kind=sp) :: bufferdvdy(0:mx1,0:mz1,myp), bufferdo2dy(0:mx1,0:mz1,myp)
    ! Buffer for Fourier space H1, H2, H3
    ! complex(kind=sp) :: bufferH1(0:mx1,0:mz1,myp), bufferH2(0:mx1,0:mz1,myp), bufferH3(0:mx1,0:mz1,myp)

    ! Buffer for Fourier space non-linear forcing
    ! complex(kind=sp) :: bufferForcingVelocity(0:mx1,0:mz1,myp), bufferForcingVorticity(0:mx1,0:mz1,myp)

    ! Buffer for Fourier space wall acceleration
    ! complex(kind=sp) :: xAccelerationBottomWall(0:mx1,0:mz1), xAccelerationTopWall(0:mx1,0:mz1), &
    !     yAccelerationBottomWall(0:mx1,0:mz1), yAccelerationTopWall(0:mx1,0:mz1), &
    !     zAccelerationBottomWall(0:mx1,0:mz1), zAccelerationTopWall(0:mx1,0:mz1)

    ! Buffer for physical space variables
    ! real(kind=sp) :: bufferUPhys(mgalx+2,mgalz,myp), bufferVPhys(mgalx+2,mgalz,myp), bufferWPhys(mgalx+2,mgalz,myp)
    ! real(kind=sp) :: buffero1Phys(mgalx+2,mgalz,myp), buffero2Phys(mgalx+2,mgalz,myp), buffero3Phys(mgalx+2,mgalz,myp)

    ! Other parameters
    logical :: collectFlowfield
    ! logical :: collectWallVelocity
    integer :: samplingFrequencyInTimesteps, fileNumber
    real(kind=sp) :: samplingPeriod, timeLastSample
    character(:), allocatable :: baseFileNameWithPath

    ! Public Variables
    public :: collectFlowfield
    ! public :: collectWallVelocity

    ! Public Subroutines
    ! Initialize
    public :: initialize_save_flowfield
    ! Determine whether to collect flowfield or not
    public :: assess_whether_to_collect_flowfield, assess_whether_to_collect_flowfield_dt
    ! public :: assess_whether_to_collect_wall_velocity
    ! Save Fourier data at a single plane to buffer
    public :: save_velocity_at_wallparallel_plane_to_buffer, save_vorticity_at_wallparallel_plane_to_buffer
    ! public :: save_phi_at_wallparallel_plane_to_buffer, save_yderivatives_at_wallparallel_plane_to_buffer, &
    !     save_H_at_wallparallel_plane_to_buffer
    ! Save Physical Data at a single plane to buffer
    ! public :: save_velocity_at_wallparallel_plane_to_buffer_phys, save_vorticity_at_wallparallel_plane_to_buffer_phys
    ! Save Fourier Forcing Data at a single plane to buffer
    ! public :: save_velocity_forcing_to_buffer, save_vorticity_forcing_at_plane_to_buffer
    ! Save top and bottom wall velocity to buffer for computing wall acceleration
    ! public :: save_velocity_bottom_wall_to_buffer, save_velocity_top_wall_to_buffer
    ! Compute wall acceleration at top and bottom wall
    ! public :: compute_acceleration_bottom_wall, compute_acceleration_top_wall
    ! Write Data to HDF file
    public :: write_flowfield_to_hdf_file

contains
    ! ***********************************************************************************************************************************
    ! Initialize
    subroutine initialize_save_flowfield(runNameWithPath)
        character(len=*), intent(in) :: runNameWithPath
        integer :: mpiRank, mpiError
        integer, parameter :: mpiMaster = 0

        ! Set Sampling parameters
        samplingFrequencyInTimesteps = 1
        !samplingPeriod = 0.0899_sp  ! corresponds to approx tUc/h = 0.1 @ ReTau=660
        !samplingPeriod = 0.09782_sp  ! corresponds to approx tUc/h = 0.1 @ Smooth wall Retau=550
        samplingPeriod = 0.09475_sp  ! corresponds to approx tUc/h = 0.1 @ kz = 3 Retau=550
        collectFlowfield = .false.
        ! collectWallVelocity = .false.
        fileNumber = 1
        timeLastSample = 0.0_sp
        baseFileNameWithPath = trim(adjustl(runNameWithPath)) // '_flowfield.'

        ! Initialize Fourier Buffers
        bufferStreamwiseVelocity = (0.0_sp, 0.0_sp)
        bufferWallnormalVelocity = (0.0_sp, 0.0_sp)
        bufferSpanwiseVelocity = (0.0_sp, 0.0_sp)
        bufferStreamwiseVorticity = (0.0_sp, 0.0_sp)
        bufferWallnormalVorticity = (0.0_sp, 0.0_sp)
        bufferSpanwiseVorticity = (0.0_sp, 0.0_sp)
        ! bufferPhi = (0.0_sp, 0.0_sp)
        ! bufferForcingVelocity = (0.0_sp, 0.0_sp)
        ! bufferForcingVorticity = (0.0_sp, 0.0_sp)
        ! bufferdvdy = (0.0_sp, 0.0_sp)
        ! bufferdo2dy = (0.0_sp, 0.0_sp)
        ! bufferH1 = (0.0_sp, 0.0_sp)
        ! bufferH2 = (0.0_sp, 0.0_sp)
        ! bufferH3 = (0.0_sp, 0.0_sp)

        ! Initialize Physical Buffers
        ! bufferUPhys = 0.0_sp;
        ! bufferVPhys = 0.0_sp;
        ! bufferWPhys = 0.0_sp;
        ! buffero1Phys = 0.0_sp;
        ! buffero2Phys = 0.0_sp;
        ! buffero3Phys = 0.0_sp;

        ! Initialize Acceleration Buffers
        ! xAccelerationBottomWall = (0.0_sp, 0.0_sp)
        ! xAccelerationTopWall = (0.0_sp, 0.0_sp)
        ! yAccelerationBottomWall = (0.0_sp, 0.0_sp)
        ! yAccelerationTopWall = (0.0_sp, 0.0_sp)
        ! zAccelerationBottomWall = (0.0_sp, 0.0_sp)
        ! zAccelerationTopWall = (0.0_sp, 0.0_sp)
        
        ! Write sampling information
        call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)
        if (mpiRank == mpiMaster) then
            write(*,'(a45,i4)') 'sampling frequency flow field in timesteps: ', samplingFrequencyInTimesteps
            write(*,'(a18,f6.4)') 'sampling period: ', samplingPeriod
        endif
    end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Determine whether to collect flowfield or not
    ! Based on step number or based on time
    subroutine assess_whether_to_collect_flowfield(timeStepNr)
        integer, intent(in) :: timeStepNr
        if (mod(timeStepNr,samplingFrequencyInTimesteps) == 0) then
            collectFlowfield = .true.
        else
            collectFlowfield = .false.
        endif
    end subroutine

    subroutine assess_whether_to_collect_flowfield_dt(time)
        real(kind=sp), intent(in) :: time
        real(kind=sp) :: deltaT
        deltaT = time - timeLastSample
        if (deltaT >= samplingPeriod) then
            collectFlowfield = .true.
            timeLastSample = time
        else
            collectFlowfield = .false.
        endif
    end subroutine

    ! subroutine assess_whether_to_collect_wall_velocity(timeStepNr)
    !     integer, intent(in) :: timeStepNr
    !     ! to compute the acceleration, we need the velocities at the data collection time step and the preceding one
    !     if (mod(timeStepNr,samplingFrequencyInTimesteps) == 0 &  ! data collection time step
    !         .or. mod(timeStepNr+1,samplingFrequencyInTimesteps) == 0) then  ! preceding timestep
    !         collectWallVelocity = .true.
    !     else
    !         collectWallVelocity = .false.
    !     endif
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Save Fourier data at a single plane to buffer
    subroutine save_velocity_at_wallparallel_plane_to_buffer(streamwiseVelocity, wallnormalVelocity, &
        spanwiseVelocity, indexDataPlane, indexStartingPlaneProcessor)
        complex(kind=sp), intent(in) :: streamwiseVelocity(:,:), wallnormalVelocity(:,:), spanwiseVelocity(:,:)
        integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
        call save_data_wallparallel_plane_to_buffer(streamwiseVelocity, indexDataPlane, &
            indexStartingPlaneProcessor, bufferStreamwiseVelocity)
        call save_data_wallparallel_plane_to_buffer(wallnormalVelocity, indexDataPlane, &
            indexStartingPlaneProcessor, bufferWallnormalVelocity)
        call save_data_wallparallel_plane_to_buffer(spanwiseVelocity, indexDataPlane, &
            indexStartingPlaneProcessor, bufferSpanwiseVelocity)
    end subroutine

    subroutine save_vorticity_at_wallparallel_plane_to_buffer(streamwiseVorticity, wallnormalVorticity, &
        spanwiseVorticity, indexDataPlane, indexStartingPlaneProcessor)
        complex(kind=sp), intent(in) :: streamwiseVorticity(:,:), wallnormalVorticity(:,:), spanwiseVorticity(:,:)
        integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
        call save_data_wallparallel_plane_to_buffer(streamwiseVorticity, indexDataPlane, &
            indexStartingPlaneProcessor, bufferStreamwiseVorticity)
        call save_data_wallparallel_plane_to_buffer(wallnormalVorticity, indexDataPlane, &
            indexStartingPlaneProcessor, bufferWallnormalVorticity)
        call save_data_wallparallel_plane_to_buffer(spanwiseVorticity, indexDataPlane, &
            indexStartingPlaneProcessor, bufferSpanwiseVorticity)
    end subroutine

    ! subroutine save_phi_at_wallparallel_plane_to_buffer(Phi, indexDataPlane, indexStartingPlaneProcessor)
    !     complex(kind=sp), intent(in) :: Phi(0:mx1,0:mz1)
    !     integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
    !     call save_data_wallparallel_plane_to_buffer(Phi, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferPhi)
    ! end subroutine

    ! subroutine save_yderivatives_at_wallparallel_plane_to_buffer(dvdy, do2dy, &
    !     indexDataPlane, indexStartingPlaneProcessor)
    !     complex(kind=sp), intent(in) :: dvdy(0:mx1,0:mz1), do2dy(0:mx1,0:mz1)
    !     integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
    !     call save_data_wallparallel_plane_to_buffer(dvdy, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferdvdy)
    !     call save_data_wallparallel_plane_to_buffer(do2dy, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferdo2dy)
    ! end subroutine

    ! subroutine save_H_at_wallparallel_plane_to_buffer(H1, H2, &
    !     H3, indexDataPlane, indexStartingPlaneProcessor)
    !     complex(kind=sp), intent(in) :: H1(:,:), H2(:,:), H3(:,:)
    !     integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
    !     call save_data_wallparallel_plane_to_buffer(H1, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferH1)
    !     call save_data_wallparallel_plane_to_buffer(H2, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferH2)
    !     call save_data_wallparallel_plane_to_buffer(H3, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferH3)
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Save Physical Data at a single plane to buffer
    ! subroutine save_velocity_at_wallparallel_plane_to_buffer_phys(streamwiseVelocity, wallnormalVelocity, &
    !     spanwiseVelocity, indexDataPlane, indexStartingPlaneProcessor)
    !     real(kind=sp), intent(in) :: streamwiseVelocity(:,:), wallnormalVelocity(:,:), spanwiseVelocity(:,:)
    !     integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
    !     call save_data_wallparallel_plane_to_buffer_phys(streamwiseVelocity, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferUPhys)
    !     call save_data_wallparallel_plane_to_buffer_phys(wallnormalVelocity, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferVPhys)
    !     call save_data_wallparallel_plane_to_buffer_phys(spanwiseVelocity, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferWPhys)
    ! end subroutine

    ! subroutine save_vorticity_at_wallparallel_plane_to_buffer_phys(streamwiseVorticity, wallnormalVorticity, &
    !     spanwiseVorticity, indexDataPlane, indexStartingPlaneProcessor)
    !     real(kind=sp), intent(in) :: streamwiseVorticity(:,:), wallnormalVorticity(:,:), spanwiseVorticity(:,:)
    !     integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
    !     call save_data_wallparallel_plane_to_buffer_phys(streamwiseVorticity, indexDataPlane, &
    !         indexStartingPlaneProcessor, buffero1Phys)
    !     call save_data_wallparallel_plane_to_buffer_phys(wallnormalVorticity, indexDataPlane, &
    !         indexStartingPlaneProcessor, buffero2Phys)
    !     call save_data_wallparallel_plane_to_buffer_phys(spanwiseVorticity, indexDataPlane, &
    !         indexStartingPlaneProcessor, buffero3Phys)
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Save Fourier Forcing Data at a single plane to buffer
    ! subroutine save_vorticity_forcing_at_plane_to_buffer(forcingVorticity, indexDataPlane, indexStartingPlaneProcessor)
    !     complex(kind=sp), intent(in) :: forcingVorticity(:,:)
    !     integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
    !     call save_data_wallparallel_plane_to_buffer(forcingVorticity, indexDataPlane, &
    !         indexStartingPlaneProcessor, bufferForcingVorticity)
    ! end subroutine

    ! subroutine save_velocity_forcing_to_buffer(forcingVelocityYXZ, indexStartingPlaneProcessor, &
    !     indexEndingPlaneProcessor)
    !     real(kind=sp), intent(in) :: forcingVelocityYXZ(:,:,:)
    !     integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
    !     real(kind=sp) :: forcingVelocityXZYAsReal(mx, mz, myp)
    !     complex(kind=sp) :: forcingVelocityXZYAsComplex(0:mx1, 0:mz1, myp)
    !     real(kind=sp), allocatable :: workBuffer(:)
    !     integer :: workBufferSize, mpiRank, mpiError, i, j, indexDataPlane
    !     workBufferSize = mx * nxymax  ! mx and nxymax defined in ctes3D
    !     allocate(workBuffer(workBufferSize))
    !     call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)
    !     call chjik2ikj(forcingVelocityYXZ, forcingVelocityXZYAsReal, workBuffer, workBuffer, mpiRank)
    !     do i = 0, mx1
    !         forcingVelocityXZYAsComplex(i,:,:) = cmplx(forcingVelocityXZYAsReal(2*i+1,:,:), &
    !             forcingVelocityXZYAsReal(2*i+2,:,:))
    !     enddo
    !     do j = 1, indexEndingPlaneProcessor-indexStartingPlaneProcessor+1
    !         indexDataPlane = indexStartingPlaneProcessor + j - 1
    !         call save_data_wallparallel_plane_to_buffer(forcingVelocityXZYAsComplex(:,:,j), indexDataPlane, &
    !             indexStartingPlaneProcessor, bufferForcingVelocity)
    !     enddo
    !     deallocate(workBuffer)
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Save data at a single wall parallel plane to buffer
    ! A complex version and a real version
    subroutine save_data_wallparallel_plane_to_buffer(dataAtWallparallelPlane, indexDataPlane, &
        indexStartingPlaneProcessor, dataBuffer)
        complex(kind=sp), intent(in) :: dataAtWallparallelPlane(:,:)
        integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
        complex(kind=sp), intent(out) :: dataBuffer(:,:,:)
        integer :: indexDataBuffer
        ! each processor operates on a subset of wallparallel planes, starting at the plane with index indexStartingPlaneProcessor
        ! to store the data in the buffer, we have to translate the global indexDataPlane to a local indexDataBuffer
        indexDataBuffer = indexDataPlane - indexStartingPlaneProcessor + 1
        dataBuffer(:,:,indexDataBuffer) = dataAtWallparallelPlane
    end subroutine

    ! subroutine save_data_wallparallel_plane_to_buffer_phys(dataAtWallparallelPlane, indexDataPlane, &
    !     indexStartingPlaneProcessor, dataBuffer)
    !     real(kind=sp), intent(in) :: dataAtWallparallelPlane(:,:)
    !     integer, intent(in) :: indexDataPlane, indexStartingPlaneProcessor
    !     real(kind=sp), intent(out) :: dataBuffer(:,:,:)
    !     integer :: indexDataBuffer
    !     ! each processor operates on a subset of wallparallel planes, starting at the plane with index indexStartingPlaneProcessor
    !     ! to store the data in the buffer, we have to translate the global indexDataPlane to a local indexDataBuffer
    !     indexDataBuffer = indexDataPlane - indexStartingPlaneProcessor + 1
    !     dataBuffer(:,:,indexDataBuffer) = dataAtWallparallelPlane
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Save top and bottom wall velocity to buffer for computing wall acceleration
    ! subroutine save_velocity_bottom_wall_to_buffer(xVelocityBottomWall, yVelocityBottomWall, &
    !     zVelocityBottomWall)
    !     complex(kind=sp), intent(in) :: xVelocityBottomWall(:,:), yVelocityBottomWall(:,:), zVelocityBottomWall(:,:)
    !     ! to compute the time derivative, we first need to store the current value in the arrays
    !     ! the derivative is taken in compute_acceleration_bottom_wall below.
    !     xAccelerationBottomWall = xVelocityBottomWall
    !     yAccelerationBottomWall = yVelocityBottomWall
    !     zAccelerationBottomWall = zVelocityBottomWall
    ! end subroutine

    ! subroutine save_velocity_top_wall_to_buffer(xVelocityTopWall, yVelocityTopWall, zVelocityTopWall)
    !     complex(kind=sp), intent(in) :: xVelocityTopWall(:,:), yVelocityTopWall(:,:), zVelocityTopWall(:,:)
    !     ! to compute the time derivative, we first need to store the current value in the arrays
    !     ! the derivative is taken in compute_acceleration_top_wall below.
    !     xAccelerationTopWall = xVelocityTopWall
    !     yAccelerationTopWall = yVelocityTopWall
    !     zAccelerationTopWall = zVelocityTopWall
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Compute wall acceleration at top and bottom wall
    ! subroutine compute_acceleration_bottom_wall(xVelocityBottomWall, yVelocityBottomWall, zVelocityBottomWall, deltaT)
    !     complex(kind=sp), intent(in) :: xVelocityBottomWall(:,:), yVelocityBottomWall(:,:), zVelocityBottomWall(:,:)
    !     real(kind=sp) :: deltaT
    !     ! note: we assume the acceleration arrays currently hold the velocity at the previous time step (see subroutine
    !     ! save_velocity_bottom_wall_to_buffer). We now compute the acceleration through finite differences
    !     xAccelerationBottomWall = (xVelocityBottomWall - xAccelerationBottomWall) / deltaT
    !     yAccelerationBottomWall = (yVelocityBottomWall - yAccelerationBottomWall) / deltaT
    !     zAccelerationBottomWall = (zVelocityBottomWall - zAccelerationBottomWall) / deltaT
    ! end subroutine

    ! subroutine compute_acceleration_top_wall(xVelocityTopWall, yVelocityTopWall, zVelocityTopWall, deltaT)
    !     complex(kind=sp), intent(in) :: xVelocityTopWall(:,:), yVelocityTopWall(:,:), zVelocityTopWall(:,:)
    !     real(kind=sp), intent(in) :: deltaT
    !     ! note: we assume the acceleration arrays currently hold the velocity at the previous time step (see subroutine
    !     ! save_velocity_top_wall_to_buffer). We now compute the acceleration through finite differences
    !     xAccelerationTopWall = (xVelocityTopWall - xAccelerationTopWall) / deltaT
    !     yAccelerationTopWall = (yVelocityTopWall - yAccelerationTopWall) / deltaT
    !     zAccelerationTopWall = (zVelocityTopWall - zAccelerationTopWall) / deltaT
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Write Data to HDF file
    subroutine write_flowfield_to_hdf_file(fundamentalWavenumberX, fundamentalWavenumberZ, yCoordinateVector, &
        indexStartingPlaneProcessor, indexEndingPlaneProcessor, wavenumberVectorX, wavenumberVectorZ, time, &
        referenceReynoldsNumber, bulkVelocity)
        real(kind=sp), intent(in) :: fundamentalWavenumberX, fundamentalWavenumberZ
        real(kind=sp), intent(in) :: yCoordinateVector(:)
        integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
        complex(kind=sp), intent(in) :: wavenumberVectorX(:), wavenumberVectorZ(:)
        real(kind=sp), intent(in) :: time, referenceReynoldsNumber
        real(kind=dp), intent(in) :: bulkVelocity
        integer :: hdfError
        integer(kind=hid_t) :: idFile
        call h5open_f(hdfError)
        call open_parallel_hdf_file(idFile)
        ! Write coordinate and wavenumber
        call write_coordinate_vectors_to_hdf_file(idFile, fundamentalWavenumberX, fundamentalWavenumberZ, &
            yCoordinateVector)
        call write_wavenumber_vectors_to_hdf_file(idFile, wavenumberVectorX, wavenumberVectorZ)
        ! Write Fourier data
        call write_velocity_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
            referenceReynoldsNumber, bulkVelocity)
        call write_vorticity_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
            referenceReynoldsNumber, bulkVelocity)
        ! call write_phi_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        !     referenceReynoldsNumber, bulkVelocity)
        ! call write_yderivative_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        !     referenceReynoldsNumber, bulkVelocity)
        ! call write_H_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        !     referenceReynoldsNumber, bulkVelocity)
        ! call write_forcing_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        !     referenceReynoldsNumber, bulkVelocity)
        ! Write Physical data
        ! call write_velocity_fields_to_hdf_file_phys(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        ! referenceReynoldsNumber, bulkVelocity)
        ! call write_vorticity_fields_to_hdf_file_phys(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        ! referenceReynoldsNumber, bulkVelocity)
        ! Write Acceleration data
        ! call write_wall_acceleration_to_hdf_file(idFile)
        ! Complete
        call close_parallel_hdf_file(idFile)
        call h5close_f(hdfError)
        call update_module_variables
    end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Open and write coordinate and wavenumbers
    subroutine open_parallel_hdf_file(idFile)
        integer(kind=hid_t), intent(out) :: idFile
        integer(kind=hid_t) :: idPropertyList
        integer :: hdfError
        character(len=4) :: fileNumberString
        character(:), allocatable :: fileName
        call h5pcreate_f(H5P_FILE_ACCESS_F, idPropertyList, hdfError)
        call h5pset_fapl_mpio_f(idPropertyList, MPI_COMM_WORLD, MPI_INFO_NULL, hdfError)
        write(fileNumberString, '(i4.4)') fileNumber
        fileName = baseFileNameWithPath // fileNumberString // '.h5'
        call h5fcreate_f(fileName, H5F_ACC_TRUNC_F, idFile, hdfError, access_prp=idPropertyList)
        call h5pclose_f(idPropertyList, hdfError)
    end subroutine

    subroutine write_coordinate_vectors_to_hdf_file(idFile, fundamentalWavenumberX, fundamentalWavenumberZ, &
        yCoordinateVector)
        integer(kind=hid_t), intent(in) :: idFile
        real(kind=sp), intent(in) :: fundamentalWavenumberX, fundamentalWavenumberZ
        real(kind=sp), intent(in) :: yCoordinateVector(:)
        integer(kind=hid_t) :: idGroup
        character(:), allocatable :: groupName
        integer :: hdfError
        groupName = 'coordinateVectors' 
        call h5gcreate_f(idFile, groupName, idGroup, hdfError)
        ! flush the file before writing the coordinate vectors. Not flushing leads to corrupted hdf files
        ! (this seems to be required only before individual writes of single processors)
        call h5fflush_f(idFile, H5F_SCOPE_GLOBAL_F, hdfError)
        call write_streamwise_coordinate_vector_to_hdf_file(fundamentalWavenumberX, idGroup)
        call write_spanwise_coordinate_vector_to_hdf_file(fundamentalWavenumberZ, idGroup)
        call write_wallnormal_coordinate_vector_to_hdf_file(yCoordinateVector, idGroup)
        call h5gclose_f(idGroup, hdfError)
    end subroutine

    subroutine write_wavenumber_vectors_to_hdf_file(idFile, wavenumberVectorX, wavenumberVectorZ)
        integer(kind=hid_t), intent(in) :: idFile
        complex(kind=sp), intent(in) :: wavenumberVectorX(:), wavenumberVectorZ(:)
        integer(kind=hid_t) :: idGroup
        character(:), allocatable :: groupName
        integer :: hdfError
        groupName = 'wavenumberVectors'
        call h5gcreate_f(idFile, groupName, idGroup, hdfError)
        ! flush the file before writing the coordinate vectors. Not flushing leads to corrupted hdf files
        ! (this seems to be required only before individual writes of single processors)
        call h5fflush_f(idFile, H5F_SCOPE_GLOBAL_F, hdfError)
        ! the dns code stores the wavenumbers as purely imaginary vector,
        ! i.e. wavenumberVectorX(k) = (0,k*fundamentalWavenumberX)
        ! in hdf5 we store the wavevector as real array, so we need to extract the imaginary part
        call write_1d_array_to_hdf_file(aimag(wavenumberVectorX), 'xWavenumberVector', idGroup)
        call write_1d_array_to_hdf_file(aimag(wavenumberVectorZ), 'zWavenumberVector', idGroup)
        call h5gclose_f(idGroup, hdfError)
    end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Write Fourier Data to HDF File
    subroutine write_velocity_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        referenceReynoldsNumber, bulkVelocity)
        integer(kind=hid_t), intent(in) :: idFile
        integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
        real(kind=sp), intent(in) :: time, referenceReynoldsNumber
        real(kind=dp), intent(in) :: bulkVelocity
        integer(kind=hid_t) :: idGroup
        character(:), allocatable :: groupName
        integer :: hdfError
        groupName = 'velocityFieldsFourier'
        call h5gcreate_f(idFile, groupName, idGroup, hdfError)
        ! the following write is parallel, therefore no flush is required
        ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each velocity
        ! component as separate data set.
        ! streamwise velocity
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferStreamwiseVelocity), 'uRealPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferStreamwiseVelocity), 'uImaginaryPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        ! wallnormal velocity
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferWallnormalVelocity), 'vRealPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferWallnormalVelocity), 'vImaginaryPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        ! spanwise velocity
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferSpanwiseVelocity), 'wRealPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferSpanwiseVelocity), 'wImaginaryPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call h5gclose_f(idGroup, hdfError)
    end subroutine

    subroutine write_vorticity_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
        referenceReynoldsNumber, bulkVelocity)
        integer(kind=hid_t), intent(in) :: idFile
        integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
        real(kind=sp), intent(in) :: time, referenceReynoldsNumber
        real(kind=dp), intent(in) :: bulkVelocity
        integer(kind=hid_t) :: idGroup
        character(:), allocatable :: groupName
        integer :: hdfError
        groupName = 'vorticityFieldsFourier'
        call h5gcreate_f(idFile, groupName, idGroup, hdfError)
        ! the following write is parallel, therefore no flush is required
        ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each vorticity
        ! component as separate data set.
        ! streamwise vorticity
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferStreamwiseVorticity), 'o1RealPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferStreamwiseVorticity), 'o1ImaginaryPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        ! wallnormal vorticity
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferWallnormalVorticity), 'o2RealPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferWallnormalVorticity), 'o2ImaginaryPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        ! spanwise vorticity
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferSpanwiseVorticity), 'o3RealPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferSpanwiseVorticity), 'o3ImaginaryPart', &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        call h5gclose_f(idGroup, hdfError)
    end subroutine

    ! subroutine write_phi_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
    !     referenceReynoldsNumber, bulkVelocity)
    !     integer(kind=hid_t), intent(in) :: idFile
    !     integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
    !     real(kind=sp), intent(in) :: time, referenceReynoldsNumber
    !     real(kind=dp), intent(in) :: bulkVelocity
    !     integer(kind=hid_t) :: idGroup
    !     character(:), allocatable :: groupName
    !     integer :: hdfError
    !     groupName = 'phiFieldsFourier'
    !     call h5gcreate_f(idFile, groupName, idGroup, hdfError)
    !     ! the following write is parallel, therefore no flush is required
    !     ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each vorticity
    !     ! component as separate data set.
    !     ! phi
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferPhi), 'phiRealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferPhi), 'phiImaginaryPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call h5gclose_f(idGroup, hdfError)
    ! end subroutine

    ! subroutine write_yderivative_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
    !     referenceReynoldsNumber, bulkVelocity)
    !     integer(kind=hid_t), intent(in) :: idFile
    !     integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
    !     real(kind=sp), intent(in) :: time, referenceReynoldsNumber
    !     real(kind=dp), intent(in) :: bulkVelocity
    !     integer(kind=hid_t) :: idGroup
    !     character(:), allocatable :: groupName
    !     integer :: hdfError
    !     groupName = 'yderivativeFieldsFourier'
    !     call h5gcreate_f(idFile, groupName, idGroup, hdfError)
    !     ! the following write is parallel, therefore no flush is required
    !     ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each velocity
    !     ! component as separate data set.
    !     ! dvdy
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferdvdy), 'dvdyRealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferdvdy), 'dvdyImaginaryPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! do2dy
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferdo2dy), 'do2dyRealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferdo2dy), 'do2dyImaginaryPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call h5gclose_f(idGroup, hdfError)
    ! end subroutine

    ! subroutine write_H_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
    !     referenceReynoldsNumber, bulkVelocity)
    !     integer(kind=hid_t), intent(in) :: idFile
    !     integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
    !     real(kind=sp), intent(in) :: time, referenceReynoldsNumber
    !     real(kind=dp), intent(in) :: bulkVelocity
    !     integer(kind=hid_t) :: idGroup
    !     character(:), allocatable :: groupName
    !     integer :: hdfError
    !     groupName = 'HFieldsFourier'
    !     call h5gcreate_f(idFile, groupName, idGroup, hdfError)
    !     ! the following write is parallel, therefore no flush is required
    !     ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each vorticity
    !     ! component as separate data set.
    !     ! H1
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferH1), 'H1RealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferH1), 'H1ImaginaryPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! H2
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferH2), 'H2RealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferH2), 'H2ImaginaryPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! H3
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferH3), 'H3RealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferH3), 'H3ImaginaryPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call h5gclose_f(idGroup, hdfError)
    ! end subroutine

    ! subroutine write_forcing_fields_to_hdf_file(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
    !     referenceReynoldsNumber, bulkVelocity)
    !     integer(kind=hid_t), intent(in) :: idFile
    !     integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
    !     real(kind=sp), intent(in) :: time, referenceReynoldsNumber
    !     real(kind=dp), intent(in) :: bulkVelocity
    !     integer(kind=hid_t) :: idGroup
    !     character(:), allocatable :: groupName
    !     integer :: hdfError
    !     groupName = 'forcingFieldsFourier'
    !     call h5gcreate_f(idFile, groupName, idGroup, hdfError)
    !     ! the following write is parallel, therefore no flush is required
    !     ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each forcing
    !     ! component as separate data set.
    !     ! forcing wall-normal velocity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferForcingVelocity), 'vForcingRealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferForcingVelocity), 'vForcingImaginaryPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! forcing wall-normal vorticity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, real(bufferForcingVorticity), 'o2ForcingRealPart', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, aimag(bufferForcingVorticity), &
    !         'o2ForcingImaginaryPart', indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
    !         referenceReynoldsNumber, bulkVelocity)
    !     call h5gclose_f(idGroup, hdfError)
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Write Physical Data to HDF File
    ! subroutine write_velocity_fields_to_hdf_file_phys(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
    !     referenceReynoldsNumber, bulkVelocity)
    !     integer(kind=hid_t), intent(in) :: idFile
    !     integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
    !     real(kind=sp), intent(in) :: time, referenceReynoldsNumber
    !     real(kind=dp), intent(in) :: bulkVelocity
    !     integer(kind=hid_t) :: idGroup
    !     character(:), allocatable :: groupName
    !     integer :: hdfError
    !     groupName = 'velocityFieldsPhysical'
    !     call h5gcreate_f(idFile, groupName, idGroup, hdfError)
    !     ! the following write is parallel, therefore no flush is required
    !     ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each velocity
    !     ! component as separate data set.
    !     ! streamwise velocity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, bufferUPhys, 'uPhys', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! wallnormal velocity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, bufferVPhys, 'vPhys', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! spanwise velocity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, bufferWPhys, 'wPhys', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call h5gclose_f(idGroup, hdfError)
    ! end subroutine

    ! subroutine write_vorticity_fields_to_hdf_file_phys(idFile, indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, &
    !     referenceReynoldsNumber, bulkVelocity)
    !     integer(kind=hid_t), intent(in) :: idFile
    !     integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
    !     real(kind=sp), intent(in) :: time, referenceReynoldsNumber
    !     real(kind=dp), intent(in) :: bulkVelocity
    !     integer(kind=hid_t) :: idGroup
    !     character(:), allocatable :: groupName
    !     integer :: hdfError
    !     groupName = 'vorticityFieldsPhysical'
    !     call h5gcreate_f(idFile, groupName, idGroup, hdfError)
    !     ! the following write is parallel, therefore no flush is required
    !     ! hdf5 has no built-in complex datatype. We therefore write the real and imaginary part of each vorticity
    !     ! component as separate data set.
    !     ! streamwise vorticity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, buffero1Phys, 'o1Phys', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! wallnormal vorticity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, buffero2Phys, 'o2Phys', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     ! spanwise vorticity
    !     call write_real_buffer_with_attributes_to_hdf_file(idGroup, buffero3Phys, 'o3Phys', &
    !         indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
    !     call h5gclose_f(idGroup, hdfError)
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Write Acceleration Data to HDF File
    ! subroutine write_wall_acceleration_to_hdf_file(idFile)
    !     integer(kind=hid_t), intent(in) :: idFile
    !     integer(kind=hid_t) :: idGroup
    !     character(:), allocatable :: groupName
    !     integer :: hdfError
    !     groupName = 'wallAccelerationsFourier'
    !     call h5gcreate_f(idFile, groupName, idGroup, hdfError)
    !     ! flush the file before writing the wall acceleration. Not flushing leads to corrupted hdf files.
    !     ! (this seems to be required only before individual writes of a single processor)
    !     call h5fflush_f(idFile, H5F_SCOPE_GLOBAL_F, hdfError)
    !     call write_acceleration_bottom_wall_to_hdf_file(idGroup)
    !     call write_acceleration_top_wall_to_hdf_file(idGroup)
    !     call h5gclose_f(idGroup, hdfError)
    ! end subroutine

    ! subroutine write_acceleration_bottom_wall_to_hdf_file(idFileOrGroup)
    !     integer(kind=hid_t), intent(in) :: idFileOrGroup
    !     call write_real_2d_buffer_to_hdf_file(real(xAccelerationBottomWall), 'xAccelerationBottomRealPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(aimag(xAccelerationBottomWall), 'xAccelerationBottomImaginaryPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(real(yAccelerationBottomWall), 'yAccelerationBottomRealPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(aimag(yAccelerationBottomWall), 'yAccelerationBottomImaginaryPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(real(zAccelerationBottomWall), 'zAccelerationBottomRealPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(aimag(zAccelerationBottomWall), 'zAccelerationBottomImaginaryPart', idFileOrGroup)
    ! end subroutine

    ! subroutine write_acceleration_top_wall_to_hdf_file(idFileOrGroup)
    !     integer(kind=hid_t), intent(in) :: idFileOrGroup
    !     call write_real_2d_buffer_to_hdf_file(real(xAccelerationTopWall), 'xAccelerationTopRealPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(aimag(xAccelerationTopWall), 'xAccelerationTopImaginaryPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(real(yAccelerationTopWall), 'yAccelerationTopRealPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(aimag(yAccelerationTopWall), 'yAccelerationTopImaginaryPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(real(zAccelerationTopWall), 'zAccelerationTopRealPart', idFileOrGroup)
    !     call write_real_2d_buffer_to_hdf_file(aimag(zAccelerationTopWall), 'zAccelerationTopImaginaryPart', idFileOrGroup)
    ! end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Close Files
    subroutine close_parallel_hdf_file(idFile)
        integer(kind=hid_t), intent(in) :: idFile
        integer :: hdfError
        call h5fclose_f(idFile, hdfError)
    end subroutine

    subroutine update_module_variables()
        fileNumber = fileNumber + 1
        collectFlowfield = .false.
    end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Low level write data to HDF file
    subroutine write_streamwise_coordinate_vector_to_hdf_file(fundamentalWavenumberX, idFileOrGroup)
        ! the streamwise coordinate is denoted by x below
        real(kind=sp), intent(in) :: fundamentalWavenumberX
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        real(kind=sp), parameter :: pi = acos(-1.0_sp)
        real(kind=sp) :: boxLengthX, deltaX
        real(kind=sp) :: xCoordinateVector(mgalx)
        integer :: i
        boxLengthX = 2.0_sp * pi / fundamentalWavenumberX
        deltaX = boxLengthX / mgalx
        do i = 1, mgalx
            xCoordinateVector(i) = (i-1) * deltaX
        enddo
        call write_1d_array_to_hdf_file(xCoordinateVector, 'xCoordinateVector', idFileOrGroup)
    end subroutine

    subroutine write_spanwise_coordinate_vector_to_hdf_file(fundamentalWavenumberZ, idFileOrGroup)
        ! the spanwise coordinate is denoted by z below
        real(kind=sp), intent(in) :: fundamentalWavenumberZ
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        real(kind=sp), parameter :: pi = acos(-1.0_sp)
        real(kind=sp) :: boxLengthZ, deltaZ
        real(kind=sp) :: zCoordinateVector(mgalz)
        integer :: i
        boxLengthZ = 2.0_sp * pi / fundamentalWavenumberZ
        deltaZ = boxLengthZ / mgalz
        do i = 1, mgalz
            zCoordinateVector(i) = (i-1) * deltaZ
        enddo
        call write_1d_array_to_hdf_file(zCoordinateVector, 'zCoordinateVector', idFileOrGroup)
    end subroutine

    subroutine write_wallnormal_coordinate_vector_to_hdf_file(yCoordinateVector, idFileOrGroup)
        ! the wallnormal coordinate is denoted by y below
        real(kind=sp), intent(in) :: yCoordinateVector(:)
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        call write_1d_array_to_hdf_file(yCoordinateVector, 'yCoordinateVector', idFileOrGroup)
    end subroutine

    subroutine write_1d_array_to_hdf_file(dataArray, nameDataArray, idFileOrGroup)
        real(kind=sp), intent(in) :: dataArray(:)
        character(len=*), intent(in) :: nameDataArray
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        integer, parameter :: rank = 1
        integer(kind=hsize_t) :: dataDimensions(rank)
        integer(kind=hid_t) :: idDataSpace, idDataSet
        integer :: hdfError, mpiRank, mpiError
        integer, parameter :: mpiMaster = 0
        dataDimensions(1) = size(dataArray)
        ! for a parallel hdf file, dataspace and dataset creation are collective functions, i.e. all processes
        ! of the communicator must call them, even if some processors don't write
        call h5screate_simple_f(rank, dataDimensions, idDataSpace, hdfError)
        call h5dcreate_f(idFileOrGroup, nameDataArray, H5T_IEEE_F32LE, idDataSpace, idDataSet, hdfError)
        call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)
        ! the data transfer can be independent, here only the master writes data
        if (mpiRank == mpiMaster) then
            call h5dwrite_f(idDataSet, H5T_NATIVE_REAL, dataArray, dataDimensions, hdfError)
        endif
        call h5dclose_f(idDataSet, hdfError)
        call h5sclose_f(idDataSpace, hdfError)
        ! flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files.
        call h5fflush_f(idFileOrGroup, H5F_SCOPE_GLOBAL_F, hdfError)
    end subroutine

    subroutine write_real_buffer_with_attributes_to_hdf_file(idFileOrGroup, dataBuffer, nameDataBuffer, &
        indexStartingPlaneProcessor, indexEndingPlaneProcessor, time, referenceReynoldsNumber, bulkVelocity)
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        character(len=*), intent(in) :: nameDataBuffer
        integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
        real(kind=sp), intent(in) :: time, referenceReynoldsNumber
        real(kind=dp), intent(in) :: bulkVelocity
        call write_real_buffer_to_hdf_file(idFileOrGroup, dataBuffer, nameDataBuffer, &
            indexStartingPlaneProcessor, indexEndingPlaneProcessor)
        call add_single_precision_attribute_to_data_set(idFileOrGroup, nameDataBuffer, time, 'time')
        call add_single_precision_attribute_to_data_set(idFileOrGroup, nameDataBuffer, &
            referenceReynoldsNumber, 'referenceReynoldsNumber')
        call add_double_precision_attribute_to_data_set(idFileOrGroup, nameDataBuffer, bulkVelocity, 'bulkVelocity')
        call add_string_attribute_to_data_set(idFileOrGroup, nameDataBuffer, &
            'Data order of array in h5 file is (yCoordinateVector, zWavenumberVector, xWavenumberVector)', 'dataLayout')
        call add_string_attribute_to_data_set(idFileOrGroup, nameDataBuffer, &
            'The data set only contains non-negative streamwise wavenumbers. ' // &
            'The negative streamwise wavenumbers can be reconstructed from Hermitian symmetry.', &
            'negativeStreamwiseWavenumbers')
    end subroutine

    subroutine write_real_buffer_to_hdf_file(idFileOrGroup, dataBuffer, nameDataBuffer, &
        indexStartingPlaneProcessor, indexEndingPlaneProcessor)
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        character(len=*), intent(in) :: nameDataBuffer
        integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
        integer, parameter :: rank = 3
        integer(kind=hsize_t) :: dimDataSet(rank), slabOffset(rank), numberOfBlocks(rank), dimMemorySpace(rank)
        integer(kind=hid_t) :: idDataSpace, idDataSet, idMemorySpace, idPropertyList
        integer(kind=hid_t) :: idDiskSpace
        integer :: hdfError
        call create_disk_dataspace(idDiskSpace, dataBuffer)
        call select_disk_hyperslab(idDiskSpace, dataBuffer, indexStartingPlaneProcessor, indexEndingPlaneProcessor)
        call create_memory_dataspace(idMemorySpace, dataBuffer)
        call select_memory_hyperslab(idMemorySpace, dataBuffer, indexStartingPlaneProcessor, indexEndingPlaneProcessor)
        call transfer_data_from_memory_to_disk(idFileOrGroup, idDiskSpace, idMemorySpace, dataBuffer, nameDataBuffer)
        call close_disk_dataspace(idDiskSpace)
        call close_memory_dataspace(idMemorySpace)
        ! flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(idFileOrGroup, H5F_SCOPE_GLOBAL_F, hdfError)
    end subroutine

    subroutine write_real_2d_buffer_to_hdf_file(dataArray, nameDataArray, idFileOrGroup)
        real(kind=sp), intent(in) :: dataArray(:,:)
        character(len=*), intent(in) :: nameDataArray
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        integer, parameter :: rank = 2
        integer(kind=hsize_t) :: dataDimensions(rank)
        integer(kind=hid_t) :: idDataSpace, idDataSet
        integer :: hdfError, mpiRank, mpiError
        integer, parameter :: mpiMaster = 0
        dataDimensions(1) = size(dataArray, 1)
        dataDimensions(2) = size(dataArray, 2)
        ! for a parallel hdf file, dataspace and dataset creation are collective functions, i.e. all processes
        ! of the communicator must call them, even if some processors don't write
        call h5screate_simple_f(rank, dataDimensions, idDataSpace, hdfError)
        call h5dcreate_f(idFileOrGroup, nameDataArray, H5T_IEEE_F32LE, idDataSpace, idDataSet, hdfError)
        call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)
        ! the data transfer can be independent, here only the master writes data
        if (mpiRank == mpiMaster) then
            call h5dwrite_f(idDataSet, H5T_NATIVE_REAL, dataArray, dataDimensions, hdfError)
        endif
        call h5dclose_f(idDataSet, hdfError)
        call h5sclose_f(idDataSpace, hdfError)
        ! flush buffer to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(idFileOrGroup, H5F_SCOPE_GLOBAL_F, hdfError)
    end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! HDF Attributes
    subroutine add_single_precision_attribute_to_data_set(idFileOrGroup, dataSetName, attribute, attributeName)
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        character(len=*), intent(in) :: dataSetName, attributeName
        real(kind=sp), intent(in) :: attribute
        integer(kind=hid_t) :: idDataSet, idAttributeSpace, idAttribute
        integer :: hdfError
        integer(kind=hsize_t) :: dimAttribute(1)
        call h5dopen_f(idFileOrGroup, dataSetName, idDataSet, hdfError)
        call h5screate_f(H5S_SCALAR_F, idAttributeSpace, hdfError)
        call h5acreate_f(idDataSet, attributeName, H5T_IEEE_F32LE, idAttributeSpace, idAttribute, hdfError)
        ! since we created a scalar dataspace, dimAttribute is ignored by h5awrite_f, but we still need to pass it
        dimAttribute(1) = 1
        call h5awrite_f(idAttribute, H5T_NATIVE_REAL, attribute, dimAttribute, hdfError)
        call h5aclose_f(idAttribute, hdfError)
        call h5sclose_f(idAttributeSpace, hdfError)
        call h5dclose_f(idDataSet, hdfError)
    end subroutine

    subroutine add_double_precision_attribute_to_data_set(idFileOrGroup, dataSetName, attribute, attributeName)
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        character(len=*), intent(in) :: dataSetName, attributeName
        real(kind=dp), intent(in) :: attribute
        integer(kind=hid_t) :: idDataSet, idAttributeSpace, idAttribute
        integer :: hdfError
        integer(kind=hsize_t) :: dimAttribute(1)
        call h5dopen_f(idFileOrGroup, dataSetName, idDataSet, hdfError)
        call h5screate_f(H5S_SCALAR_F, idAttributeSpace, hdfError)
        call h5acreate_f(idDataSet, attributeName, H5T_IEEE_F64LE, idAttributeSpace, idAttribute, hdfError)
        ! since we created a scalar dataspace, dimAttribute is ignored by h5awrite_f, but we still need to pass it
        dimAttribute(1) = 1
        call h5awrite_f(idAttribute, H5T_NATIVE_DOUBLE, attribute, dimAttribute, hdfError)
        call h5aclose_f(idAttribute, hdfError)
        call h5sclose_f(idAttributeSpace, hdfError)
        call h5dclose_f(idDataSet, hdfError)
    end subroutine

    subroutine add_string_attribute_to_data_set(idFileOrGroup, dataSetName, stringAttribute, attributeName)
        integer(kind=hid_t), intent(in) :: idFileOrGroup
        character(len=*), intent(in) :: dataSetName, attributeName
        character(len=*), intent(in) :: stringAttribute
        integer(kind=hid_t) :: idDataSet, idAttributeSpace, idDataType, idAttribute
        integer(kind=hsize_t) :: dimAttribute(1)
        integer :: hdfError
        call h5dopen_f(idFileOrGroup, dataSetName, idDataSet, hdfError)
        ! to deal with the variable length of the string attribute, we create a scalar dataspace and define a new
        ! character datatype of appropriate length
        call h5screate_f(H5S_SCALAR_F, idAttributeSpace, hdfError)
        call h5tcopy_f(H5T_FORTRAN_S1, idDataType, hdfError)
        call h5tset_size_f(idDataType, len(stringAttribute, kind=size_t), hdfError)
        call h5acreate_f(idDataSet, attributeName, idDataType, idAttributeSpace, idAttribute, hdfError)
        ! since we created a scalar dataspace, dimAttribute is ignored by h5awrite_f, but we still need to pass it
        dimAttribute(1) = 1
        call h5awrite_f(idAttribute, idDataType, stringAttribute, dimAttribute, hdfError)
        call h5aclose_f(idAttribute, hdfError)
        call h5tclose_f(idDataType, hdfError)
        call h5sclose_f(idAttributeSpace, hdfError)
        call h5dclose_f(idDataSet, hdfError)
    end subroutine
    ! ***********************************************************************************************************************************


    ! ***********************************************************************************************************************************
    ! Manage HDF file
    subroutine create_disk_dataspace(idDiskSpace, dataBuffer)
        integer(kind=hid_t), intent(out) :: idDiskSpace
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        integer, parameter :: rank = 3
        integer(kind=hsize_t) :: dimDiskSpace(rank)
        integer :: hdfError
        ! size of the dataset on disk corresponds to the size of the computational domain. The size of the domain in
        ! the streamwise (1st) and spanwise (2nd) direction depends on whether real or complex data is written, but
        ! in both cases, the databuffer has the size of the domain.
        ! The size in the wallnormal (3rd) direction is always my. Taking this dimension from the databuffer would be
        ! incorrect, sice each processor only stores a subset of the wallnormal planes (of size myp instead of my)
        dimDiskSpace(1) = size(dataBuffer, 1)
        dimDiskSpace(2) = size(dataBuffer, 2)
        dimDiskSpace(3) = my
        call h5screate_simple_f(rank, dimDiskSpace, idDiskSpace, hdfError)
    end subroutine

    subroutine select_disk_hyperslab(idDiskSpace, dataBuffer, indexStartingPlaneProcessor, indexEndingPlaneProcessor)
        ! select that hyperlsab of the disk dataspace, which corresponds to the processor's slice
        ! of the computational domain
        integer(kind=hid_t), intent(in) :: idDiskSpace
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
        integer, parameter :: rank = 3
        integer(kind=hsize_t) :: slabOffset(rank), numberOfBlocks(rank)
        integer :: hdfError
        slabOffset(1) = 0
        slabOffset(2) = 0
        slabOffset(3) = indexStartingPlaneProcessor - 1
        numberOfBlocks(1) = size(dataBuffer, 1)
        numberOfBlocks(2) = size(dataBuffer, 2)
        numberOfBlocks(3) = indexEndingPlaneProcessor - indexStartingPlaneProcessor + 1
        call h5sselect_hyperslab_f(idDiskSpace, H5S_SELECT_SET_F, slabOffset, numberOfBlocks, hdfError)
    end subroutine

    subroutine create_memory_dataspace(idMemorySpace, dataBuffer)
        integer(kind=hid_t), intent(out) :: idMemorySpace
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        integer, parameter :: rank = 3
        integer(kind=hsize_t) :: dimMemorySpace(rank)
        integer :: hdfError
        ! the size of the dataspace in memory corresponds to the size of the procesor's slice of the domain
        dimMemorySpace(1) = size(dataBuffer, 1)
        dimMemorySpace(2) = size(dataBuffer, 2)
        dimMemorySpace(3) = size(dataBuffer, 3)
        call h5screate_simple_f(rank, dimMemorySpace, idMemorySpace, hdfError)
    end subroutine

    subroutine select_memory_hyperslab(idMemorySpace, dataBuffer, indexStartingPlaneProcessor, indexEndingPlaneProcessor)
        integer(kind=hid_t), intent(in) :: idMemorySpace
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        integer, intent(in) :: indexStartingPlaneProcessor, indexEndingPlaneProcessor
        integer, parameter :: rank = 3
        integer(kind=hsize_t) :: slabOffset(rank), numberOfBlocks(rank)
        integer :: hdfError
        slabOffset(1) = 0
        slabOffset(2) = 0
        slabOffset(3) = 0
        numberOfBlocks(1) = size(dataBuffer, 1)
        numberOfBlocks(2) = size(dataBuffer, 2)
        numberOfBlocks(3) = indexEndingPlaneProcessor - indexStartingPlaneProcessor + 1
        call h5sselect_hyperslab_f(idMemorySpace, H5S_SELECT_SET_F, slabOffset, numberOfBlocks, hdfError)
    end subroutine

    subroutine transfer_data_from_memory_to_disk(idFileOrGroup, idDiskSpace, idMemorySpace, dataBuffer, nameDataBuffer)
        integer(kind=hid_t), intent(in) :: idFileOrGroup, idDiskSpace, idMemorySpace
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        character(len=*), intent(in) :: nameDataBuffer
        integer(kind=hid_t) :: idDataSet, idPropertyList
        integer :: hdfError
        integer, parameter :: rank = 3
        integer(kind=hsize_t) :: dimDataBuffer(rank)
        call h5dcreate_f(idFileOrGroup, nameDataBuffer, H5T_IEEE_F32LE, idDiskSpace, idDataSet, hdfError)
        call h5pcreate_f(H5P_DATASET_XFER_F, idPropertyList, hdfError)
        call h5pset_dxpl_mpio_f(idPropertyList, H5FD_MPIO_COLLECTIVE_F, hdfError)
        dimDataBuffer(1) = size(dataBuffer, 1)
        dimDataBuffer(2) = size(dataBuffer, 2)
        dimDataBuffer(3) = size(dataBuffer, 3)
        call h5dwrite_f(idDataSet, H5T_NATIVE_REAL, dataBuffer, dimDataBuffer, hdfError, &
                mem_space_id=idMemorySpace, file_space_id=idDiskSpace, xfer_prp=idPropertyList)
        call h5pclose_f(idPropertyList, hdfError)
        call h5dclose_f(idDataSet, hdfError)
    end subroutine

    subroutine close_disk_dataspace(idDiskSpace)
        integer(kind=hid_t), intent(in) :: idDiskSpace
        integer :: hdfError
        call h5sclose_f(idDiskSpace, hdfError)
    end subroutine

    subroutine close_memory_dataspace(idMemorySpace)
        integer(kind=hid_t), intent(in) :: idMemorySpace
        integer :: hdfError
        call h5sclose_f(idMemorySpace, hdfError)
    end subroutine
    ! ***********************************************************************************************************************************

end module save_flowfield
