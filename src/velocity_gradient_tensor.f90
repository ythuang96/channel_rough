module velocity_gradient_tensor
    use types, only: sp
    implicit none
    private
    include 'ctes3D'
    logical :: collectData = .false.
    integer, parameter :: collectionFrequencyInTimesteps = 5000
    integer :: fileNrVisualization = 1
    real(kind=sp), allocatable :: dudx(:,:,:), dudy(:,:,:), dudz(:,:,:), &
                                  dvdx(:,:,:), dvdy(:,:,:), dvdz(:,:,:), &
                                  dwdx(:,:,:), dwdy(:,:,:), dwdz(:,:,:), &
                                  vSlice(:,:,:), vEntireDomain(:,:,:)
    public :: collectData, allocate_data_buffers, deallocate_data_buffers, &
              compute_velocity_gradient_tensor, write_velocity_gradient_tensor, &
              save_velocity_fields, write_velocity_fields, &
              collectionFrequencyInTimesteps, fileNrVisualization

contains
    subroutine allocate_data_buffers(startingPlaneIndex, endingPlaneIndex, mpiRank)
        integer, intent(in) :: startingPlaneIndex, endingPlaneIndex, mpiRank
        integer, parameter :: mpiMaster = 0
        integer :: nPlanes
        nPlanes = endingPlaneIndex - startingPlaneIndex + 1
        allocate(dudx(mgalx, mgalz, nPlanes))
        if (.not. allocated(dudx)) stop 'Error in allocate_data_buffers'
        allocate(dudy(mgalx, mgalz, nPlanes))
        if (.not. allocated(dudy)) stop 'Error in allocate_data_buffers'
        allocate(dudz(mgalx, mgalz, nPlanes))
        if (.not. allocated(dudz)) stop 'Error in allocate_data_buffers'
        allocate(dvdx(mgalx, mgalz, nPlanes))
        if (.not. allocated(dvdx)) stop 'Error in allocate_data_buffers'
        allocate(dvdy(mgalx, mgalz, nPlanes))
        if (.not. allocated(dvdy)) stop 'Error in allocate_data_buffers'
        allocate(dvdz(mgalx, mgalz, nPlanes))
        if (.not. allocated(dvdz)) stop 'Error in allocate_data_buffers'
        allocate(dwdx(mgalx, mgalz, nPlanes))
        if (.not. allocated(dwdx)) stop 'Error in allocate_data_buffers'
        allocate(dwdy(mgalx, mgalz, nPlanes))
        if (.not. allocated(dwdy)) stop 'Error in allocate_data_buffers'
        allocate(dwdz(mgalx, mgalz, nPlanes))
        if (.not. allocated(dwdz)) stop 'Error in allocate_data_buffers'
        allocate(vSlice(mgalx, mgalz, nPlanes))
        if (.not. allocated(vSlice)) stop 'Error in allocate_data_buffers'
        if (mpiRank == mpiMaster) then
            allocate(vEntireDomain(mgalx, mgalz, my))
            if (.not. allocated(vEntireDomain)) stop 'Error in allocate_data_buffers'
        endif
    end subroutine allocate_data_buffers

    subroutine deallocate_data_buffers(mpiRank)
        integer, intent(in) :: mpiRank
        integer, parameter :: mpiMaster = 0
        deallocate(dudx)
        if (allocated(dudx)) stop 'Error in deallocate_data_buffers'
        deallocate(dudy)
        if (allocated(dudy)) stop 'Error in deallocate_data_buffers'
        deallocate(dudz)
        if (allocated(dudz)) stop 'Error in deallocate_data_buffers'
        deallocate(dvdx)
        if (allocated(dvdx)) stop 'Error in deallocate_data_buffers'
        deallocate(dvdy)
        if (allocated(dvdy)) stop 'Error in deallocate_data_buffers'
        deallocate(dvdz)
        if (allocated(dvdz)) stop 'Error in deallocate_data_buffers'
        deallocate(dwdx)
        if (allocated(dwdx)) stop 'Error in deallocate_data_buffers'
        deallocate(dwdy)
        if (allocated(dwdy)) stop 'Error in deallocate_data_buffers'
        deallocate(dwdz)
        if (allocated(dwdz)) stop 'Error in deallocate_data_buffers'
        deallocate(vSlice)
        if (allocated(vSlice)) stop 'Error in deallocate_data_buffers'
        if (mpiRank == mpiMaster) then
            deallocate(vEntireDomain)
        endif
        if (allocated(vEntireDomain)) stop 'Error in deallocate_data_buffers'
    end subroutine deallocate_data_buffers

    subroutine compute_velocity_gradient_tensor(uFourier, vFourier, wFourier, &
        omega1Fourier, omega3Fourier, waveVectorX, waveVectorZ, planeIndex, &
        startingPlaneIndex)
        complex(kind=sp), intent(in), dimension(0:mx1,0:mz1) :: uFourier, &
            vFourier, wFourier, omega1Fourier, omega3Fourier
        complex(kind=sp), intent(in), dimension(0:mx1) :: waveVectorX
        complex(kind=sp), intent(in), dimension(0:mz1) :: waveVectorZ
        integer, intent(in) :: planeIndex, startingPlaneIndex
        complex(kind=sp), dimension(0:mx1, 0:mz1) :: dudxFourier, dudyFourier, &
            dudzFourier, dvdxFourier, dvdyFourier, dvdzFourier, dwdxFourier, &
            dwdyFourier, dwdzFourier
        real(kind=sp), dimension(mgalx+2, mgalz) :: dudxPhysical, dudyPhysical, &
            dudzPhysical, dvdxPhysical, dvdyPhysical, dvdzPhysical, &
            dwdxPhysical, dwdyPhysical, dwdzPhysical
        integer :: i, k, arrayIndex
        ! spectral derivatives in x and z
        do k = 0, mz1
            dudxFourier(:,k) = waveVectorX * uFourier(:,k)
            dvdxFourier(:,k) = waveVectorX * vFourier(:,k)
            dwdxFourier(:,k) = waveVectorX * wFourier(:,k)
        enddo
        do i = 0, mx1
            dudzFourier(i,:) = waveVectorZ * uFourier(i,:)
            dvdzFourier(i,:) = waveVectorZ * vFourier(i,:)
            dwdzFourier(i,:) = waveVectorZ * wFourier(i,:)
        enddo
        ! dvdy from continuity
        dvdyFourier = -1.0_sp * (dudxFourier + dwdzFourier)
        ! dudy and dwdy from the vorticity
        dwdyFourier = omega1Fourier + dvdzFourier
        dudyFourier = dvdxFourier - omega3Fourier
        ! transform to physical space
        call fourxz(dudxFourier, dudxPhysical, 1, 1)
        call fourxz(dudyFourier, dudyPhysical, 1, 1)
        call fourxz(dudzFourier, dudzPhysical, 1, 1)
        call fourxz(dvdxFourier, dvdxPhysical, 1, 1)
        call fourxz(dvdyFourier, dvdyPhysical, 1, 1)
        call fourxz(dvdzFourier, dvdzPhysical, 1, 1)
        call fourxz(dwdxFourier, dwdxPhysical, 1, 1)
        call fourxz(dwdyFourier, dwdyPhysical, 1, 1)
        call fourxz(dwdzFourier, dwdzPhysical, 1, 1)
        ! save results
        arrayIndex = planeIndex - startingPlaneIndex + 1
        dudx(:,:,arrayIndex) = dudxPhysical(1:mgalx,:)
        dudy(:,:,arrayIndex) = dudyPhysical(1:mgalx,:)
        dudz(:,:,arrayIndex) = dudzPhysical(1:mgalx,:)
        dvdx(:,:,arrayIndex) = dvdxPhysical(1:mgalx,:)
        dvdy(:,:,arrayIndex) = dvdyPhysical(1:mgalx,:)
        dvdz(:,:,arrayIndex) = dvdzPhysical(1:mgalx,:)
        dwdx(:,:,arrayIndex) = dwdxPhysical(1:mgalx,:)
        dwdy(:,:,arrayIndex) = dwdyPhysical(1:mgalx,:)
        dwdz(:,:,arrayIndex) = dwdzPhysical(1:mgalx,:)
    end subroutine compute_velocity_gradient_tensor

    subroutine save_velocity_fields(planeIndex, startingPlaneIndex, vPlane)
        integer, intent(in) :: planeIndex, startingPlaneIndex
        real(kind=sp), intent(in) :: vPlane(:,:)
        integer :: arrayIndex, upperBound
        integer, parameter :: fftOffset = 2  ! fft requires 2 additional elements in x direction
        upperBound = ubound(vPlane, 1) - fftOffset
        arrayIndex = planeIndex - startingPlaneIndex + 1
        vSlice(:, :, arrayIndex) = vPlane(1:upperBound, :)
    end subroutine

    subroutine write_velocity_gradient_tensor(mpiId, runName, fileNr, &
        fundamentalWavenumberX, fundamentalWavenumberZ, yCoordinateVector, &
        startingPlaneIndices, endingPlaneIndices)
        use hdf5
        integer, intent(in) :: mpiId, fileNr
        character(len=*), intent(in) :: runName
        real(kind=sp), intent(in) :: fundamentalWavenumberX, &
            fundamentalWavenumberZ
        real(kind=sp), intent(in) :: yCoordinateVector(my)
        integer, intent(in) :: startingPlaneIndices(0:numerop-1), &
            endingPlaneIndices(0:numerop-1)
        character(len=3) :: fileNrString
        character(len=256) :: fileName
        character(len=*), parameter :: groupName = "velocity_gradient_tensor"
        integer, parameter :: mpiMaster = 0
        integer(kind=hid_t) :: fileId, groupId
        integer :: hdfError, processId
        call h5open_f(hdfError)  ! initialize hdf interface
        if(mpiId /= mpiMaster) then  ! slaves send data to master
            call send_velocity_gradient_tensor_data(mpiId)
        else  ! master
            ! create file
            write(fileNrString, '(i3.3)') fileNr
            fileName = trim(adjustl(runName)) // '_du.' // fileNrString // '.h5'
            call h5fcreate_f(trim(adjustl(fileName)), H5F_ACC_TRUNC_F, &
                             fileId, hdfError)
            ! write grid (required for visualization)
            call write_grid(fileId, fundamentalWavenumberX, &
                fundamentalWavenumberZ, yCoordinateVector)
            ! create group (aka directory) for velocity gradient tensor
            call h5gcreate_f(fileId, groupName, groupId, hdfError)
            call write_velocity_gradient_tensor_master(groupId, &
                startingPlaneIndices(mpiMaster), endingPlaneIndices(mpiMaster))
            do processId = 1, numerop-1
                call receive_velocity_gradient_tensor_data(processId, &
                    startingPlaneIndices(processId), endingPlaneIndices(processId))
                call write_velocity_gradient_tensor_slave(groupId, &
                    startingPlaneIndices(processId), endingPlaneIndices(processId))
            enddo
            ! close group and file
            call h5gclose_f(groupId, hdfError)
            call h5fclose_f(fileId, hdfError)
        endif
        call h5close_f(hdfError)  ! close hdf interface
    end subroutine write_velocity_gradient_tensor

    subroutine send_velocity_gradient_tensor_data(mpiId)
        include 'mpif.h'
        integer, intent(in) :: mpiId
        integer, parameter :: mpiMaster = 0
        integer :: errorStatus
        call MPI_SEND(dudx, size(dudx), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dudy, size(dudy), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dudz, size(dudz), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dvdx, size(dvdx), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dvdy, size(dvdy), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dvdz, size(dvdz), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dwdx, size(dwdx), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dwdy, size(dwdy), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
        call MPI_SEND(dwdz, size(dwdz), MPI_REAL, mpiMaster, mpiId, &
                      MPI_COMM_WORLD, errorStatus)
    end subroutine send_velocity_gradient_tensor_data

    subroutine receive_velocity_gradient_tensor_data(mpiId, &
            startingPlaneIndex, endingPlaneIndex)
        include 'mpif.h'
        integer, intent(in) :: mpiId, startingPlaneIndex, endingPlaneIndex
        integer :: numberOfPlanes, errorStatus, mpiStatus(MPI_STATUS_SIZE)
        integer, parameter :: mpiMaster = 0
        numberOfPlanes = endingPlaneIndex - startingPlaneIndex + 1
        ! reallocate data buffers in case the size of the domain slice
        ! of the current MPI process is different
        ! the data buffers (dudx, etc.) all have the same dimensions, we only
        ! need to compare the size of one (here we arbitrarily choose dudx)
        if(size(dudx, 3) /= numberOfPlanes) then
            call deallocate_data_buffers(mpiMaster)
            call allocate_data_buffers(startingPlaneIndex, endingPlaneIndex, mpiMaster)
        endif
        ! receive data from slave
        call MPI_RECV(dudx, size(dudx), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dudy, size(dudy), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dudz, size(dudz), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dvdx, size(dvdx), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dvdy, size(dvdy), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dvdz, size(dvdz), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dwdx, size(dwdx), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dwdy, size(dwdy), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
        call MPI_RECV(dwdz, size(dwdz), MPI_REAL, mpiId, mpiId, &
                      MPI_COMM_WORLD, mpiStatus, errorStatus)
    end subroutine receive_velocity_gradient_tensor_data

    subroutine write_grid(fileId, fundamentalWavenumberX, &
        fundamentalWavenumberZ, yCoordinateVector)
        use hdf5
        integer(kind=hid_t), intent(in) :: fileId
        real(kind=sp), intent(in) :: fundamentalWavenumberX, &
            fundamentalWavenumberZ
        real(kind=sp), intent(in) :: yCoordinateVector(:)
        real(kind=sp), parameter :: pi = acos(-1.0_sp)
        real(kind=sp) :: Lx, Lz, dx, dz
        integer :: i, hdfError
        real(kind=sp) :: xCoordinateVector(mgalx), zCoordinateVector(mgalz)
        character(len=*), parameter :: groupName = "grid"
        integer(kind=hid_t) :: groupId, dataSpaceId, dataSetId
        integer, parameter :: rank = 1
        integer(kind=hsize_t) :: dataDimensions(rank)
        ! generate coordinate vectors in x and z
        Lx = 2.0_sp * pi / fundamentalWavenumberX
        dx = Lx / mgalx
        do i = 1, mgalx
            xCoordinateVector(i) = (i-1) * dx
        enddo
        Lz = 2.0_sp * pi / fundamentalWavenumberZ
        dz = Lz / mgalz
        do i = 1, mgalz
            zCoordinateVector(i) = (i-1) * dz
        enddo
        ! write to file
        ! create group (aka folder) for coordinate vectors
        call h5gcreate_f(fileId, groupName, groupId, hdfError)
        ! write x coordinate vector
        dataDimensions(1) = mgalx
        call h5screate_simple_f(rank, dataDimensions, dataSpaceId, hdfError)
        call h5dcreate_f(groupId, "x", H5T_IEEE_F32LE, dataSpaceId, &
                         dataSetId, hdfError)
        call h5dwrite_f(dataSetId, H5T_NATIVE_REAL, xCoordinateVector, &
                        dataDimensions, hdfError)
        call h5dclose_f(dataSetId, hdfError)
        call h5sclose_f(dataSpaceId, hdfError)
        ! write y coordinate vector
        dataDimensions(1) = my
        call h5screate_simple_f(rank, dataDimensions, dataSpaceId, hdfError)
        call h5dcreate_f(groupId, "y", H5T_IEEE_F32LE, dataSpaceId, &
                         dataSetId, hdfError)
        call h5dwrite_f(dataSetId, H5T_NATIVE_REAL, yCoordinateVector, &
                        dataDimensions, hdfError)
        call h5dclose_f(dataSetId, hdfError)
        call h5sclose_f(dataSpaceId, hdfError)
        ! write z coordinate vector
        dataDimensions(1) = mgalz
        call h5screate_simple_f(rank, dataDimensions, dataSpaceId, hdfError)
        call h5dcreate_f(groupId, "z", H5T_IEEE_F32LE, dataSpaceId, &
                         dataSetId, hdfError)
        call h5dwrite_f(dataSetId, H5T_NATIVE_REAL, zCoordinateVector, &
                        dataDimensions, hdfError)
        call h5dclose_f(dataSetId, hdfError)
        call h5sclose_f(dataSpaceId, hdfError)
        ! close group
        call h5gclose_f(groupId, hdfError)
    end subroutine write_grid

    subroutine write_velocity_gradient_tensor_master(groupId, &
        startingPlaneIndex, endingPlaneIndex)
        use hdf5
        integer(kind=hid_t), intent(in) :: groupId
        integer, intent(in) :: startingPlaneIndex, endingPlaneIndex
        call write_tensor_component_master(groupId, "dudx", dudx, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dudy", dudy, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dudz", dudz, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dvdx", dvdx, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dvdy", dvdy, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dvdz", dvdz, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dwdx", dwdx, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dwdy", dwdy, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_master(groupId, "dwdz", dwdz, &
            startingPlaneIndex, endingPlaneIndex)
    end subroutine write_velocity_gradient_tensor_master

    subroutine write_velocity_gradient_tensor_slave(groupId, &
        startingPlaneIndex, endingPlaneIndex)
        use hdf5
        integer(kind=hid_t), intent(in) :: groupId
        integer, intent(in) :: startingPlaneIndex, endingPlaneIndex
        call write_tensor_component_slave(groupId, "dudx", dudx, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dudy", dudy, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dudz", dudz, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dvdx", dvdx, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dvdy", dvdy, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dvdz", dvdz, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dwdx", dwdx, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dwdy", dwdy, &
            startingPlaneIndex, endingPlaneIndex)
        call write_tensor_component_slave(groupId, "dwdz", dwdz, &
            startingPlaneIndex, endingPlaneIndex)
    end subroutine write_velocity_gradient_tensor_slave

    subroutine write_tensor_component_master(locationId, dataSetName, &
        dataBuffer, startingPlaneIndex, endingPlaneIndex)
        use hdf5
        integer(kind=hid_t), intent(in) :: locationId
        character(len=*), intent(in) :: dataSetName
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        integer, intent(in) :: startingPlaneIndex, endingPlaneIndex
        integer, parameter :: rank = 3
        integer(kind=hsize_t), dimension(rank) :: dimDataArray, &
            dimMemorySpace, offset, numberOfBlocks
        integer(kind=hid_t) :: dataSpaceId, dataSetId, memorySpaceId
        integer :: hdfError
        ! disk data corresponds to computational grid
        dimDataArray(1) = mgalx
        dimDataArray(2) = mgalz
        dimDataArray(3) = my
        ! memory data corresponds to process slice
        dimMemorySpace(1) = mgalx
        dimMemorySpace(2) = mgalz
        dimMemorySpace(3) = endingPlaneIndex - startingPlaneIndex + 1
        ! create disk dataspace (data array)
        call h5screate_simple_f(rank, dimDataArray, dataSpaceId, hdfError)
        ! create memory dataspace
        call h5screate_simple_f(rank, dimMemorySpace, memorySpaceId, hdfError)
        ! create disk dataset
        call h5dcreate_f(locationId, dataSetName, H5T_IEEE_F32LE, dataSpaceId, &
                         dataSetId, hdfError)
        ! select subset of disk dataset
        offset(1) = 0
        offset(2) = 0
        offset(3) = startingPlaneIndex - 1
        numberOfBlocks = dimMemorySpace
        call h5sselect_hyperslab_f(dataSpaceId, H5S_SELECT_SET_F, offset, &
                                   numberOfBlocks, hdfError)
        ! write data
        call h5dwrite_f(dataSetId, H5T_NATIVE_REAL, dataBuffer, dimMemorySpace, &
                        hdfError, memorySpaceId, dataSpaceId)
        ! close interfaces
        call h5dclose_f(dataSetId, hdfError)
        call h5sclose_f(memorySpaceId, hdfError)
        call h5sclose_f(dataSpaceId, hdfError)
    end subroutine write_tensor_component_master

    subroutine write_tensor_component_slave(locationId, dataSetName, &
        dataBuffer, startingPlaneIndex, endingPlaneIndex)
        use hdf5
        integer(kind=hid_t), intent(in) :: locationId
        character(len=*), intent(in) :: dataSetName
        real(kind=sp), intent(in) :: dataBuffer(:,:,:)
        integer, intent(in) :: startingPlaneIndex, endingPlaneIndex
        integer, parameter :: rank = 3
        integer(kind=hsize_t), dimension(rank) :: dimMemorySpace, &
            offset, numberOfBlocks
        integer(kind=hid_t) :: dataSetId, dataSpaceId, memorySpaceId
        integer :: hdfError
        ! memory data corresponds to process slice
        dimMemorySpace(1) = mgalx
        dimMemorySpace(2) = mgalz
        dimMemorySpace(3) = endingPlaneIndex - startingPlaneIndex + 1
        ! create memory dataspace
        call h5screate_simple_f(rank, dimMemorySpace, memorySpaceId, hdfError)
        ! open existing disk dataset
        call h5dopen_f(locationId, dataSetName, dataSetId, hdfError)
        ! get disk dataspace identifier
        call h5dget_space_f(dataSetId, dataSpaceId, hdfError)
        ! select subset of disk dataset
        offset(1) = 0
        offset(2) = 0
        offset(3) = startingPlaneIndex - 1
        numberOfBlocks = dimMemorySpace
        call h5sselect_hyperslab_f(dataSpaceId, H5S_SELECT_SET_F, offset, &
                                   numberOfBlocks, hdfError)
        ! write data
        call h5dwrite_f(dataSetId, H5T_NATIVE_REAL, dataBuffer, dimMemorySpace, &
                      hdfError, memorySpaceId, dataSpaceId)
        ! close interfaces
        call h5dclose_f(dataSetId, hdfError)
        call h5sclose_f(dataSpaceId, hdfError)
        call h5sclose_f(memorySpaceId, hdfError)
    end subroutine write_tensor_component_slave

    subroutine write_velocity_fields(mpiRank, startingPlaneIndices, endingPlaneIndices, runName, &
        fileNr, fundamentalWavenumberX, fundamentalWavenumberZ, yCoordinateVector, time)
        integer, intent(in) :: mpiRank
        integer, intent(in) :: startingPlaneIndices(0:), endingPlaneIndices(0:)
        character(len=*), intent(in) :: runName
        integer, intent(in) :: fileNr
        real(kind=sp), intent(in) :: fundamentalWavenumberX, fundamentalWavenumberZ
        real(kind=sp), intent(in) :: yCoordinateVector(:), time
        integer, parameter :: mpiMaster = 0
        call send_receive_velocity_fields(mpiRank, startingPlaneIndices, endingPlaneIndices)
        if (mpiRank == mpiMaster) then
            call write_hdf_file_velocity_fields(runName, fileNr, fundamentalWavenumberX, &
                fundamentalWavenumberZ, yCoordinateVector)
            call write_xml_file_velocity_fields(runName, fileNr, time, mgalx, my, mgalz)
        endif
    end subroutine

    subroutine send_receive_velocity_fields(mpiRank, startingPlaneIndices, endingPlaneIndices)
        use mpi
        integer, intent(in) :: mpiRank
        integer, intent(in) :: startingPlaneIndices(0:), endingPlaneIndices(0:)
        integer, parameter :: mpiMaster = 0
        integer :: errorStatus, idxStart, idxEnd, numberOfPlanes, processId
        if (mpiRank /= mpiMaster) then  ! slaves send data to master
            call MPI_SEND(vSlice, size(vSlice), MPI_REAL, mpiMaster, mpiRank, &
                          MPI_COMM_WORLD, errorStatus)
        else  ! master
            ! stores own data
            idxStart = startingPlaneIndices(mpiMaster)
            idxEnd = endingPlaneIndices(mpiMaster)
            vEntireDomain(:, :, idxStart:idxEnd) = vSlice(:, :, :)
            ! and receives data from slaves and stores it
            do processId = 1, numerop-1
                idxStart = startingPlaneIndices(processId)
                idxEnd = endingPlaneIndices(processId)
                numberOfPlanes = idxEnd - idxStart + 1
                ! reallocate data buffer in case the size of the domain slice of
                ! the current MPI process is different
                if (size(vSlice, 3) /= numberOfPlanes) then
                    deallocate(vSlice)
                    allocate(vSlice(mgalx, mgalz, numberOfPlanes))
                endif
                call MPI_RECV(vSlice, size(vSlice), MPI_REAL, processId, processId, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, errorStatus)
                vEntireDomain(:, :, idxStart:idxEnd) = vSlice(:, :, :)
            enddo
        endif
    end subroutine

    subroutine write_hdf_file_velocity_fields(runName, fileNr, fundamentalWavenumberX, &
        fundamentalWavenumberZ, yCoordinateVector)
        use hdf5
        character(len=*), intent(in) :: runName
        integer, intent(in) :: fileNr
        real(kind=sp), intent(in) :: fundamentalWavenumberX, fundamentalWavenumberZ
        real(kind=sp), intent(in) :: yCoordinateVector(:)
        character(:), allocatable :: fileName
        character(len=3) :: fileNrString
        real(kind=sp) :: pi, domainSizeX, domainSizeZ, deltaX, deltaZ
        integer :: nGridpointsX, nGridpointsZ, i
        real(kind=sp) :: xCoordinateVector(1:ubound(vEntireDomain, 1)), &
                         zCoordinateVector(1:ubound(vEntireDomain, 2))
        integer :: error, hdfError
        integer(kind=hid_t) :: fileId, fspaceId, dsetId
        integer(kind=HSIZE_T), dimension(1) :: dim1d
        integer(kind=HSIZE_T), dimension(3) :: dim3d

        ! set filename
        write(fileNrString, '(i3.3)') fileNr
        fileName = trim(adjustl(runName)) // '_vel.' // fileNrString // '.h5'

        ! generate grid coordinates
        pi = acos(-1.0_sp)  ! as in cfdiff8.v7.f
        nGridpointsX = ubound(vEntireDomain, 1)
        domainSizeX = 2.0_sp * pi / fundamentalWavenumberX
        deltaX = domainSizeX / real(nGridpointsX)
        do i = 1, nGridpointsX
            xCoordinateVector(i) = (i-1) * deltaX
        enddo
        nGridpointsZ = ubound(vEntireDomain, 2)
        domainSizeZ = 2.0_sp * pi / fundamentalWavenumberZ
        deltaZ = domainSizeZ / real(nGridpointsZ)
        do i = 1, nGridpointsZ
            zCoordinateVector(i) = (i-1) * deltaZ
        enddo

        ! write hdf file
        call h5open_f(error) ! initialize interface
        call h5fcreate_f(fileName, H5F_ACC_TRUNC_F, fileId, hdfError) ! create hdf5 file
        ! write x coordinate
        dim1d(1) = nGridpointsX
        call h5screate_simple_f(1, dim1d, fspaceId, hdfError) ! create dataspace
        call h5dcreate_f(fileId, 'x', H5T_NATIVE_REAL, fspaceId, dsetId, hdfError) ! create dataset
        call h5dwrite_f(dsetId, H5T_NATIVE_REAL, xCoordinateVector, dim1d, hdfError) ! write data
        call h5dclose_f(dsetId, hdfError) ! close dataset
        call h5sclose_f (fspaceId, hdfError) ! close dataspace
        ! write y coordinate
        dim1d(1) = size(yCoordinateVector)
        call h5screate_simple_f(1, dim1d, fspaceId, hdfError)
        call h5dcreate_f(fileId, 'y', H5T_NATIVE_REAL, fspaceId, dsetId, hdfError)
        call h5dwrite_f(dsetId, H5T_NATIVE_REAL, yCoordinateVector, dim1d, hdfError)
        call h5dclose_f(dsetId, hdfError)
        call h5sclose_f (fspaceId, hdfError)     
        ! write z coordinate
        dim1d(1) = nGridpointsZ
        call h5screate_simple_f(1, dim1d, fspaceId, hdfError)
        call h5dcreate_f(fileId, 'z', H5T_NATIVE_REAL, fspaceId, dsetId, hdfError)
        call h5dwrite_f(dsetId, H5T_NATIVE_REAL, zCoordinateVector, dim1d, hdfError)
        call h5dclose_f(dsetId, hdfError)
        call h5sclose_f (fspaceId, hdfError)
        ! write wall normal velocity
        dim3d(1) = nGridpointsX
        dim3d(2) = nGridpointsZ
        dim3d(3) = size(yCoordinateVector)
        call h5screate_simple_f(3, dim3d, fspaceId, hdfError)
        call h5dcreate_f(fileId, 'v', H5T_NATIVE_REAL, fspaceId, dsetId, hdfError)
        call h5dwrite_f(dsetId, H5T_NATIVE_REAL, vEntireDomain, dim3d, hdfError)
        call h5dclose_f(dsetId, hdfError)
        call h5sclose_f(fspaceId, hdfError)
        ! close interfaces
        call h5fclose_f(fileId, hdfError) ! close hdf5 file
        call h5close_f(error) ! close interface
    end subroutine

    subroutine write_xml_file_velocity_fields(runName, fileNr, time, nPointsX, nPointsY, &
        nPointsZ)
        character(len=*), intent(in) :: runName
        integer, intent(in) :: fileNr
        real(kind=sp), intent(in) :: time
        integer, intent(in) :: nPointsX, nPointsY, nPointsZ
        character(:), allocatable :: xmlFile, hdfFileNoPath
        character(len=3) :: fileNrString
        character(len=512) :: buffer
        integer :: fileId

        ! set filename
        write(fileNrString, '(i3.3)') fileNr
        xmlFile = trim(adjustl(runName)) // '_vel.' // fileNrString // '.xmf'
        hdfFileNoPath = trim(adjustl(runName(index(runName, '/', .true.)+1 :))) // '_vel.' // fileNrString // '.h5'

        ! write xmf file
        open(newunit=fileId, file=xmlFile, status='new')
        ! file header
        buffer = '<?xml version="1.0" ?>'
        write(fileId,'(a)') trim(buffer)
        buffer = '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
        write(fileId,'(a)') trim(buffer)
        buffer = '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">'
        write(fileId,'(a)') trim(buffer)
        buffer = ' <Domain>'
        write(fileId,'(a)') trim(buffer)
        ! grid information
        buffer = '   <Grid GridType="Uniform">'
        write(fileId,'(a)') trim(buffer)
        write(buffer,'(a,f15.5,a)') '     <Time Value="',time,'"/>'
        write(fileId,'(a)') trim(buffer)
        write(buffer,'(a,3(i8),a)') &
            '     <Topology TopologyType="3DRectMesh" Dimensions="',nPointsY,nPointsZ,nPointsX,'"/>'
        write(fileId,'(a)') trim(buffer)
        ! coordinates vectors
        buffer = '     <Geometry GeometryType="VXVYVZ">'
        write(fileId,'(a)') trim(buffer)
        ! x
        write(buffer,'(a,i8,a)') &
            '       <DataItem Name="X" Dimensions="',nPointsX,'" NumberType="Float" Precision="4" Format="HDF">'
        write(fileId,'(a)') trim(buffer)
        write(buffer,'(3(a))') '       ',trim(hdfFileNoPath),':/x'
        write(fileId,'(a)') trim(buffer)
        buffer = '       </DataItem>'
        write(fileId,'(a)') trim(buffer)
        ! z
        write(buffer,'(a,i8,a)') &
            '       <DataItem Name="Z" Dimensions="',nPointsZ,'" NumberType="Float" Precision="4" Format="HDF">'
        write(fileId,'(a)') trim(buffer)
        write(buffer,'(3(a))') '       ',trim(hdfFileNoPath),':/z'
        write(fileId,'(a)') trim(buffer)
        buffer = '       </DataItem>'
        write(fileId,'(a)') trim(buffer)
        ! y
        write(buffer,'(a,i8,a)') &
            '       <DataItem Name="Y" Dimensions="',nPointsY,'" NumberType="Float" Precision="4" Format="HDF">'
        write(fileId,'(a)') trim(buffer)
        write(buffer,'(3(a))') '       ',trim(hdfFileNoPath),':/y'
        write(fileId,'(a)') trim(buffer)
        buffer = '       </DataItem>'
        write(fileId,'(a)') trim(buffer)
        buffer = '     </Geometry>'
        write(fileId,'(a)') trim(buffer)
        ! data
        buffer = '     <Attribute Name="v" Active="1" AttributeType="Scalar" Center="Node">'
        write(fileId,'(a)') trim(buffer)
        write(buffer,'(a,3(i8),a)') &
            '       <DataItem Dimensions="',nPointsY,nPointsZ,nPointsX,'" NumberType="Float" Precision="4" Format="HDF">'
        write(fileId,'(a)') trim(buffer)
        write(buffer,'(3(a))') '       ', trim(hdfFileNoPath),':/v'
        write(fileId,'(a)') trim(buffer)
        buffer = '       </DataItem>'
        write(fileId,'(a)') trim(buffer)
        buffer = '     </Attribute>'
        ! close xmf file header
        write(fileId,'(a)') trim(buffer)
        buffer = '   </Grid>'
        write(fileId,'(a)') trim(buffer)
        buffer = ' </Domain>'
        write(fileId,'(a)') trim(buffer)
        buffer = '</Xdmf>'
        write(fileId,'(a)') trim(buffer)
        ! close xmf file
        close(fileId)
    end subroutine

end module velocity_gradient_tensor
