module h5save
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private

    public :: check_filename
    public :: h5save_R_dp, h5save_R1_dp, h5save_R2_dp
    public :: h5save_CPartial_Init
    public :: h5save_C3Partial_SingleDim3_sp, h5save_C3Partial_SingleDim2_sp
    public :: h5save_C4Partial_SingleDim3_dp

    ! The following standard for complex variables are used:
    ! if varibale 'var' is complex, save as '/var/var_REAL' and '/var/var_IMAG'
    ! (one group with two datasets)
    ! Same as my matlab h5 libaries

    ! NOTE on precision:
    ! This code is tailor made for the DNS, all scalars, such as time and Re
    ! are saved in double precision, flow fields are saved in single precision
    !
    ! A more general version of this code that also includes more functions can
    ! be found in the DNSFortranLibraries

contains
    ! subroutine check_filename( filename, myid )
    ! This function check if the file with filename already exist
    ! if it exist, will append the current date and time to the filename to
    ! to create a new filename and send to all processors
    ! Arguments:
    !   filename: [string, Input/Output] h5 filename with path, it is appended
    !                                    with date and time if file already exist
    !   myid    : [integer, Input] processor ID
    subroutine check_filename( filename, myid )
        character(len=*), intent(inout) :: filename
        integer, intent(in) :: myid

        logical :: file_exists

        integer :: values(8)
        character(len=10) :: current_data_and_time

        integer :: loc, count

        integer :: ierr


        if ( myid .eq. 0 ) then
            ! master checks if file exist
            INQUIRE(FILE=filename, EXIST=file_exists)

            if ( file_exists ) then
                ! get current date and time
                call date_and_time(VALUES=values)
                write(current_data_and_time,"(5I2.2)") values(2), values(3), values(5), values(6), values(7)

                ! append date and time to current filename if file already exist
                loc = INDEX( filename, ".h5" )
                filename = filename(1:loc-1) // '_' // current_data_and_time // '.h5'
            endif
        endif

        ! send filename to all processors regardless of whether it is appended
        ! with the date and time or not
        count = LEN( filename )
        call MPI_BCAST(filename, count, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        ! write to output file
        if ( myid .eq. 0 ) then
            write(*,*) "  Saving data to file: ", filename
        endif

    end subroutine check_filename


    ! subroutine h5save_R_dp( filename, varname, scalar )
    ! save a real number to h5 file
    !
    ! DOUBLE PRECISION ONLY!
    !
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name saved in h5
    !   scalar  : [double/single scalar, Input] data to be saved
    subroutine h5save_R_dp( filename, varname, scalar )
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), intent(in) :: scalar

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(1) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists


        data_dim = 1

        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        else
            ! Create new file
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error)
        endif

        dset_name = varname
        ! Create dataspace with rank 0 and size 1
        CALL h5screate_simple_f(0, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var' and write data
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, scalar, data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_R_dp


    ! subroutine h5save_R1_dp( filename, varname, vector )
    ! save a real rank 1 vector to h5 file
    !
    ! DOUBLE PRECISION ONLY!
    !
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (real 1d vector)
    !   vector  : [double/single 1d vector, Input] data to be saved
    subroutine h5save_R1_dp( filename, varname, vector )
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), intent(in), dimension(:) :: vector

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(1) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists


        ! get vector dimensions
        data_dim = shape(vector)

        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        else
            ! Create new file
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error)
        endif

        dset_name = varname
        ! Create dataspace with rank 1 and size data_dim
        CALL h5screate_simple_f(1, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var' and write data
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, vector, data_dim, error)

        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_R1_dp


    ! subroutine h5save_R2_dp( filename, varname, matrix )
    ! save a real rank 2 matrix to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (real 2d matrix)
    !   matrix  : [double 2d matrix, Input] data to be saved
    subroutine h5save_R2_dp( filename, varname, matrix )
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), intent(in), dimension(:,:) :: matrix

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(2) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists


        ! get matrix dimensions
        data_dim = shape(matrix)

        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        else
            ! Create new file
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error)
        endif

        dset_name = varname
        ! Create dataspace with rank 2 and size data_dim
        CALL h5screate_simple_f(2, data_dim, dspace_id, error)
        ! Create double precision dataset with path '/var' and write data
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, matrix, data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_R2_dp


    ! subroutine h5save_CPartial_Init( filename, varname, full_data_dim, precision )
    ! initialize a complex rank N matrix to h5 file in preparation for partial data saving
    !
    ! SINGLE PRECISION ONLY!
    !
    ! Arguments:
    !   filename     : [string, Input] h5 filename with path
    !   varname      : [string, Input] variable name
    !   full_data_dim: [integer, size N, Input] the dimension of the full data matrix
    !   precision    : [integer, sp of dp]
    !
    ! Note:
    !   This function needs to be called once and only once before any partial saving is performed
    subroutine h5save_CPartial_Init( filename, varname, full_data_dim, precision)
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in), dimension(:) :: full_data_dim
        integer, intent(in) :: precision

        integer(HSIZE_T), dimension(:), allocatable :: full_data_dim2
        character(len=100) :: dset_name ! dataset name
        integer :: matrix_rank

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: group_id  ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id   ! dataset id

        logical :: file_exists


        ! Get matrix rank
        matrix_rank = size(full_data_dim)
        ! allocate full_data_dim2
        ALLOCATE(full_data_dim2(matrix_rank))
        ! Integer type conversion for the full data dimension
        full_data_dim2 = INT( full_data_dim, HSIZE_T)

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        else
            ! Create new file
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error)
        endif

        ! Create a group in the HDF5 file with the variable name
        CALL h5gcreate_f(file_id, varname, group_id, error)
        ! Close the group
        CALL h5gclose_f(group_id, error)

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! Create disk dataspace with rank N and size full_data_dim
        call h5screate_simple_f(matrix_rank, full_data_dim2, dspace_id, error)

        ! Create double/single precision dataset with path '/var/var_REAL'
        if (precision .eq. sp) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
        else if (precision .eq. dp) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        endif
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"

        ! Create disk dataspace with rank N and size full_data_dim
        call h5screate_simple_f(matrix_rank, full_data_dim2, dspace_id, error)

        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        if (precision .eq. sp) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
        else if (precision .eq. dp) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        endif
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
        ! Deallocate memory
        DEALLOCATE(full_data_dim2)

    end subroutine h5save_CPartial_Init


    ! subroutine h5save_C3Partial_SingleDim3_sp( filename, varname, matrix, dim3index )
    ! save a complex rank 2 matrix to h5 file as a single plane in the complex rank 3 matrix
    !
    ! SINGLE PRECISION ONLY!
    !
    ! Arguments:
    !   filename : [string, Input] h5 filename with path
    !   varname  : [string, Input] variable name
    !   matrix   : [double/single complex 2d matrix, Input] data to be saved
    !   dim3index: [integer, Input] the index for the third dimension
    !
    ! Note:
    !   The variable in this h5 file must already be initilaized using the function h5save_C3Partial_Init
    subroutine h5save_C3Partial_SingleDim3_sp( filename, varname, matrix, dim3index)
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        complex(kind=sp), intent(in), dimension(:,:) :: matrix
        integer, intent(in) :: dim3index

        real(kind=sp), dimension(:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(2) :: matrix_dim
        integer(HSIZE_T), dimension(3) :: slice_data_dim
        integer(HSIZE_T), dimension(3) :: slabOffset ! hyperslab offset

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id


        ! get matrix dimension
        matrix_dim = shape(matrix)
        slice_data_dim(1) = matrix_dim(1)
        slice_data_dim(2) = matrix_dim(2)
        slice_data_dim(3) = 1
        ! Allocate temp buffer
        allocate(temp( slice_data_dim(1), slice_data_dim(2)))

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Open file with read and write access
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        ! Group is already created by h5save_C3Partial_Init

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim3index (offset dim3index-1)
        slabOffset = (/ 0, 0, dim3index - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, sp)
        call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim3index (offset dim3index-1)
        slabOffset = (/ 0, 0, dim3index - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp = aimag(matrix)
        call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
        ! Deallocate temp buffer
        DEALLOCATE(temp)
    end subroutine h5save_C3Partial_SingleDim3_sp


    ! subroutine h5save_C3Partial_SingleDim2_sp( filename, varname, matrix, dim2index )
    ! save a complex rank 2 matrix to h5 file as a single plane of dimension 2
    ! of the full complex rank 3 matrix
    !
    ! SINGLE PRECISION ONLY!
    !
    ! Arguments:
    !   filename : [string, Input] h5 filename with path
    !   varname  : [string, Input] variable name
    !   matrix   : [single complex 2d matrix, size (dim1size, dim3size) Input] data to be saved
    !   dim2index: [integer, Input] the index for the second dimension
    !
    ! Note:
    !   The variable in this h5 file must already be initilaized using the function h5save_CPartial_Init
    subroutine h5save_C3Partial_SingleDim2_sp( filename, varname, matrix, dim2index)
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        complex(kind=sp), intent(in), dimension(:,:) :: matrix
        integer, intent(in) :: dim2index

        real(kind=sp), dimension(:,:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(2) :: matrix_dim
        integer(HSIZE_T), dimension(3) :: slice_data_dim
        integer(HSIZE_T), dimension(3) :: slabOffset ! hyperslab offset

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id


        ! get matrix dimension
        matrix_dim = shape(matrix)
        slice_data_dim(1) = matrix_dim(1)
        slice_data_dim(2) = 1
        slice_data_dim(3) = matrix_dim(2)
        ! Allocate temp buffer
        allocate(temp( slice_data_dim(1), slice_data_dim(2), slice_data_dim(3)))

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Open file with read and write access
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        ! Group is already created by h5save_CPartial_Init

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim2index (offset dim2index-1)
        slabOffset = (/ 0, dim2index - 1, 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp(:,1,:) = real(matrix, sp)
        call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim2index (offset dim2index-1)
        slabOffset = (/ 0, dim2index - 1, 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp(:,1,:) = aimag(matrix)
        call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
        ! Deallocate temp buffer
        DEALLOCATE(temp)
    end subroutine h5save_C3Partial_SingleDim2_sp


    ! subroutine h5save_C4Partial_SingleDim3_dp( filename, varname, matrix, dim3index )
    ! save a complex rank 3 matrix to h5 file as a single plane in the complex rank 4 matrix
    !
    ! Arguments:
    !   filename : [string, Input] h5 filename with path
    !   varname  : [string, Input] variable name
    !   matrix   : [double complex 3d matrix, Input] data to be saved
    !   dim3index: [integer, Input] the index for the third dimension
    !
    ! Note:
    !   The variable in this h5 file must already be initilaized using the function h5save_CPartial_Init
    !
    ! DOUBLE PRECISION ONLY!
    subroutine h5save_C4Partial_SingleDim3_dp( filename, varname, matrix, dim3index)
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        complex(kind=dp), intent(in), dimension(:,:,:) :: matrix
        integer, intent(in) :: dim3index

        real(kind=dp), dimension(:,:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(3) :: matrix_dim
        integer(HSIZE_T), dimension(4) :: slice_data_dim
        integer(HSIZE_T), dimension(4) :: slabOffset ! hyperslab offset

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id


        ! get matrix dimension
        matrix_dim = shape(matrix)
        slice_data_dim(1) = matrix_dim(1)
        slice_data_dim(2) = matrix_dim(2)
        slice_data_dim(3) = 1
        slice_data_dim(4) = matrix_dim(3)
        ! Allocate temp buffer
        allocate(temp( slice_data_dim(1), slice_data_dim(2), slice_data_dim(4)))

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Open file with read and write access
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        ! Group is already created by h5save_CPartial_Init

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim3index (offset dim3index-1)
        slabOffset = (/ 0, 0, dim3index - 1, 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 4 and size slice_data_dim
        call h5screate_simple_f(4, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, dp)
        call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim3index (offset dim3index-1)
        slabOffset = (/ 0, 0, dim3index - 1 , 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 4 and size slice_data_dim
        call h5screate_simple_f(4, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp = aimag(matrix)
        call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
        ! Deallocate temp buffer
        DEALLOCATE(temp)
    end subroutine h5save_C4Partial_SingleDim3_dp

end module h5save

