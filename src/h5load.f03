module h5load
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private

    public :: h5load_R_dp, h5load_R1_dp, h5load_R2_dp, h5load_C2_sp, h5load_C3_sp
    public :: h5load_C4_partial_dp
    public :: h5load_check_variable_existence


    ! Note:
    ! h5load_R  reads a real scalar
    ! h5load_R1 reads a real rank 1 matrix (a vector)
    !           This will not work if the h5 data size is scalar
    !           It will work if the h5 data size is a 1X1 matrix (technically a scalar), but it will be easier to use h5load_R
    ! h5load_C2 reads a complex rank 2 matrix
    ! h5load_C3 reads a complex rank 3 matrix


contains
    ! function output = h5load_R_dp(filename, varname )
    ! Read real scalar from h5 file
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname : [string, Input]
    !             variable name in h5 file, must be real numerical scalar
    ! Output:
    !   output:   [double/single scalar]
    function h5load_R_dp(filename, varname ) result(scalar)
        character(len=*), intent(in) :: filename, varname
        real(kind=dp) :: scalar

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier

        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims

        INTEGER :: error ! Error flag


        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, varname, dset_id, error)

        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of buffer
        data_dims(1) = 1
        CALL h5dread_f(dset_id, H5T_IEEE_F64LE, scalar, data_dims, error)

        ! close dataset
        CALL h5dclose_f(dset_id, error)
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
    end function h5load_R_dp


    ! function ouput = h5load_R1_dp(filename, varname )
    ! Read real rank 1 matrix (vector) from h5 file
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname : [string, Input]
    !             variable name in h5 file, must be real numerical vector
    ! Output:
    !   output:   [double vector]
    function h5load_R1_dp(filename, varname ) result(vector)
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), dimension(:,:), allocatable :: matrix
        real(kind=dp), dimension(:), allocatable :: vector

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier
        INTEGER(HID_T) :: space_id       ! Dataspace identifier

        INTEGER :: error ! Error flag
        INTEGER :: dim1, dim2 ! matrix dimensions

        INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
        INTEGER(HSIZE_T), DIMENSION(2) :: max_dims


        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, varname, dset_id, error)
        !Get dataspace ID
        CALL h5dget_space_f(dset_id, space_id,error)
        ! Initialize this to 0, otherwise if the data is rank 1, data_dim(2) will
        ! remain the random initial value, which might be problematic
        data_dims(1) = 0
        data_dims(2) = 0
        !Get dataspace dims
        CALL h5sget_simple_extent_dims_f(space_id, data_dims, max_dims, error)
        dim1 = data_dims(1)
        dim2 = data_dims(2)
        CALL h5sclose_f(space_id, error)

        ! Check the dimensions, if it is a scalar, then throw an error
        if ((dim1 .eq. 0) .and. (dim2 .eq. 0)) then
            write(*,*) "Error: reading a scalar with h5load_R1, use h5load_R instead"
            ! close dataset
            CALL h5dclose_f(dset_id, error)
            ! close file
            CALL h5fclose_f(file_id, error)
            ! close h5 interface
            CALL h5close_f(error)
            stop 2
        endif
        ! Check the dimensions, if it is not a vector or column vector, then throw an error
        if ( dim2 .gt. 1) then
            write(*,*) "Error: variable is not a vector or a column vector. h5load_R1 failed."
            ! close dataset
            CALL h5dclose_f(dset_id, error)
            ! close file
            CALL h5fclose_f(file_id, error)
            ! close h5 interface
            CALL h5close_f(error)
            stop 2
        endif

        ! Allocate dimensions to dset_data for reading
        if (dim2 .eq. 0) then
            dim2 = 1
            ! if data is rank 1, dim2 is 0 (initial value), we will set it to 1 for
            ! ALLOCATE(matrix(dim1,dim2))
            ! we keep the this to rank 2 so that it can properly read NX1 vectors
        endif
        ALLOCATE(matrix(dim1,dim2))
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F64LE, matrix, data_dims, error)

        ! Convert to dp and reshape matrix into a vector
        vector = reshape( real(matrix, dp), (/ dim1 /) )

        ! Deallocate matrix
        DEALLOCATE(matrix)

        ! close dataset
        CALL h5dclose_f(dset_id, error)
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
    end function h5load_R1_dp


    ! function ouput = h5load_R2_dp( filename, varname )
    ! Read real rank 2 matrix from h5 file
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname : [string, Input]
    !             variable name in h5 file, must be real numerical rank 2 matrix
    ! Output:
    !   output:   [dp percision real matrix]
    function h5load_R2_dp(filename, varname ) result(matrix)
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), dimension(:,:), allocatable :: matrix
        real(kind=dp), dimension(:,:), allocatable :: temp

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier
        INTEGER(HID_T) :: space_id       ! Dataspace identifier

        INTEGER :: error ! Error flag
        INTEGER :: dim1, dim2 ! matrix dimensions

        INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
        INTEGER(HSIZE_T), DIMENSION(2) :: max_dims


        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

        ! ----------------------- Get Matrix Dimensions -----------------------
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, varname, dset_id, error)
        !Get dataspace ID
        CALL h5dget_space_f(dset_id, space_id, error)
        !Get dataspace dims
        CALL h5sget_simple_extent_dims_f(space_id, data_dims, max_dims, error)
        dim1 = data_dims(1)
        dim2 = data_dims(2)
        CALL h5sclose_f(space_id, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! -------------------------- Allocate Matrix --------------------------
        ALLOCATE( temp(dim1,dim2))

        ! ----------------------------- Read Data -----------------------------
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, varname, dset_id, error)
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F64LE, temp, data_dims, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! ------------------------- Build Real Output -------------------------
        matrix = real( temp, dp)

        ! ------------------------------ Clean Up ------------------------------
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
        ! Deallocate
        DEALLOCATE( temp )

    end function h5load_R2_dp


    ! function ouput = h5load_C2_sp( filename, varname )
    ! Read complex rank 2 matrix from h5 file
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname : [string, Input]
    !             variable name in h5 file, must be complex numerical rank 2 matrix
    ! Output:
    !   output:   [sp percision complex matrix]
    function h5load_C2_sp(filename, varname ) result(matrix)
        character(len=*), intent(in) :: filename, varname
        complex(kind=sp), dimension(:,:), allocatable :: matrix
        real(kind=sp), dimension(:,:), allocatable :: matrix_r, matrix_i

        character(len=100) :: dset_name ! dataset name

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier
        INTEGER(HID_T) :: space_id       ! Dataspace identifier

        INTEGER :: error ! Error flag
        INTEGER :: dim1, dim2 ! matrix dimensions

        INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
        INTEGER(HSIZE_T), DIMENSION(2) :: max_dims


        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

        ! ----------------------- Get Matrix Dimensions -----------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_REAL"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        !Get dataspace ID
        CALL h5dget_space_f(dset_id, space_id, error)
        !Get dataspace dims
        CALL h5sget_simple_extent_dims_f(space_id, data_dims, max_dims, error)
        dim1 = data_dims(1)
        dim2 = data_dims(2)
        CALL h5sclose_f(space_id, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! -------------------------- Allocate Matrix --------------------------
        ALLOCATE( matrix_r(dim1,dim2))
        ALLOCATE( matrix_i(dim1,dim2))

        ! ----------------------------- Real Part -----------------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_REAL"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F32LE, matrix_r, data_dims, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! --------------------------- Imaginary Part ---------------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_IMAG"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F32LE, matrix_i, data_dims, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! ------------------------ Build Complex Output ------------------------
        matrix = CMPLX( matrix_r, matrix_i, sp)

        ! ------------------------------ Clean Up ------------------------------
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
        ! Deallocate
        DEALLOCATE( matrix_r )
        DEALLOCATE( matrix_i )

    end function h5load_C2_sp


    ! function ouput = h5load_C3_sp( filename, varname )
    ! Read complex rank 3 matrix from h5 file
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname : [string, Input]
    !             variable name in h5 file, must be complex numerical rank 3 matrix
    ! Output:
    !   output:   [sp percision complex matrix]
    function h5load_C3_sp(filename, varname ) result(matrix)
        character(len=*), intent(in) :: filename, varname
        complex(kind=sp), dimension(:,:,:), allocatable :: matrix
        real(kind=sp), dimension(:,:,:), allocatable :: matrix_r, matrix_i

        character(len=100) :: dset_name ! dataset name

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier
        INTEGER(HID_T) :: space_id       ! Dataspace identifier

        INTEGER :: error ! Error flag
        INTEGER :: dim1, dim2, dim3 ! matrix dimensions

        INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
        INTEGER(HSIZE_T), DIMENSION(3) :: max_dims


        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

        ! ----------------------- Get Matrix Dimensions -----------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_REAL"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        !Get dataspace ID
        CALL h5dget_space_f(dset_id, space_id, error)
        !Get dataspace dims
        CALL h5sget_simple_extent_dims_f(space_id, data_dims, max_dims, error)
        dim1 = data_dims(1)
        dim2 = data_dims(2)
        dim3 = data_dims(3)
        CALL h5sclose_f(space_id, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! -------------------------- Allocate Matrix --------------------------
        ALLOCATE( matrix_r(dim1,dim2,dim3))
        ALLOCATE( matrix_i(dim1,dim2,dim3))

        ! ----------------------------- Real Part -----------------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_REAL"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F32LE, matrix_r, data_dims, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! --------------------------- Imaginary Part ---------------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_IMAG"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F32LE, matrix_i, data_dims, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! ------------------------ Build Complex Output ------------------------
        matrix = CMPLX( matrix_r, matrix_i, sp)

        ! ------------------------------ Clean Up ------------------------------
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
        ! Deallocate
        DEALLOCATE( matrix_r )
        DEALLOCATE( matrix_i )

    end function h5load_C3_sp


    ! function output = h5load_C4_partial_dp( filename, varname, read_offset, read_dims )
    ! Read a partial complex rank 4 matrix from h5 file
    ! Arguments:
    !   filename   : [string, Input]
    !                h5 filename with path
    !   varname    : [string, Input]
    !                variable name in h5 file, must be complex numerical rank 4 matrix
    !   read_offset: [integer, size(4), Input]
    !                the offset in each dimension, offset = 0 indicates no offset
    !   read_dims  : [integer, size(4), Input]
    !                the dimension of the partial data in each dimension
    ! Output:
    !   output:   [dp percision complex matrix of size read_dims]
    !
    ! Double Precision Only!
    function h5load_C4_partial_dp(filename, varname, read_offset, read_dims ) result(matrix)
        character(len=*), intent(in) :: filename, varname
        integer, dimension(4), intent(in) :: read_offset, read_dims
        complex(kind=dp), dimension(:,:,:,:), allocatable :: matrix
        real(kind=dp), dimension(:,:,:,:), allocatable :: matrix_r, matrix_i

        character(len=100) :: dset_name ! dataset name

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier
        INTEGER(HID_T) :: space_id       ! Dataspace identifier
        INTEGER(HID_T) :: mspace_id      ! Memory space identifier

        INTEGER :: error ! Error flag

        INTEGER(HSIZE_T), DIMENSION(4) :: partial_offset ! offset for the partial data, not offset = 0 is no offset
        INTEGER(HSIZE_T), DIMENSION(4) :: partial_dims ! dimension for the partial data


        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

        ! ------------------- Partial Dimensions and Offset -------------------
        ! Integer type conversion for the partial data dimensions
        partial_dims = INT( read_dims, HSIZE_T)
        partial_offset = INT( read_offset, HSIZE_T)
        ! Create a memory dataspace matching the dimensions of the partial subset
        call h5screate_simple_f(4, partial_dims, mspace_id, error)

        ! -------------------------- Allocate Matrix --------------------------
        ALLOCATE( matrix_r(read_dims(1),read_dims(2),read_dims(3),read_dims(4)))
        ALLOCATE( matrix_i(read_dims(1),read_dims(2),read_dims(3),read_dims(4)))

        ! ----------------------------- Real Part -----------------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_REAL"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        ! Get dataspace ID
        call h5dget_space_f(dset_id, space_id, error)
        ! Select hyperslab with offset
        call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, partial_offset, partial_dims, error)

        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F64LE, matrix_r, partial_dims, error, mspace_id, space_id)

        ! clsoe dataspace
        CALL h5sclose_f(space_id, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! --------------------------- Imaginary Part ---------------------------
        ! dataset name for real part
        dset_name = varname // "/" // varname // "_IMAG"
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, dset_name, dset_id, error)
        ! Get dataspace ID
        call h5dget_space_f(dset_id, space_id, error)
        ! Select hyperslab with offset
        call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, partial_offset, partial_dims, error)

        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F64LE, matrix_i, partial_dims, error, mspace_id, space_id)

        ! clsoe dataspace
        CALL h5sclose_f(space_id, error)
        ! close dataset
        CALL h5dclose_f(dset_id, error)

        ! ------------------------ Build Complex Output ------------------------
        matrix = CMPLX( matrix_r, matrix_i, dp)

        ! ------------------------------ Clean Up ------------------------------
        ! close memory space
        call h5sclose_f(mspace_id, error)
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
        ! Deallocate
        DEALLOCATE( matrix_r )
        DEALLOCATE( matrix_i )

    end function h5load_C4_partial_dp


    ! function output = h5load_check_variable_existence( filename, varname )
    ! Check if a variable exist in the hdf5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name in h5 file
    ! Output:
    !   output:   [boolean] ture if the variable exist
    function h5load_check_variable_existence( filename, varname ) result( existence )
        character(len=*), intent(in) :: filename, varname
        LOGICAL :: existence

        ! h5 variables
        INTEGER(HID_T) :: file_id ! File identifier
        INTEGER :: error ! Error flag
        character(len=100) :: dset_name ! dataset name

        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
        if ( error .eq. -1) then
            write(*,*) "Failed to open file: ", filename
            write(*,*)
            stop
        endif

        ! dset name for the variable
        dset_name = "/" // varname
        ! check if it exist
        CALL h5lexists_f(file_id, dset_name, existence, error)

        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
    end function h5load_check_variable_existence


end module h5load

