module io

use global
use hdf5

implicit none

!------------------------------------------------
contains
!------------------------------------------------

! Simple function to check the error status of HDF5 API calls
! ALEXEI: add debugging flag to only use this if debugging
subroutine check_hdf_error(hdferr, message)
  integer :: hdferr
  character(len = *), intent(in) :: message

  if (hdferr == -1) then
    print *, 'Error ', message
    call abort
  endif
end subroutine check_hdf_error

subroutine init_particle_io(particle_id, n_fields, file_id, dset_ids, &
                            write_size)
  integer, intent(in)                       :: particle_id, write_size
  integer, intent(inout)                    :: n_fields
  character(len = 30)                       :: filename 
  integer(hid_t)                            :: file_id, dataspace, crp_list
  integer(hid_t), dimension(:), allocatable :: dset_ids
  integer                                   :: hdferr
  integer, parameter                        :: space_rank = 2
  !integer(hsize_t), dimension(1)           :: max_dim = (/ h5s_unlimited_f /)
  integer(hsize_t), dimension(1)            :: chunk_size = (/ 1 /)
  integer(hsize_t), dimension(2)            :: dims, maxdims
  integer(hsize_t)                          :: offset, count, space_dims
  integer(hid_t)                            :: data_type 

  call h5open_f(hdferr)
  call check_hdf_error(hdferr, 'initialising hdf5')

  write(filename,"('p',I0.10,'.hdf5')") particle_id

  ! Create output file for this particle
  call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, hdferr)
  call check_hdf_error(hdferr, 'opening file')

  ! Create data space
  space_dims = 3
  dims = (/ 3, 0 /)
  maxdims = (/ h5s_unlimited_f, h5s_unlimited_f /)
  call h5screate_simple_f(space_rank, dims, dataspace, hdferr, maxdims)
  call check_hdf_error(hdferr, 'creating dataspace')

  ! Specify chunk size
  call h5pcreate_f(h5p_dataset_create_f, crp_list, hdferr)
  call check_hdf_error(hdferr, 'creating property list')
  call h5pset_chunk_f(crp_list, space_rank, chunk_size, hdferr)
  call check_hdf_error(hdferr, 'creating chunks')

  ! Create datasets
  ! ALEXEI: temporarily specify the number of fields here.
  n_fields = 1
  allocate(dset_ids(n_fields))
  data_type = h5t_native_real
  do i = 1, n_fields
    call h5dcreate_f(file_id, 'particle_data', data_type, dataspace, dset_ids(i), hdferr, crp_list)
    call check_hdf_error(hdferr, 'creating dataset')
    call h5dclose_f(dset_ids(i), hdferr)
    call check_hdf_error(hdferr, 'closing dataset')
  end do
  call h5sclose_f(dataspace, hdferr)
  call check_hdf_error(hdferr, 'closing dataspace')

end subroutine init_particle_io

subroutine write_particle_data(file_id, offset, write_size, &
                               data_R)
  integer(hid_t)    :: file_id
  integer(hid_t)    :: memspace
  integer(hid_t)    :: dataspace
  integer           :: hdferr
  integer(hid_t)    :: dset_id
  integer, intent(inout)  :: offset ! number of elements alread written
  integer, intent(in)  :: write_size ! number of elements being written
  integer(hsize_t), dimension(2)  :: offset_2d ! number of elements alread written (2d)
  integer(hsize_t), dimension(2) :: size_2d ! dimensions of data to write (2d)
  integer(hsize_t), dimension(2) :: extended_size_2d ! dimensions of dataset after extension (2d)
  integer(hsize_t)  :: size_1d  ! dimensions of data to write (1d)
  integer :: space_rank ! number of dimensions in dataspace
  integer(hsize_t), dimension(:), allocatable :: space_dims ! dataspace dimensions
  real(num), dimension(3, write_size), optional :: data_R
 
  ! Reopen the dataset we're writing to
  call h5dopen_f(file_id, 'particle_data', dset_id, hdferr)
  call check_hdf_error(hdferr, 'opening dataset')

  if (present(data_R)) then
    
    extended_size_2d = (/ 3, offset+write_size /)
    size_2d = (/ 3, write_size /)
    offset_2d = (/ 0, offset /)
    call h5dset_extent_f(dset_id, extended_size_2d, hdferr)
    call check_hdf_error(hdferr, 'extending dataset')
    
    space_rank = 2
    allocate(space_dims(space_rank)) ! ALEXEI: do we really need to allocate this every time?
    space_dims = (/ 3, write_size /)
    call h5screate_simple_f(space_rank, space_dims, memspace, hdferr)
    call check_hdf_error(hdferr, 'creating memspace')

    call h5dget_space_f(dset_id, dataspace, hdferr)
    call check_hdf_error(hdferr, 'getting dataspace')

    call h5sselect_hyperslab_f(dataspace, h5s_select_set_f, offset_2d, &
                               size_2d, hdferr)
    call check_hdf_error(hdferr, 'selecting hyperslab')

    call h5dwrite_f(dset_id, h5t_native_real, data_R, size_2d,  &
                    hdferr, memspace, dataspace)
    call check_hdf_error(hdferr, 'writing data')

    deallocate(space_dims)
  end if
  offset = offset + write_size

end subroutine write_particle_data

subroutine close_file(file_id, dset_ids)
  integer :: hdferr
  integer(hid_t) :: file_id
  integer(hid_t), dimension(:), allocatable :: dset_ids

  call h5fclose_f(file_id, hdferr)
  call check_hdf_error(hdferr, 'closing hdf file')

end subroutine close_file

end module io
