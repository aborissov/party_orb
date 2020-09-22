module io

use global
use hdf5

implicit none

integer, parameter :: n_datasets = 2
character(len=5), dimension(n_datasets), parameter :: dataset_names = (/ "R____", "v_par" /)
integer, dimension(n_datasets), parameter :: dataset_dim = (/ 3, 1 /)

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
  !integer(hsize_t), dimension(1)           :: max_dim = (/ h5s_unlimited_f /)
  integer(hsize_t), dimension(1)            :: chunk_size = (/ 1 /)
  integer(hsize_t), dimension(2)            :: dims, maxdims
  integer(hsize_t)                          :: offset, count, space_dims
  integer(hid_t)                            :: data_type 

  ! ALEXEI: Don't do anything here (for profiling)
  !return

  call h5open_f(hdferr)
  call check_hdf_error(hdferr, 'initialising hdf5')

  write(filename,"('p',I0.10,'.hdf5')") particle_id

  ! Create output file for this particle
  call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, hdferr)
  call check_hdf_error(hdferr, 'opening file')
    
  ! Allocate array for dataset ids and specify data type
  allocate(dset_ids(n_datasets))
  data_type = h5t_native_double

  do i = 1, n_datasets

    ! Create data space
    dims = (/ dataset_dim(i), 0 /)
    maxdims = (/ h5s_unlimited_f, h5s_unlimited_f /)
    call h5screate_simple_f(2, dims, dataspace, hdferr, maxdims)
    call check_hdf_error(hdferr, 'creating dataspace')

    ! Specify chunk size
    call h5pcreate_f(h5p_dataset_create_f, crp_list, hdferr)
    call check_hdf_error(hdferr, 'creating property list')
    call h5pset_chunk_f(crp_list, 2, chunk_size, hdferr)
    call check_hdf_error(hdferr, 'creating chunks')

    ! Create the datasets and close them
    call h5dcreate_f(file_id, dataset_names(i), data_type, dataspace, dset_ids(i), hdferr, crp_list)
    call check_hdf_error(hdferr, 'creating dataset')
    call h5dclose_f(dset_ids(i), hdferr)
    call check_hdf_error(hdferr, 'closing dataset')

  end do
  call h5sclose_f(dataspace, hdferr)
  call check_hdf_error(hdferr, 'closing dataspace')

end subroutine init_particle_io

subroutine write_particle_data(file_id, offset, write_size, &
                               data_R, data_v_par)
  integer(hid_t)    :: file_id
  integer(hid_t)    :: memspace
  integer(hid_t)    :: dataspace
  integer           :: hdferr
  integer(hid_t)    :: dset_id
  integer, intent(inout)  :: offset ! number of elements alread written
  integer, intent(in)  :: write_size ! number of elements being written
  real(num), dimension(3, write_size), optional :: data_R
  real(num), dimension(1, write_size), optional :: data_v_par

  ! Debugging
  integer :: i

  ! ALEXEI: Don't do anything here (for profiling)
  !return

  ! Reopen and write to the individual datasets
  if (present(data_R)) then
    call write_array(file_id, 'R____', data_R, offset, write_size, 1)
  end if
  if (present(data_v_par)) then
    call write_array(file_id, 'v_par', data_v_par, offset, write_size, 2)
  end if

  ! Update offset 
  offset = offset + write_size

end subroutine write_particle_data

subroutine write_array(file_id, dset_name, data, offset, n_entries, dataset_index)
  integer(hid_t), intent(in)    :: file_id
  character(len=*), intent(in) :: dset_name
  real(num), dimension(:,:), intent(in) :: data
  integer, intent(in) :: offset
  integer, intent(in) :: n_entries
  integer, intent(in) :: dataset_index
  integer(hsize_t), dimension(2) :: extended_size_2d
  integer(hsize_t), dimension(2) :: size_2d
  integer(hsize_t), dimension(2) :: offset_2d
  integer(hsize_t), dimension(:), allocatable :: space_dims ! dataspace dimensions
  integer(hid_t) :: dset_id
  integer(hid_t) :: dataspace
  integer(hid_t) :: memspace
  integer :: hdferr
    
  call h5dopen_f(file_id, dset_name, dset_id, hdferr)
  call check_hdf_error(hdferr, 'opening dataset')

  extended_size_2d = (/ dataset_dim(dataset_index), offset+n_entries /)
  size_2d = (/ dataset_dim(dataset_index), n_entries /)
  offset_2d = (/ 0, offset /)
  call h5dset_extent_f(dset_id, extended_size_2d, hdferr)
  call check_hdf_error(hdferr, 'extending dataset')
  
  allocate(space_dims(2)) ! ALEXEI: do we really need to allocate this every time?
  space_dims = (/ dataset_dim(dataset_index), n_entries /)
  call h5screate_simple_f(2, space_dims, memspace, hdferr)
  call check_hdf_error(hdferr, 'creating memspace')

  call h5dget_space_f(dset_id, dataspace, hdferr)
  call check_hdf_error(hdferr, 'getting dataspace')

  call h5sselect_hyperslab_f(dataspace, h5s_select_set_f, offset_2d, &
                             size_2d, hdferr)
  call check_hdf_error(hdferr, 'selecting hyperslab')

  call h5dwrite_f(dset_id, h5t_native_double, data, size_2d,  &
                  hdferr, memspace, dataspace)
  call check_hdf_error(hdferr, 'writing data')

  deallocate(space_dims)
  
  ! Close things we opened here
  call h5dclose_f(dset_id, hdferr)
  call check_hdf_error(hdferr, 'closing dataset')
  call h5sclose_f(dataspace, hdferr)
  call check_hdf_error(hdferr, 'closing dataspace')
  call h5sclose_f(memspace, hdferr)
  call check_hdf_error(hdferr, 'closing memspace')

end subroutine write_array

subroutine close_file(file_id, dset_ids, offset, write_size, data_R, data_v_par)
  integer :: hdferr
  integer(hid_t) :: file_id
  integer(hid_t), dimension(:), allocatable :: dset_ids
  integer, intent(in) :: write_size
  integer :: offset
  real(num), dimension(3, write_size), optional :: data_R
  real(num), dimension(1, write_size), optional :: data_v_par

  ! ALEXEI: Don't do anything here (for profiling)
  !return

  ! Write remaining data
  call write_particle_data(file_id, offset, write_size, data_R, data_v_par)

  ! Close the file
  call h5fclose_f(file_id, hdferr)
  call check_hdf_error(hdferr, 'closing hdf file')

end subroutine close_file

end module io
