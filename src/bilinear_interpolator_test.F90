program bilinear_interpolator_test

  use fiona
  use proj
  use bilinear_interpolator_mod

  implicit none

  type(proj_type) p
  type(bilinear_interpolator_type) interp

  integer wrfout_num_lon, wrfout_num_lat
  real(8) moad_cen_lat, stand_lon, truelat1, truelat2
  real(8), allocatable, dimension(:,:) :: wrfout_xlon, wrfout_xlat
  real(8), allocatable, dimension(:,:) :: wrfout_x, wrfout_y
  integer i, j

  call fiona_init()

  call fiona_create_dataset('wrfout', mode='r', file_path='wrfout')
  call fiona_get_dim('wrfout', 'west_east', size=wrfout_num_lon)
  call fiona_get_dim('wrfout', 'south_north', size=wrfout_num_lat)
  call fiona_get_att('wrfout', 'MOAD_CEN_LAT', moad_cen_lat)
  call fiona_get_att('wrfout', 'STAND_LON', stand_lon)
  call fiona_get_att('wrfout', 'TRUELAT1', truelat1)
  call fiona_get_att('wrfout', 'TRUELAT2', truelat2)
  allocate(wrfout_xlon(wrfout_num_lon,wrfout_num_lat))
  allocate(wrfout_xlat(wrfout_num_lon,wrfout_num_lat))
  allocate(wrfout_x   (wrfout_num_lon,wrfout_num_lat))
  allocate(wrfout_y   (wrfout_num_lon,wrfout_num_lat))
  call fiona_start_input('wrfout')
  call fiona_input('wrfout', 'XLONG', wrfout_xlon)
  call fiona_input('wrfout', 'XLAT' , wrfout_xlat)
  call fiona_end_input('wrfout')

  ! Calculate projected XY coordinates for each WRF grid.
  call p%init(latlon_crs(), lcc_crs(moad_cen_lat, stand_lon, truelat1, truelat2))
  do j = 1, wrfout_num_lat
    do i = 1, wrfout_num_lon
      call p%transform(wrfout_xlon(i,j), wrfout_xlat(i,j), wrfout_x(i,j), wrfout_y(i,j))
    end do
  end do

  deallocate(wrfout_xlon)
  deallocate(wrfout_xlat)
  deallocate(wrfout_x)
  deallocate(wrfout_y)

end program bilinear_interpolator_test
