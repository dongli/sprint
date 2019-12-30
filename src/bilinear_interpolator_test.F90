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
  real(8) lon, lat, x, y, wrfout_values(2,2)
  integer i, j, grid_idx(2,4)

  call fiona_init()

  ! Input test data.
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

  ! Calculate projected XY coordinates for each WRF grid.
  call p%init(latlon_crs(), lcc_crs(moad_cen_lat, stand_lon, truelat1, truelat2))
  do j = 1, wrfout_num_lat
    do i = 1, wrfout_num_lon
      call p%transform(wrfout_xlon(i,j), wrfout_xlat(i,j), wrfout_x(i,j), wrfout_y(i,j))
    end do
  end do

  ! Initialize bilinear interpolator.
  call interp%init(wrfout_x, wrfout_y)

  ! Prepare bilinear interpolation weights.
  lon = 108
  lat = 35
  call p%transform(lon, lat, x, y)
  call interp%prepare(x, y)

  ! Apply bilinear interpolation.
  grid_idx = interp%get_enclose_grid_idx(x, y)
  call fiona_input('wrfout', 'U10', wrfout_values, start=grid_idx(:,1), count=[2,2])
  print *, interp%apply(wrfout_values, x, y)

  call fiona_end_input('wrfout')

  deallocate(wrfout_xlon)
  deallocate(wrfout_xlat)
  deallocate(wrfout_x)
  deallocate(wrfout_y)

end program bilinear_interpolator_test
