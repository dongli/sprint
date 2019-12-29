module bilinear_interpolator_mod

  use kdtree
  use container
  use string

  implicit none

  type bilinear_interpolator_type
    integer grid_num_x
    integer grid_num_y
    real(8), allocatable :: grid_x(:,:)
    real(8), allocatable :: grid_y(:,:)
    type(kdtree_type) tree
    type(hash_table_type) point_cache
  contains
    procedure, private :: bilinear_interpolator_init_2d_r8
    generic :: init => bilinear_interpolator_init_2d_r8
    procedure :: prepare => bilinear_interpolator_prepare
    final :: bilinear_interpolator_final
  end type bilinear_interpolator_type

contains

  subroutine bilinear_interpolator_init_2d_r8(this, x, y)

    class(bilinear_interpolator_type), intent(out) :: this
    real(8), intent(in) :: x(:,:)
    real(8), intent(in) :: y(:,:)

    this%grid_num_x = size(x, 1)
    this%grid_num_y = size(x, 2)

    this%grid_x = x
    this%grid_y = y

    call this%tree%build(                              &
      reshape(x, [this%grid_num_x * this%grid_num_y]), &
      reshape(y, [this%grid_num_x * this%grid_num_y])  &
    )

    this%point_cache = hash_table(1000)

  end subroutine bilinear_interpolator_init_2d_r8

  subroutine bilinear_interpolator_prepare(this, x, y)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in) :: y

    character(30) key
    type(hash_table_type) cache
    integer ngb_idx(4)

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_cache%hashed(key)) then
      cache = hash_table(10)
      call this%tree%search([x, y], ngb_idx)
    end if

  end subroutine bilinear_interpolator_prepare

  subroutine bilinear_interpolator_final(this)

    type(bilinear_interpolator_type), intent(inout) :: this

    if (allocated(this%grid_x)) deallocate(this%grid_x)
    if (allocated(this%grid_y)) deallocate(this%grid_y)

  end subroutine bilinear_interpolator_final

end module bilinear_interpolator_mod
