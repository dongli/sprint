module bilinear_interpolator_mod

  ! TODO:
  !
  ! - Handle periodic boundary condition.

  use kdtree
  use container
  use string
  use math_mod

  implicit none

  type point_cache_type
    logical :: is_outside = .false.
    integer :: enclose_grid_idx(2,4) = 0
    real(8) :: weights(4) = 0.0d0
  end type point_cache_type

  type bilinear_interpolator_type
    logical :: curvilinear = .false.
    integer grid_nx
    integer grid_ny
    real(8), allocatable :: grid_x(:,:)
    real(8), allocatable :: grid_y(:,:)
    ! Only use KDTree for curvilinear grids.
    type(kdtree_type) tree
    type(hash_table_type) point_caches
  contains
    procedure, private :: bilinear_interpolator_init_2d_r4
    procedure, private :: bilinear_interpolator_init_2d_r8
    procedure, private :: bilinear_interpolator_init_curv2d_r4
    procedure, private :: bilinear_interpolator_init_curv2d_r8
    generic :: init => bilinear_interpolator_init_2d_r4,     &
                       bilinear_interpolator_init_2d_r8,     &
                       bilinear_interpolator_init_curv2d_r4, &
                       bilinear_interpolator_init_curv2d_r8
    procedure, private :: bilinear_interpolator_prepare_r4
    procedure, private :: bilinear_interpolator_prepare_r8
    generic :: prepare => bilinear_interpolator_prepare_r4, &
                          bilinear_interpolator_prepare_r8
    procedure, private :: bilinear_interpolator_is_outside_r4
    procedure, private :: bilinear_interpolator_is_outside_r8
    generic :: is_outside => bilinear_interpolator_is_outside_r4, &
                             bilinear_interpolator_is_outside_r8
    procedure, private :: bilinear_interpolator_get_enclose_grid_idx_r4
    procedure, private :: bilinear_interpolator_get_enclose_grid_idx_r8
    generic :: get_enclose_grid_idx => bilinear_interpolator_get_enclose_grid_idx_r4, &
                                       bilinear_interpolator_get_enclose_grid_idx_r8
    procedure, private :: bilinear_interpolator_apply_2d_r4
    procedure, private :: bilinear_interpolator_apply_2d_r8
    procedure, private :: bilinear_interpolator_apply_3d_r8
    procedure, private :: bilinear_interpolator_apply_4d_r8
    generic :: apply => bilinear_interpolator_apply_2d_r4, &
                        bilinear_interpolator_apply_2d_r8, &
                        bilinear_interpolator_apply_3d_r8, &
                        bilinear_interpolator_apply_4d_r8
    final :: bilinear_interpolator_final
  end type bilinear_interpolator_type

contains

  subroutine bilinear_interpolator_init_2d_r4(this, x, y, cache_size)

    class(bilinear_interpolator_type), intent(out) :: this
    real(4), intent(in) :: x(:)
    real(4), intent(in) :: y(:)
    integer, intent(in), optional :: cache_size

    this%grid_nx = size(x)
    this%grid_ny = size(y)

    if (allocated(this%grid_x)) deallocate(this%grid_x)
    if (allocated(this%grid_y)) deallocate(this%grid_y)
    allocate(this%grid_x(this%grid_nx,1))
    allocate(this%grid_y(this%grid_ny,1))
    this%grid_x(:,1) = x
    this%grid_y(:,1) = y

    if (present(cache_size)) then
      this%point_caches = hash_table(cache_size)
    else
      this%point_caches = hash_table(1000)
    end if

  end subroutine bilinear_interpolator_init_2d_r4

  subroutine bilinear_interpolator_init_2d_r8(this, x, y, cache_size)

    class(bilinear_interpolator_type), intent(out) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(:)
    integer, intent(in), optional :: cache_size

    this%grid_nx = size(x)
    this%grid_ny = size(y)

    if (allocated(this%grid_x)) deallocate(this%grid_x)
    if (allocated(this%grid_y)) deallocate(this%grid_y)
    allocate(this%grid_x(this%grid_nx,1))
    allocate(this%grid_y(this%grid_ny,1))
    this%grid_x(:,1) = x
    this%grid_y(:,1) = y

    if (present(cache_size)) then
      this%point_caches = hash_table(cache_size)
    else
      this%point_caches = hash_table(1000)
    end if

  end subroutine bilinear_interpolator_init_2d_r8

  subroutine bilinear_interpolator_init_curv2d_r4(this, x, y, cache_size)

    class(bilinear_interpolator_type), intent(out) :: this
    real(4), intent(in) :: x(:,:)
    real(4), intent(in) :: y(:,:)
    integer, intent(in), optional :: cache_size

    this%curvilinear = .true.
    this%grid_nx = size(x, 1)
    this%grid_ny = size(x, 2)

    if (allocated(this%grid_x)) deallocate(this%grid_x)
    if (allocated(this%grid_y)) deallocate(this%grid_y)
    allocate(this%grid_x(this%grid_nx,this%grid_ny))
    allocate(this%grid_y(this%grid_nx,this%grid_ny))
    this%grid_x = x
    this%grid_y = y

    call this%tree%build(                                  &
      reshape(this%grid_x, [this%grid_nx * this%grid_ny]), &
      reshape(this%grid_y, [this%grid_nx * this%grid_ny])  &
    )

    if (present(cache_size)) then
      this%point_caches = hash_table(cache_size)
    else
      this%point_caches = hash_table(1000)
    end if

  end subroutine bilinear_interpolator_init_curv2d_r4

  subroutine bilinear_interpolator_init_curv2d_r8(this, x, y, cache_size)

    class(bilinear_interpolator_type), intent(out) :: this
    real(8), intent(in) :: x(:,:)
    real(8), intent(in) :: y(:,:)
    integer, intent(in), optional :: cache_size

    this%curvilinear = .true.
    this%grid_nx = size(x, 1)
    this%grid_ny = size(x, 2)

    if (allocated(this%grid_x)) deallocate(this%grid_x)
    if (allocated(this%grid_y)) deallocate(this%grid_y)
    allocate(this%grid_x(this%grid_nx,this%grid_ny))
    allocate(this%grid_y(this%grid_nx,this%grid_ny))
    this%grid_x = x
    this%grid_y = y

    call this%tree%build(                        &
      reshape(this%grid_x, [this%grid_nx * this%grid_ny]), &
      reshape(this%grid_y, [this%grid_nx * this%grid_ny])  &
    )

    if (present(cache_size)) then
      this%point_caches = hash_table(cache_size)
    else
      this%point_caches = hash_table(1000)
    end if

  end subroutine bilinear_interpolator_init_curv2d_r8

  subroutine bilinear_interpolator_prepare_r4(this, x, y)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(4), intent(in) :: x
    real(4), intent(in) :: y

    character(30) key
    type(point_cache_type) point_cache
    integer i, j, grid_i, grid_j
    real(8) ab(2) ! Normalized coordinate
    ! For curvilinear grids
    integer ngb_idx(12)
    real(8) ngb_dist(12)
    real(8) sorted_ngb_dist(4)
    real(8) grid_x(4), grid_y(4)
    real(8) jacob(2,2), jacobi(2,2)
    ! For equidistant grids
    real(8) dx, dy

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      if (this%curvilinear) then
        call this%tree%search([dble(x), dble(y)], ngb_idx, ngb_dist_=ngb_dist)
        sorted_ngb_dist = -999
        ! Select the enclosing grids and sort them in anti-clockwise order.
        do i = 1, 12
          grid_i = mod(ngb_idx(i) - 1, this%grid_nx) + 1
          grid_j = (ngb_idx(i) - 1) / this%grid_nx + 1
          j = 0
          if (this%grid_x(grid_i,grid_j) <= x .and. this%grid_y(grid_i,grid_j) <= y) then
            if (sorted_ngb_dist(1) < 0 .or. sorted_ngb_dist(1) > ngb_dist(i)) j = 1
          else if (this%grid_x(grid_i,grid_j) >= x .and. this%grid_y(grid_i,grid_j) <= y) then
            if (sorted_ngb_dist(2) < 0 .or. sorted_ngb_dist(2) > ngb_dist(i)) j = 2
          else if (this%grid_x(grid_i,grid_j) >= x .and. this%grid_y(grid_i,grid_j) >= y) then
            if (sorted_ngb_dist(3) < 0 .or. sorted_ngb_dist(3) > ngb_dist(i)) j = 3
          else if (this%grid_x(grid_i,grid_j) <= x .and. this%grid_y(grid_i,grid_j) >= y) then
            if (sorted_ngb_dist(4) < 0 .or. sorted_ngb_dist(4) > ngb_dist(i)) j = 4
          end if
          if (j /= 0) then
            point_cache%enclose_grid_idx(1,j) = grid_i
            point_cache%enclose_grid_idx(2,j) = grid_j
            sorted_ngb_dist(j) = ngb_dist(i)
            grid_x(j) = this%grid_x(grid_i,grid_j)
            grid_y(j) = this%grid_y(grid_i,grid_j)
          end if
        end do
        if (any(point_cache%enclose_grid_idx == 0)) then
          point_cache%is_outside = .true.
        else
          ! Calculate interpolation weights.
          jacob(1,1) = ((grid_x(2) + grid_x(3)) - (grid_x(4) + grid_x(1))) * 0.5 ! dxda
          jacob(1,2) = ((grid_x(3) + grid_x(4)) - (grid_x(1) + grid_x(2))) * 0.5 ! dxdb
          jacob(2,1) = ((grid_y(2) + grid_y(3)) - (grid_y(4) + grid_y(1))) * 0.5 ! dyda
          jacob(2,2) = ((grid_y(3) + grid_y(4)) - (grid_y(1) + grid_y(2))) * 0.5 ! dydb
          jacobi = mat_inv(jacob)
          ab = matmul(jacobi, [x, y]) - matmul(jacobi, [grid_x(1), grid_y(1)])
        end if
      else
        ! We assume grids are equidistant.
        dx = this%grid_x(2,1) - this%grid_x(1,1)
        dy = this%grid_y(2,1) - this%grid_y(1,1)
        grid_i = ceiling((x - this%grid_x(1,1)) / dx)
        grid_j = ceiling((y - this%grid_y(1,1)) / dy)
        if (grid_i < 1 .or. grid_i > this%grid_nx .or. grid_j < 1 .or. grid_j > this%grid_ny) then
          point_cache%is_outside = .true.
        else
          point_cache%enclose_grid_idx(1,1) = grid_i    ; point_cache%enclose_grid_idx(2,1) = grid_j
          point_cache%enclose_grid_idx(1,2) = grid_i + 1; point_cache%enclose_grid_idx(2,2) = grid_j
          point_cache%enclose_grid_idx(1,3) = grid_i + 1; point_cache%enclose_grid_idx(2,3) = grid_j + 1
          point_cache%enclose_grid_idx(1,4) = grid_i    ; point_cache%enclose_grid_idx(2,4) = grid_j + 1
          ab(1) = (x - this%grid_x(grid_i,1)) / dx
          ab(2) = (y - this%grid_y(grid_j,1)) / dy
        end if
      end if
      if (.not. point_cache%is_outside) then
        point_cache%weights(1) = (1 - ab(1)) * (1 - ab(2))
        point_cache%weights(2) = ab(1) * (1 - ab(2))
        point_cache%weights(3) = ab(1) * ab(2)
        point_cache%weights(4) = (1 - ab(1)) * ab(2)
      end if
      call this%point_caches%insert(key, point_cache)
    end if

  end subroutine bilinear_interpolator_prepare_r4

  subroutine bilinear_interpolator_prepare_r8(this, x, y)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in) :: y

    character(30) key
    type(point_cache_type) point_cache
    integer i, j, grid_i, grid_j
    real(8) ab(2)
    ! For curvilinear grids
    integer ngb_idx(12)
    real(8) ngb_dist(12)
    real(8) sorted_ngb_dist(4)
    real(8) grid_x(4), grid_y(4)
    real(8) jacob(2,2), jacobi(2,2)
    ! For equidistant grids
    real(8) dx, dy

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      if (this%curvilinear) then
        call this%tree%search([x, y], ngb_idx, ngb_dist_=ngb_dist)
        sorted_ngb_dist = -999
        ! Select the enclosing grids and sort them in anti-clockwise order.
        do i = 1, 12
          grid_i = mod(ngb_idx(i) - 1, this%grid_nx) + 1
          grid_j = (ngb_idx(i) - 1) / this%grid_nx + 1
          j = 0
          if (this%grid_x(grid_i,grid_j) <= x .and. this%grid_y(grid_i,grid_j) <= y) then
            if (sorted_ngb_dist(1) < 0 .or. sorted_ngb_dist(1) > ngb_dist(i)) j = 1
          else if (this%grid_x(grid_i,grid_j) >= x .and. this%grid_y(grid_i,grid_j) <= y) then
            if (sorted_ngb_dist(2) < 0 .or. sorted_ngb_dist(2) > ngb_dist(i)) j = 2
          else if (this%grid_x(grid_i,grid_j) >= x .and. this%grid_y(grid_i,grid_j) >= y) then
            if (sorted_ngb_dist(3) < 0 .or. sorted_ngb_dist(3) > ngb_dist(i)) j = 3
          else if (this%grid_x(grid_i,grid_j) <= x .and. this%grid_y(grid_i,grid_j) >= y) then
            if (sorted_ngb_dist(4) < 0 .or. sorted_ngb_dist(4) > ngb_dist(i)) j = 4
          end if
          if (j /= 0) then
            point_cache%enclose_grid_idx(1,j) = grid_i
            point_cache%enclose_grid_idx(2,j) = grid_j
            sorted_ngb_dist(j) = ngb_dist(i)
            grid_x(j) = this%grid_x(grid_i,grid_j)
            grid_y(j) = this%grid_y(grid_i,grid_j)
          end if
        end do
        if (any(point_cache%enclose_grid_idx == 0)) then
          point_cache%is_outside = .true.
        else
          ! Calculate interpolation weights.
          jacob(1,1) = ((grid_x(2) + grid_x(3)) - (grid_x(4) + grid_x(1))) * 0.5 ! dxda
          jacob(1,2) = ((grid_x(3) + grid_x(4)) - (grid_x(1) + grid_x(2))) * 0.5 ! dxdb
          jacob(2,1) = ((grid_y(2) + grid_y(3)) - (grid_y(4) + grid_y(1))) * 0.5 ! dyda
          jacob(2,2) = ((grid_y(3) + grid_y(4)) - (grid_y(1) + grid_y(2))) * 0.5 ! dydb
          jacobi = mat_inv(jacob)
          ab = matmul(jacobi, [x, y]) - matmul(jacobi, [grid_x(1), grid_y(1)])
        end if
      else
        ! We assume grids are equidistant.
        dx = this%grid_x(2,1) - this%grid_x(1,1)
        dy = this%grid_y(2,1) - this%grid_y(1,1)
        grid_i = floor((x - this%grid_x(1,1)) / dx)
        grid_j = floor((y - this%grid_y(1,1)) / dy)
        if (grid_i < 1 .or. grid_i > this%grid_nx .or. grid_j < 1 .or. grid_j > this%grid_ny) then
          point_cache%is_outside = .true.
        else
          point_cache%enclose_grid_idx(1,1) = grid_i    ; point_cache%enclose_grid_idx(2,1) = grid_j
          point_cache%enclose_grid_idx(1,2) = grid_i + 1; point_cache%enclose_grid_idx(2,2) = grid_j
          point_cache%enclose_grid_idx(1,3) = grid_i + 1; point_cache%enclose_grid_idx(2,3) = grid_j + 1
          point_cache%enclose_grid_idx(1,4) = grid_i    ; point_cache%enclose_grid_idx(2,3) = grid_j + 1
          ab(1) = (x - this%grid_x(grid_i,1)) / dx
          ab(2) = (y - this%grid_y(grid_j,1)) / dy
        end if
      end if
      if (.not. point_cache%is_outside) then
        point_cache%weights(1) = (1 - ab(1)) * (1 - ab(2))
        point_cache%weights(2) = ab(1) * (1 - ab(2))
        point_cache%weights(3) = ab(1) * ab(2)
        point_cache%weights(4) = (1 - ab(1)) * ab(2)
      end if
      call this%point_caches%insert(key, point_cache)
    end if

  end subroutine bilinear_interpolator_prepare_r8

  logical function bilinear_interpolator_is_outside_r4(this, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(4), intent(in) :: x
    real(4), intent(in) :: y

    character(30) key

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    end if
    select type (point_cache => this%point_caches%value(key))
    type is (point_cache_type)
      res = point_cache%is_outside
    end select

  end function bilinear_interpolator_is_outside_r4

  logical function bilinear_interpolator_is_outside_r8(this, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in) :: y

    character(30) key

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    end if
    select type (point_cache => this%point_caches%value(key))
    type is (point_cache_type)
      res = point_cache%is_outside
    end select

  end function bilinear_interpolator_is_outside_r8

  function bilinear_interpolator_get_enclose_grid_idx_r4(this, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(4), intent(in) :: x
    real(4), intent(in) :: y
    integer res(2,4)

    character(30) key

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    end if
    select type (point_cache => this%point_caches%value(key))
    type is (point_cache_type)
      res = point_cache%enclose_grid_idx
    end select

  end function bilinear_interpolator_get_enclose_grid_idx_r4

  function bilinear_interpolator_get_enclose_grid_idx_r8(this, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    integer res(2,4)

    character(30) key

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    end if
    select type (point_cache => this%point_caches%value(key))
    type is (point_cache_type)
      res = point_cache%enclose_grid_idx
    end select

  end function bilinear_interpolator_get_enclose_grid_idx_r8

  function bilinear_interpolator_apply_2d_r4(this, enclose_values, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(4), intent(in) :: enclose_values(:,:)
    real(4), intent(in) :: x
    real(4), intent(in) :: y
    real(4) res

    character(30) key

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    else
      select type (point_cache => this%point_caches%value(key))
      type is (point_cache_type)
        res = enclose_values(1,1) * point_cache%weights(1) + &
              enclose_values(2,1) * point_cache%weights(2) + &
              enclose_values(1,2) * point_cache%weights(3) + &
              enclose_values(2,2) * point_cache%weights(4)
      end select
    end if

  end function bilinear_interpolator_apply_2d_r4

  function bilinear_interpolator_apply_2d_r8(this, enclose_values, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(8), intent(in) :: enclose_values(:,:)
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8) res

    character(30) key

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    else
      select type (point_cache => this%point_caches%value(key))
      type is (point_cache_type)
        res = enclose_values(1,1) * point_cache%weights(1) + &
              enclose_values(2,1) * point_cache%weights(2) + &
              enclose_values(1,2) * point_cache%weights(3) + &
              enclose_values(2,2) * point_cache%weights(4)
      end select
    end if

  end function bilinear_interpolator_apply_2d_r8

  function bilinear_interpolator_apply_3d_r8(this, enclose_values, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(8), intent(in) :: enclose_values(:,:,:)
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8) res(size(enclose_values, 3))

    character(30) key
    integer k

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    else
      select type (point_cache => this%point_caches%value(key))
      type is (point_cache_type)
        do k = 1, size(enclose_values, 3)
          res(k) = enclose_values(1,1,k) * point_cache%weights(1) + &
                   enclose_values(2,1,k) * point_cache%weights(2) + &
                   enclose_values(1,2,k) * point_cache%weights(3) + &
                   enclose_values(2,2,k) * point_cache%weights(4)
        end do
      end select
    end if

  end function bilinear_interpolator_apply_3d_r8

  function bilinear_interpolator_apply_4d_r8(this, enclose_values, x, y) result(res)

    class(bilinear_interpolator_type), intent(inout) :: this
    real(8), intent(in) :: enclose_values(:,:,:,:)
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8) res(size(enclose_values, 3),size(enclose_values, 4))

    character(30) key
    integer k, l

    key = to_string(x, 4) // ':' // to_string(y, 4)
    if (.not. this%point_caches%hashed(key)) then
      call this%prepare(x, y)
    else
      select type (point_cache => this%point_caches%value(key))
      type is (point_cache_type)
        do l = 1, size(enclose_values, 4)
          do k = 1, size(enclose_values, 3)
            res(k,l) = enclose_values(1,1,k,l) * point_cache%weights(1) + &
                       enclose_values(2,1,k,l) * point_cache%weights(2) + &
                       enclose_values(1,2,k,l) * point_cache%weights(3) + &
                       enclose_values(2,2,k,l) * point_cache%weights(4)
          end do
        end do
      end select
    end if

  end function bilinear_interpolator_apply_4d_r8

  subroutine bilinear_interpolator_final(this)

    type(bilinear_interpolator_type), intent(inout) :: this

    if (allocated(this%grid_x)) deallocate(this%grid_x)
    if (allocated(this%grid_y)) deallocate(this%grid_y)

  end subroutine bilinear_interpolator_final

end module bilinear_interpolator_mod
