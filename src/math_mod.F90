module math_mod

  implicit none

  interface mat_inv
    module procedure mat_inv_2x2
  end interface mat_inv

contains

  pure function mat_inv_2x2(A) result(B)

    real(8), intent(in) :: A(2,2)
    real(8) B(2,2)

    real(8) detinv

    detinv = 1.0d0 / (A(1,1) * A(2,2) - A(1,2) * A(2,1))

    B(1,1) =  detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) =  detinv * A(1,1)

  end function mat_inv_2x2

end module math_mod
