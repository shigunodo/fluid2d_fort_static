module fundamentals
  implicit none
  private
  !----------------------------publication--------------------------
  public :: MUSCL3, MP5, MP5_sub, central_diff3, central_diff4
  !-----------------------------------------------------------------
  contains
  !=================================================================
  !----------------3rd order MUSCL reconstruction-------------------
  ! using minmod limitter
  ! need Nb>=2
  !=================================================================
  ! minmod limitter
  subroutine minmod(x, y, ret)
    implicit none
    real(8), intent(in) :: x, y
    real(8), intent(out) :: ret
    ret = sign(max(min(abs(x), sign(y, x*y)), 0.0_8), x)
  end subroutine minmod
  ! 3rd order MUSCL interpolation with minmod limitter
  ! evaluating L/R at cell-boundary between q and qp
  subroutine MUSCL3(qm, q, qp, q2p, qL, qR)
    implicit none
    real(8), intent(in)  :: qm, q, qp, q2p
    real(8), intent(out) :: qL, qR
    real(8) :: b, k
    real(8) :: dp, dm, ddp, ddm
    ! parameters
    k   = 1.0_8/3.0_8
    b   = (3.0_8 - k)/(1.0_8 - k)
    ! calculation of L
    dp = qp - q
    dm = q - qm
    call minmod(dp, b*dm, ddp)
    call minmod(dm, b*dp, ddm)
    qL = q + 0.25_8*((1.0_8 - k)*ddm + (1.0_8 + k)*ddp)
    ! calculation of R
    dp = q2p - qp
    dm = qp - q
    call minmod(dp, b*dm, ddp)
    call minmod(dm, b*dp, ddm)
    qR = qp - 0.25_8*((1.0_8 - k)*ddp + (1.0_8 + k)*ddm)
  end subroutine MUSCL3
  !==================================================================
  !--------5th order Monotonicity-Preserving reconstruction----------
  ! need Nb>=3
  !==================================================================
  ! evaluating qL/R between q and qp
  subroutine minmod4(x1, x2, x3, x4, ret)
    implicit none
    real(8), intent(in)  :: x1, x2, x3, x4
    real(8), intent(out) :: ret
    ret = (sign(0.5_8,x1)+sign(0.5_8,x2)) &
    & *abs((sign(0.5_8,x1)+sign(0.5_8,x3)) &
    & *(sign(0.5_8,x1)+sign(0.5_8,x4))) &
    & *min(abs(x1),abs(x2),abs(x3),abs(x4))
  end subroutine minmod4
  subroutine median(x, y, z, ret)
    implicit none
    real(8), intent(in)  :: x, y, z
    real(8), intent(out) :: ret
    call minmod(y-x,z-x, ret)
    ret = x + ret
  end subroutine median
  subroutine MP5(q2m, qm, q, qp, q2p, q3p, qL, qR)
    implicit none
    real(8), intent(in)  :: q2m, qm, q, qp, q2p, q3p
    real(8), intent(out) :: qL, qR
    call MP5_sub(q2m, qm, q, qp, q2p, qL)
    call MP5_sub(q3p, q2p, qp, q, qm, qR)
  end subroutine MP5
  subroutine MP5_sub(q2m, qm, q, qp, q2p, qL)
    implicit none
    real(8), intent(in)  :: q2m, qm, q, qp, q2p
    real(8), intent(out) :: qL
    real(8), parameter   :: alp = 2.0_8
    real(8)              :: qLL, qMP, dm, d, dp, dMm, dMp
    real(8)              :: qUL, qMD, qLC, qmin, qmax
    qL     = (2.0_8*q2m - 13.0_8*qm + 47.0_8*q + 27.0_8*qp - 3.0_8*q2p)/60.0_8
    call minmod(qp-q, alp*(q-qm), qMP)
    qMP     = q + qMP
    if( (qL-q)*(qL-qMP)>1.d-10 ) then
      qLL     = qL
      dm      = q2m + q   - 2.0_8*qm
      d       = qm  + qp  - 2.0_8*q
      dp      = q   + q2p - 2.0_8*qp
      call minmod4(4.0_8*dm-d, 4.0_8*d-dm, dm, d, dMm)
      call minmod4(4.0_8*d-dp, 4.0_8*dp-d, d, dp, dMp)
      qUL     = q + alp*(q - qm)
      qMD     = 0.5_8*(q+qp) - 0.5_8*dMp
      qLC     = q + 0.5_8*(q-qm) + 4.0_8/3.0_8*dMm
      qmin    = max(min(q,qp,qMD), min(q,qUL,qLC))
      qmax    = min(max(q,qp,qMD), max(q,qUL,qLC))
      call median(qLL, qmin, qmax, qL)
    end if
  end subroutine MP5_sub
  !==================================================================
  !-------------------central difference-----------------------------
  !==================================================================
  ! 3rd order central difference
  ! returning the first derivative at i=1.5
  ! from values at i=0, 1, 2, 3
  subroutine central_diff3(q0, q1, q2, q3, diff)
    implicit none
    real(8), intent(in)  :: q0, q1, q2, q3
    real(8), intent(out) :: diff
    diff = q0/24.0_8 - 9.0_8/8.0_8*q1 + 9.0_8/8.0_8*q2 - q3/24.0_8
  end subroutine central_diff3
  ! 4th order central differende
  ! returning the first derivative at i=2
  ! from values at i=0, 1, 3, 4
  subroutine central_diff4(q0, q1, q3, q4, diff)
    implicit none
    real(8), intent(in)  :: q0, q1, q3, q4
    real(8), intent(out) :: diff
    diff = q0/12.0_8 - 2.0_8/3.0_8*q1 + 2.0_8/3.0_8*q3 - q4/12.0_8
  end subroutine central_diff4
end module fundamentals