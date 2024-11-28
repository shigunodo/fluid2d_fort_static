module ideal_eos
  implicit none
  !=====================================================================
  !------------------ideal gas eos with gamma = 1.4---------------------
  !=====================================================================
  real(8), parameter :: eos_gamma = 1.4_8
  !
  !
  contains
  !
  subroutine calc_p(rho, u, v, e, p)
    implicit none
    real(8), intent(in) :: rho, u, v, e
    real(8), intent(out):: p
    p = (eos_gamma-1.0_8)*(e - 0.5_8*rho*(u*u+v*v))
  end subroutine calc_p
  !
  subroutine calc_T(rho, u, v, e, T)
    implicit none
    real(8), intent(in) :: rho, u, v, e
    real(8), intent(out):: T
    T = e/rho - 0.5_8*(u*u+v*v)
  end subroutine calc_T
  !
  subroutine calc_cs(rho, u, v, e, cs)
    implicit none
    real(8), intent(in) :: rho, u, v, e
    real(8), intent(out):: cs
    cs = sqrt(eos_gamma*(eos_gamma-1.0_8)*(e/rho - 0.5_8*(u*u+v*v)))
  end subroutine calc_cs
  !
  subroutine calc_v(u, v, velo)
    implicit none
    real(8), intent(in) :: u, v
    real(8), intent(out):: velo
    velo = sqrt(u*u+v*v)
  end subroutine calc_v
  !
  subroutine calc_e(rho, u, v, h, e)
    implicit none
    real(8), intent(in)  :: rho, u, v, h
    real(8), intent(out) :: e
    e = rho*(h + (eos_gamma-1.0_8)*0.5_8*(u*u+v*v)) / eos_gamma
  end subroutine calc_e
  !
  subroutine calc_rho_e(p, T, u, v, rho, e)
    implicit none
    real(8), intent(in)  :: p, T, u, v
    real(8), intent(out) :: rho, e
    rho = p/T/(eos_gamma-1.0_8)
    e   = rho*(T + 0.5_8*(u*u+v*v))
  end subroutine calc_rho_e
  !
  subroutine calc_e_wp(rho, u, v, p, e)
    implicit none
    real(8), intent(in)  :: rho, p, u, v
    real(8), intent(out) :: e
    e   = p/(eos_gamma-1.0_8) + 0.5_8*rho*(u*u+v*v)
  end subroutine calc_e_wp
end module ideal_eos