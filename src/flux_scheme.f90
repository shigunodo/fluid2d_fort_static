module flux_scheme
  use ideal_eos
  implicit none
  private
  !------------------------------publication---------------------------
  public :: Roe_average, calc_conservative, calc_eigen
  public :: calc_flux_diff, calc_basic
  public :: Roe_FDS ! classical Roe-type FDS scheme
  public :: RoeMAS  ! RoeMAS scheme Qu et al. (2015,2022)
  !--------------------------------------------------------------------
  contains
  !====================================================================
  !----------------------------Roe average-----------------------------
  ! h = (e + p)/rho: enthalpy per mass, cs: local speed of sound
  !====================================================================
  subroutine Roe_average(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                          & rhoA, uA, vA, eA)
    implicit none
    real(8), intent(in)  :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8), intent(out) :: rhoA, uA, vA, eA
    real(8)              :: pL, hL, pR, hR, hA, sqrL, sqrR
    ! preparation
    call calc_p(rhoL, uL, vL, eL, pL)
    call calc_p(rhoR, uR, vR, eR, pR)
    hL   = (eL + pL) / rhoL
    hR   = (eR + pR) / rhoR
    sqrL = sqrt(rhoL)
    sqrR = sqrt(rhoR)
    ! averaging
    rhoA = sqrL*sqrR
    uA   = (sqrL*uL + sqrR*uR) / (sqrL + sqrR)
    vA   = (sqrL*vL + sqrR*vR) / (sqrL + sqrR)
    hA   = (sqrL*hL + sqrR*hR) / (sqrL + sqrR)
    call calc_e(rhoA, uA, vA, hA, eA)
  end subroutine Roe_average
  !====================================================================
  !--------------construction of conservative variables----------------
  ! 2d viscous Navier-Stokes equations system
  !====================================================================
  ! conservative variables
  subroutine calc_conservative(rho, u, v, e, S, Q)
    implicit none
    real(8), intent(in)  :: rho, u, v, e, S
    real(8), intent(out) :: Q(4)
    Q(1) = rho*S
    Q(2) = rho*u*S
    Q(3) = rho*v*S
    Q(4) = e*S
  end subroutine calc_conservative
  ! advective flux
  subroutine calc_flux_conv(rho, u, v, e, ixS, iyS, Fc)
    implicit none
    real(8), intent(in)  :: rho, u, v, e, ixS, iyS
    real(8), intent(out) :: Fc(4)
    real(8)              :: bigU, p
    bigU  = ixS*u + iyS*v
    call calc_p(rho, u, v, e, p)
    Fc(1) = rho*bigU
    Fc(2) = rho*u*bigU + ixS*p
    Fc(3) = rho*v*bigU + iyS*p
    Fc(4) = (e + p)*bigU
  end subroutine calc_flux_conv
  ! diffuison flux
  ! ixS_conv: metrix/J for the calculation of flux
  ! ix_tau: metrix for the calculation of stress tensor
  subroutine calc_flux_diff(Re, Pr, u, v, ixS_conv, iyS_conv, &
                            & ix_tau, iy_tau, jx_tau, jy_tau, &
                            & ui, uj, vi, vj, Ti, Tj, Fd)
    implicit none
    real(8), intent(in)  :: Re, Pr, u, v, ixS_conv, iyS_conv
    real(8), intent(in)  :: ix_tau, iy_tau, jx_tau, jy_tau
    real(8), intent(in)  :: ui, uj, vi, vj, Ti, Tj
    real(8), intent(out) :: Fd(4)
    real(8)              :: tau_xx, tau_xy, tau_yy, b_x, b_y
    real(8)              :: ux, uy, vx, vy, Tx, Ty
    ux     = ui*ix_tau + uj*jx_tau
    uy     = ui*iy_tau + uj*jy_tau
    vx     = vi*ix_tau + vj*jx_tau
    vy     = vi*iy_tau + vj*jy_tau
    Tx     = Ti*ix_tau + Tj*jx_tau
    Ty     = Ti*iy_tau + Tj*jy_tau
    tau_xx = 2.0_8/3.0_8/Re*(2.0_8*ux - vy)
    tau_yy = 2.0_8/3.0_8/Re*(2.0_8*vy - ux)
    tau_xy = (uy + vx)/Re
    b_x    = tau_xx*u + tau_xy*v + Tx/Re/Pr
    b_y    = tau_xy*u + tau_yy*v + Ty/Re/Pr
    Fd(1)  = 0.0_8
    Fd(2)  = ixS_conv*tau_xx + iyS_conv*tau_xy
    Fd(3)  = ixS_conv*tau_xy + iyS_conv*tau_yy
    Fd(4)  = ixS_conv*b_x + iyS_conv*b_y
  end subroutine calc_flux_diff
  ! LAM:  abs of the diagonalized matrix for the flux Jacobian matrix
  ! R:    right eigen matrix of the flux Jacobian matrix
  ! Rinv: inverse of eigen matrix of the flux Jacobian matrix
  ! using ideal gas eos
  subroutine calc_eigen(rho, u, v, e, ix, iy, LAM, R, Rinv)
    implicit none
    real(8), intent(in)  :: rho, u, v, e, ix, iy
    real(8), intent(out) :: LAM(4), R(4,4), Rinv(4,4)
    real(8)              :: ixb, iyb, bigU, bigUb, cs, h, p, sqr, b1, b2
    ! preparaiton
    call calc_cs(rho, u, v, e, cs)
    call calc_p(rho, u, v, e, p)
    h      = (e + p) / rho
    sqr    = sqrt(ix*ix+iy*iy)
    ixb    = ix / sqr
    iyb    = iy / sqr
    bigU   = ix*u + iy*v
    bigUb  = bigU / sqr
    b1     = 0.5*(u*u+v*v)*(eos_gamma-1.0_8)/cs/cs
    b2     = (eos_gamma-1.0_8)/cs/cs
    !
    LAM(1) = bigU - cs*sqr
    LAM(2) = bigU
    LAM(3) = bigU + cs*sqr
    LAM(4) = bigU
    !
    R(1,1) = 1.0_8
    R(1,2) = 1.0_8
    R(1,3) = 1.0_8
    R(1,4) = 0.0_8
    R(2,1) = u - ixb*cs
    R(2,2) = u
    R(2,3) = u + ixb*cs
    R(2,4) = -iyb
    R(3,1) = v - iyb*cs
    R(3,2) = v
    R(3,3) = v + iyb*cs
    R(3,4) = ixb
    R(4,1) = h - cs*bigUb
    R(4,2) = 0.5_8*(u*u+v*v)
    R(4,3) = h + cs*bigUb
    R(4,4) = - (iyb*u - ixb*v)
    !
    Rinv(1,1) = 0.5_8*(b1 + bigUb/cs)
    Rinv(1,2) = -0.5_8*(ixb/cs + b2*u)
    Rinv(1,3) = -0.5_8*(iyb/cs + b2*v)
    Rinv(1,4) = 0.5_8*b2
    Rinv(2,1) = 1.0_8 - b1
    Rinv(2,2) = b2*u
    Rinv(2,3) = b2*v
    Rinv(2,4) = -b2
    Rinv(3,1) = 0.5_8*(b1 - bigUb/cs)
    Rinv(3,2) = 0.5_8*(ixb/cs - b2*u)
    Rinv(3,3) = 0.5_8*(iyb/cs - b2*v)
    Rinv(3,4) = 0.5_8*b2
    Rinv(4,1) = iyb*u - ixb*v
    Rinv(4,2) = -iyb
    Rinv(4,3) = ixb
    Rinv(4,4) = 0.0_8
  end subroutine calc_eigen
  !====================================================================
  !-----------construction of basic vaiables from conservatives--------
  !====================================================================
  subroutine calc_basic(Q, S, rho, u, v, e)
    implicit none
    real(8), intent(in)  :: Q(4), S
    real(8), intent(out) :: rho, u, v, e
    rho = Q(1)/S
    u   = Q(2)/rho/S
    v   = Q(3)/rho/S
    e   = Q(4)/S
  end subroutine calc_basic
  !====================================================================
  !-----------------construction of flux Jacobian----------------------
  !====================================================================
  subroutine calc_flux_jacobian(rho, u, v, e, ix, iy, A)
    implicit none
    real(8), intent(in)  :: rho, u, v, e, ix, iy
    real(8), intent(out) :: A(4,4)
    real(8)              :: bigU, b1, cs
    call calc_cs(rho, u, v, e, cs)
    bigU      = ix*u+iy*v
    b1        = 0.5_8*(eos_gamma-1.0_8)/cs/cs*(u*u+v*v)
    !
    A(1,1)    = 0.0_8
    A(1,2)    = ix
    A(1,3)    = iy
    A(1,4)    = 0.0_8
    A(2,1)    = - u*bigU + ix*b1*cs*cs
    A(2,2)    = bigU - (eos_gamma-2.0_8)*ix*u
    A(2,3)    = iy*u - (eos_gamma-1.0_8)*ix*v
    A(2,4)    = (eos_gamma-1.0_8)*ix
    A(3,1)    = - v*bigU + iy*b1*cs*cs
    A(3,2)    = ix*v - (eos_gamma-1.0_8)*iy*u
    A(3,3)    = bigU - (eos_gamma-2.0_8)*iy*v
    A(3,4)    = (eos_gamma-1.0_8)*iy
    A(4,1)    = - eos_gamma*e/rho*bigU + 2.0_8*b1*cs*cs*bigU
    A(4,2)    = ix*(eos_gamma*e/rho-b1*cs*cs) - (eos_gamma-1.0_8)*u*bigU
    A(4,3)    = iy*(eos_gamma*e/rho-b1*cs*cs) - (eos_gamma-1.0_8)*v*bigU
    A(4,4)    = eos_gamma*bigU
  end subroutine calc_flux_jacobian
    !
  !====================================================================
  !----------------calculation of numerical flux-----------------------
  ! using classical Roe-type FDS scheme
  !====================================================================
  ! ixS, iyS, S should be evaluated at cell-boundaries
  ! by using some averaging method
  subroutine Roe_FDS(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                            & ixS, iyS, S, Fc)
    real(8), intent(in)  :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8), intent(in)  :: ixS, iyS, S
    real(8), intent(out) :: Fc(4)
    real(8)              :: ix, iy, rhoA, uA, vA, eA
    real(8)              :: QR(4), QL(4), FR(4), FL(4)
    real(8)              :: LAM(4), R(4,4), Rinv(4,4)
    integer              :: i, j, l
    ! preparation
    ix     = ixS/S
    iy     = iyS/S
    ! cocnstructing some vectors and matrices
    call Roe_average(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                      & rhoA, uA, vA, eA)
    call calc_conservative(rhoL, uL, vL, eL, S, QL)
    call calc_conservative(rhoR, uR, vR, eR, S, QR)
    call calc_flux_conv(rhoL, uL, vL, eL, ixS, iyS, FL)
    call calc_flux_conv(rhoR, uR, vR, eR, ixS, iyS, FR)
    call calc_eigen(rhoA, uA, vA, eA, ix, iy, LAM, R, Rinv)
    ! calculating numerical flux at cell-boundaries
    do i = 1, 4
      Fc(i) = FR(i) + FL(i)
      do l = 1, 4
        do j = 1, 4
          Fc(i) = Fc(i) - R(i,j)*abs(LAM(j))*Rinv(j,l)*(QR(l) - QL(l))
        end do
      end do
      Fc(i) = Fc(i)*0.5_8
    end do
  end subroutine Roe_FDS
  !====================================================================
  !----------------calculation of numerical flux-----------------------
  ! using RoeMAS scheme, Qu et al. (2015, 2022)
  !====================================================================
  ! ixS, iyS, S should be evaluated at cell-boundaries
  ! by using some averaging method
  subroutine RoeMAS(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                    & ixS, iyS, S, P00, P0p, P0m, Pp0, Ppp, Ppm, Fc)
    real(8), intent(in)  :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8), intent(in)  :: ixS, iyS, S
    real(8), intent(in)  :: P00, P0p, P0m, Pp0, Ppp, Ppm
    real(8), intent(out) :: Fc(4)
    real(8)              :: ix, iy, rhoA, uA, vA, eA
    real(8)              :: sqr, ixb, iyb, M, fM
    real(8)              :: csA, csL, csR, pL, pR, pA
    real(8)              :: h1, P1, l2, l3, Drho, dp, dU
    real(8)              :: FR(4), FL(4)
    real(8), parameter   :: tor = 1.0d-10
    ! preparation
    ix     = ixS/S
    iy     = iyS/S
    sqr    = sqrt(ix*ix+iy*iy)
    ixb    = ix/sqr
    iyb    = iy/sqr
    call Roe_average(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                      & rhoA, uA, vA, eA)
    call calc_cs(rhoA, uA, vA, eA, csA)
    call calc_cs(rhoL, uL, vL, eL, csL)
    call calc_cs(rhoR, uR, vR, eR, csR)
    call calc_p(rhoA, uA, vA, eA, pA)
    call calc_p(rhoL, uL, vL, eL, pL)
    call calc_p(rhoR, uR, vR, eR, pR)
    M      = abs(ixb*uA+iyb*vA)/csA
    !
    P1     = min(P00/Pp0, Pp0/P00)
    h1     = min(P1, min(P00/P0m,P0m/P00), min(P00/P0p,P0p/P00), &
                    & min(Pp0/Ppm,Ppm/Pp0), min(Pp0/Ppp,Ppp/Pp0))
    if(abs(P1-1.0_8)>tor .and. uA*uA+vA*vA>tor .and. M<1.0_8) then
      Drho = csA*sqr*S*(M**(M**(1-h1)))
    else
      Drho = csA*sqr*S*M
    end if
    fM     = min(M, 1.0_8)
    l2     = max(abs(ixS*uA+iyS*vA + fM*sqr*S*csA), &
                & abs(ixS*uR+iyS*vR + fM*sqr*S*csR))
    l3     = min(abs(ixS*uA+iyS*vA - fM*sqr*S*csA), &
                & abs(ixS*uL+iyS*vL - fM*sqr*S*csL))
    dp     = (l2-l3)*0.5_8*(pR-pL) &
          & + (l2+l3)*0.5_8*rhoA*(ixb*uR + iyb*vR - ixb*uL - iyb*vL)
    dU     = ((l2+l3)*0.5_8-abs(ixS*uA+iyS*vA))*(pR-pL)/rhoA/csA/csA &
            & + (l2-l3)*0.5_8*(ixb*uR + iyb*vR - ixb*uL - iyb*vL)/csA
    call calc_flux_conv(rhoL, uL, vL, eL, ixS, iyS, FL)
    call calc_flux_conv(rhoR, uR, vR, eR, ixS, iyS, FR)
    Fc     = 0.5_8*(FR + FL)
    Fc(1)  = Fc(1) - 0.5_8*(Drho*(rhoR-rhoL) + dU*rhoA)
    Fc(2)  = Fc(2) - 0.5_8*(Drho*(rhoR*uR-rhoL*uL) + dU*rhoA*uA + dp*ixb)
    Fc(3)  = Fc(3) - 0.5_8*(Drho*(rhoR*vR-rhoL*vL) + dU*rhoA*vA + dp*iyb)
    Fc(4)  = Fc(4) - 0.5_8*(Drho*(eR + pR - eL - pL) &
                            & + dU*(eA+pA) + dp*(ixb*uA+iyb*vA))
  end subroutine RoeMAS


end module flux_scheme