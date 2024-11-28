module proceed
  use parameters_io
  use ideal_eos
  use flux_scheme
  use fundamentals
  use bc
  implicit none
  private
  !---------------------public subroutines-----------------------------
  public :: proceed_time
  !
  contains
  !====================================================================
  !--------------------time integrating by Tout------------------------
  ! t: current time
  ! integration is conducted until t = Tout*tstep
  !====================================================================
  subroutine proceed_time(Tout, tstep, t, iter, rho, u, v, e)
    implicit none
    real(8), intent(in)    :: Tout
    integer, intent(in)    :: tstep
    real(8), intent(inout) :: rho(Ni,Nj), u(Ni,Nj), v(Ni,Nj), e(Ni,Nj)
    real(8), intent(inout) :: t
    integer, intent(out)   :: iter
    real(8)                :: t0, dt, CFL, cs, diff, conv
    integer                :: i,j
    t0 = t
    iter = 0
    do while(t < Tout*real(tstep,8))
      ! CFL condition
      ! dx: (Ni-2*Nb, Nj-2*Nb) array
      CFL = 1.0d10
      do j = 1, Nj-2*Nb
        do i = 1, Ni-2*Nb
          call calc_cs(rho(i+Nb, j+Nb), &
                  & u(i+Nb, j+Nb), &
                  & v(i+Nb, j+Nb), &
                  & e(i+Nb, j+Nb), conv)
          conv = dx(i,j) / (conv + sqrt(u(i+Nb, j+Nb)*u(i+Nb, j+Nb) &
                            & + v(i+Nb, j+Nb)*v(i+Nb, j+Nb)))
          diff = 0.1_8*Re*dx(i,j)*dx(i,j)
          CFL = min(CFL, conv)
        end do
      end do
      dt  = 0.7_8 * CFL
      ! integration with 3rd order SSP Runge-Kutta method
      call SSPRK3(dt, rho, u, v, e)
      ! update time variable
      t = t + dt
      iter = iter + 1
    end do
  end subroutine proceed_time
  !====================================================================
  !------------------3rd order SSP Runge-Kutta method------------------
  !====================================================================
  subroutine SSPRK3(dt, rho, u, v, e)
    implicit none
    real(8), intent(in)    :: dt
    real(8), intent(inout) :: rho(Ni,Nj), u(Ni,Nj), v(Ni,Nj), e(Ni,Nj)
    integer                :: i, j
    real(8)                :: SA
    real(8), save :: Q0(4,Ni-2*Nb,Nj-2*Nb), Q1(4,Ni-2*Nb,Nj-2*Nb), Q2(4,Ni-2*Nb,Nj-2*Nb)
    ! constructing conservative variable
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb
        ! Jacobian is evaluated at (i, j) by average of adjacent 4 cells
        ! parameter S is (Ni-1, Nj-1) array
        SA = 0.25_8*(S(Nb+i-1,Nb+j-1) + S(Nb+i,Nb+j-1) + S(Nb+i-1,Nb+j) + S(Nb+i,Nb+j))
        call calc_conservative(rho(Nb+i,Nb+j), u(Nb+i,Nb+j), &
                      & v(Nb+i,Nb+j), e(Nb+i,Nb+j), SA, Q0(:,i,j))
      end do
    end do
    ! 1st stage
    call RHS(rho, u, v, e, Q1)
    Q1 = Q0 + dt*Q1
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb
        SA = 0.25_8*(S(Nb+i-1,Nb+j-1) + S(Nb+i,Nb+j-1) + S(Nb+i-1,Nb+j) + S(Nb+i,Nb+j))
        call calc_basic(Q1(:,i,j), SA, rho(Nb+i,Nb+j), u(Nb+i,Nb+j), &
                        & v(Nb+i,Nb+j), e(Nb+i,Nb+j))
      end do
    end do
    ! reflecting boundary condition
    call bc_periodic(rho, u, v, e)
    ! 2nd stage
    call RHS(rho, u, v, e, Q2)
    Q2 = 0.75_8*Q0 + 0.25_8*(Q1 + dt*Q2)
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb
        SA = 0.25_8*(S(Nb+i-1,Nb+j-1) + S(Nb+i,Nb+j-1) + S(Nb+i-1,Nb+j) + S(Nb+i,Nb+j))
        call calc_basic(Q2(:,i,j), SA, rho(Nb+i,Nb+j), u(Nb+i,Nb+j), &
                        & v(Nb+i,Nb+j), e(Nb+i,Nb+j))
      end do
    end do
    ! reflecting boundary condition
    call bc_periodic(rho, u, v, e)
    ! 3rd stage
    call RHS(rho, u, v, e, Q1)
    Q0 = Q0/3.0_8 + 2.0_8/3.0_8*(Q2 + dt*Q1)
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb
        SA = 0.25_8*(S(Nb+i-1,Nb+j-1) + S(Nb+i,Nb+j-1) + S(Nb+i-1,Nb+j) + S(Nb+i,Nb+j))
        call calc_basic(Q0(:,i,j), SA, rho(Nb+i,Nb+j), u(Nb+i,Nb+j), &
                        & v(Nb+i,Nb+j), e(Nb+i,Nb+j))
      end do
    end do
    ! reflecting boundary condition
    call bc_periodic(rho, u, v, e)
  end subroutine SSPRK3
  !====================================================================
  !------------------2nd order SSP Runge-Kutta method------------------
  !====================================================================
  subroutine SSPRK2(dt, rho, u, v, e)
    implicit none
    real(8), intent(in)    :: dt
    real(8), intent(inout) :: rho(Ni,Nj), u(Ni,Nj), v(Ni,Nj), e(Ni,Nj)
    integer                :: i, j
    real(8)                :: SA
    real(8), save :: Q0(4,Ni-2*Nb,Nj-2*Nb), Q1(4,Ni-2*Nb,Nj-2*Nb), Q2(4,Ni-2*Nb,Nj-2*Nb)
    ! constructing conservative variable
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb
        ! Jacobian is evaluated at (i, j) by average of adjacent 4 cells
        ! parameter S is (Ni-1, Nj-1) array
        SA = 0.25_8*(S(Nb+i-1,Nb+j-1) + S(Nb+i,Nb+j-1) + S(Nb+i-1,Nb+j) + S(Nb+i,Nb+j))
        call calc_conservative(rho(Nb+i,Nb+j), u(Nb+i,Nb+j), &
                      & v(Nb+i,Nb+j), e(Nb+i,Nb+j), SA, Q0(:,i,j))
      end do
    end do
    ! 1st stage
    call RHS(rho, u, v, e, Q1)
    Q1 = Q0 + dt*Q1
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb
        SA = 0.25_8*(S(Nb+i-1,Nb+j-1) + S(Nb+i,Nb+j-1) + S(Nb+i-1,Nb+j) + S(Nb+i,Nb+j))
        call calc_basic(Q1(:,i,j), SA, rho(Nb+i,Nb+j), u(Nb+i,Nb+j), &
                        & v(Nb+i,Nb+j), e(Nb+i,Nb+j))
      end do
    end do
    ! 2nd stage
    call RHS(rho, u, v, e, Q2)
    Q2 = 0.5_8*Q0 + 0.5_8*(Q1 + dt*Q2)
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb
        SA = 0.25_8*(S(Nb+i-1,Nb+j-1) + S(Nb+i,Nb+j-1) + S(Nb+i-1,Nb+j) + S(Nb+i,Nb+j))
        call calc_basic(Q2(:,i,j), SA, rho(Nb+i,Nb+j), u(Nb+i,Nb+j), &
                        & v(Nb+i,Nb+j), e(Nb+i,Nb+j))
      end do
    end do
  end subroutine SSPRK2
  !====================================================================
  !-------------------------calculating RHS----------------------------
  !====================================================================
  subroutine RHS(rho, u, v, e, ret)
    implicit none
    integer,parameter    :: Ns = 6
    real(8), intent(in)  :: rho(Ni,Nj), u(Ni,Nj), v(Ni, Nj), e(Ni,Nj)
    real(8), intent(out) :: ret(4,Ni-2*Nb,Nj-2*Nb)
    integer              :: i, j, l
    real(8)              :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8)              :: ixSA, iySA, jxSA, jySA, SA
    real(8), save        :: Fi(4,Ni-2*Nb+1,Nj-2*Nb), Fj(4,Ni-2*Nb,Nj-2*Nb+1)
    ! for reconstruction by conservative or characteristic variables
    real(8)              :: S_s(Ns)
    ! for diffuison flux
    real(8), save        :: temp(Ni,Nj)
    real(8)              :: Fd(4), flag_i, flag_j
    ! for RoeMAS scheme
    real(8), save        :: p(Ni,Nj)
    ! for test problems
    integer              :: con_sch, ma_sch
    ! calculating numerical fluxes
    ! returning (Ni-2*Nb+1, Nj-2*Nb) array
    ! first, calculating temperature and pressure
    do j = 1, Nj
      do i = 1, Ni
        call calc_T(rho(i,j), u(i,j), v(i,j), e(i,j), temp(i,j)) !for diffusion
        call calc_p(rho(i,j), u(i,j), v(i,j), e(i,j), p(i,j)) ! for RoeMAS
      end do
    end do
    !
    con_sch = 4
    ma_sch  = 1
    !$omp parallel private(SA, ixSA, iySA, jxSA, jySA, & !S_s, &
    !$omp & rhoL, uL, vL, eL, rhoR, uR, vR, eR, flag_i, flag_j, Fd)
    !$omp do
    do j = 1, Nj-2*Nb
      do i = 1, Ni-2*Nb+1
        !----------------------------------------------------------------
        !----------------------------i-direction-------------------------
        !----------------------------------------------------------------
        ! evaluating flux at (i+0.5,j) from stencil [(i-1,j),...,(i+2,j)]
        !
        !
        !-------------------calculating metrices------------------------
        ! inverse of Jacobian
        ! evaluated at (i+0.5,j)
        ! by averaging the values at (i+0.5,j+0.5) and (i+0.5,j-0.5)
        SA     = 0.5_8*(S(Nb+i-1, Nb+j-1) + S(Nb+i-1, Nb+j))
        ! metices/J
        ! evaluated at (i+0.5,j)
        ! by averaging the values at (i,j) and (i+1,j)
        ixSA   = 0.5_8*(ixS(i, j+1) + ixS(i+1, j+1))
        iySA   = 0.5_8*(iyS(i, j+1) + iyS(i+1, j+1))
        jxSA   = 0.5_8*(jxS(i, j+1) + jxS(i+1, j+1))
        jySA   = 0.5_8*(jyS(i, j+1) + jyS(i+1, j+1))
        ! Jacobian stencil for reconstruction evaluated at (i,j)
        ! for reconstruction by conservative and characteristic variables
        do l = 1, Ns
          S_s(l) = 0.25_8*(S(Nb-Ns/2+i-2+l,Nb+j-1) + &
                          & S(Nb-Ns/2+i-2+l,Nb+j) + &
                          & S(Nb-Ns/2+i-1+l,Nb+j-1) + &
                          & S(Nb-Ns/2+i-1+l,Nb+j))
        end do
        !
        !---------------------------reconstruction-----------------------
        ! reconstruction scheme
        !==================3rd order MUSCL interpolation=================
        if(con_sch == 1) then
          call reconstruction_basic_MUSCL3( &
                  & rho(Nb+i-2:Nb+i+1,Nb+j), u(Nb+i-2:Nb+i+1,Nb+j), &
                  & v(Nb+i-2:Nb+i+1,Nb+j), e(Nb+i-2:Nb+i+1,Nb+j), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 2) then
          call reconstruction_cons_MUSCL3( &
                  & rho(Nb+i-2:Nb+i+1,Nb+j), u(Nb+i-2:Nb+i+1,Nb+j), &
                  & v(Nb+i-2:Nb+i+1,Nb+j), e(Nb+i-2:Nb+i+1,Nb+j), &
                  & S_s(Ns/2-1:Ns/2+2), SA, &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 3) then
          call reconstruction_char_MUSCL3( &
                  & rho(Nb+i-2:Nb+i+1,Nb+j), u(Nb+i-2:Nb+i+1,Nb+j), &
                  & v(Nb+i-2:Nb+i+1,Nb+j), e(Nb+i-2:Nb+i+1,Nb+j), &
                  & S_s(Ns/2-1:Ns/2+2), SA, &
                  & ixS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & iyS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 7) then
          call reconstruction_flux_MUSCL3( &
                  & rho(Nb+i-2:Nb+i+1,Nb+j), u(Nb+i-2:Nb+i+1,Nb+j), &
                  & v(Nb+i-2:Nb+i+1,Nb+j), e(Nb+i-2:Nb+i+1,Nb+j), &
                  & S_s(Ns/2-1:Ns/2+2), SA, &
                  & ixS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & iyS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        !===================5th order MP scheme==========================
        else if(con_sch == 4) then
          call reconstruction_basic_MP5( &
                  & rho(Nb+i-3:Nb+i+2,Nb+j), u(Nb+i-3:Nb+i+2,Nb+j), &
                  & v(Nb+i-3:Nb+i+2,Nb+j), e(Nb+i-3:Nb+i+2,Nb+j), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 5) then
          call reconstruction_cons_MP5( &
                  & rho(Nb+i-3:Nb+i+2,Nb+j), u(Nb+i-3:Nb+i+2,Nb+j), &
                  & v(Nb+i-3:Nb+i+2,Nb+j), e(Nb+i-3:Nb+i+2,Nb+j), &
                  & S_s(Ns/2-2:Ns/2+3), SA, &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 6) then
          call reconstruction_char_MP5( &
                  & rho(Nb+i-3:Nb+i+2,Nb+j), u(Nb+i-3:Nb+i+2,Nb+j), &
                  & v(Nb+i-3:Nb+i+2,Nb+j), e(Nb+i-3:Nb+i+2,Nb+j), &
                  & S_s(Ns/2-2:Ns/2+3), SA, &
                  & ixS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & iyS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        end if
        !================================================================
        !
        !------------------calculating convection flux---------------
        ! numerical flux scheme
        !======================Roe-type FDS=========================
        if( ma_sch == 1) then
          call Roe_FDS(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                            & ixSA, iySA, SA, Fi(:,i,j))
        !=========================RoeMAS============================
        else if( ma_sch== 2) then
          call RoeMAS(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                      & ixSA, iySA, SA, p(Nb+i,Nb+j), p(Nb+i,Nb+j+1), &
                      & p(Nb+i,Nb+j-1), p(Nb+i+1,Nb+j), &
                      & p(Nb+i+1,Nb+j+1), p(Nb+i+1,Nb+j-1), Fi(:,i,j))
        end if
        !===========================================================
        ! Qiita の記事ではオイラー方程式を解くため、粘性項は計算しない
        !------------------calculating diffusion flux-----------------
        !flag_i = 1.0_8
        !flag_j = 0.0_8
        !call construct_numflux_diff(Re, Pr, rhoL, uL, vL, eL, &
        !    & rhoR, uR, vR, eR, ixSA, iySA, jxSA, jySA, SA, &
        !    & u(Nb+i-2:Nb+i+1,Nb+j), v(Nb+i-2:Nb+i+1,Nb+j), &
        !    & temp(Nb+i-2:Nb+i+1,Nb+j), &
        !    & u(Nb+i-1:Nb+i,Nb+j-2:Nb+j+2), v(Nb+i-1:Nb+i,Nb+j-2:Nb+j+2), &
        !    & temp(Nb+i-1:Nb+i,Nb+j-2:Nb+j+2), &
        !    & flag_i, flag_j, Fd)
        !Fi(:,i,j) = Fi(:,i,j) - Fd
      end do
    end do
    !$omp end do
    !
    !$omp do
    do j = 1, Nj-2*Nb+1
      do i = 1, Ni-2*Nb
        !----------------------------------------------------------------
        !----------------------------j-direction-------------------------
        !----------------------------------------------------------------
        ! evaluating flux at (i,j+0.5) from stencil [(i,j-1),...,(i,j+2)]
        !
        !
        !---------------------calculating metrices-------------------
        ! inverse of Jacobian
        ! evaluated at (i,j+0.5)
        ! by averaging the values at (i+0.5,j+0.5) and (i-0.5,j+0.5)
        SA     = 0.5_8*(S(Nb+i-1, Nb+j-1) + S(Nb+i, Nb+j-1))
        ! metices/J
        ! evaluated at (i,j+0.5)
        ! by averaging the values at (i,j) and (i,j+1)
        ixSA   = 0.5_8*(ixS(i+1, j) + ixS(i+1, j+1))
        iySA   = 0.5_8*(iyS(i+1, j) + iyS(i+1, j+1))
        jxSA   = 0.5_8*(jxS(i+1, j) + jxS(i+1, j+1))
        jySA   = 0.5_8*(jyS(i+1, j) + jyS(i+1, j+1))
        ! Jacobian stencil for reconstruction evaluated at (i,j)
        ! for reconstruction by conservative and characteristic variables
        do l = 1, Ns
          S_s(l) = 0.25_8*(S(Nb+i-1,Nb-Ns/2+j-2+l) + &
                          & S(Nb+i,Nb-Ns/2+j-2+l) + &
                          & S(Nb+i-1,Nb-Ns/2+j-1+l) + &
                          & S(Nb+i,Nb-Ns/2+j-1+l))
        end do
        !
        !------------------------reconstruction--------------------------
        ! reconstruction scheme
        !==================3rd order MUSCL interpolation=================
        if(con_sch == 1) then
          call reconstruction_basic_MUSCL3( &
                  & rho(Nb+i,Nb+j-2:Nb+j+1), u(Nb+i,Nb+j-2:Nb+j+1), &
                  & v(Nb+i,Nb+j-2:Nb+j+1), e(Nb+i,Nb+j-2:Nb+j+1), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 2) then
          call reconstruction_cons_MUSCL3( &
                  & rho(Nb+i,Nb+j-2:Nb+j+1), u(Nb+i,Nb+j-2:Nb+j+1), &
                  & v(Nb+i,Nb+j-2:Nb+j+1), e(Nb+i,Nb+j-2:Nb+j+1), &
                  & S_s(Ns/2-1:Ns/2+2), SA, &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 3) then
          call reconstruction_char_MUSCL3( &
                  & rho(Nb+i,Nb+j-2:Nb+j+1), u(Nb+i,Nb+j-2:Nb+j+1), &
                  & v(Nb+i,Nb+j-2:Nb+j+1), e(Nb+i,Nb+j-2:Nb+j+1), &
                  & S_s(Ns/2-1:Ns/2+2), SA, &
                  & jxS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & jyS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 7) then
          call reconstruction_flux_MUSCL3( &
                  & rho(Nb+i,Nb+j-2:Nb+j+1), u(Nb+i,Nb+j-2:Nb+j+1), &
                  & v(Nb+i,Nb+j-2:Nb+j+1), e(Nb+i,Nb+j-2:Nb+j+1), &
                  & S_s(Ns/2-1:Ns/2+2), SA, &
                  & jxS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & jyS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        !===================5th order MP scheme==========================
        else if(con_sch == 4) then
          call reconstruction_basic_MP5( &
                  & rho(Nb+i,Nb+j-3:Nb+j+2), u(Nb+i,Nb+j-3:Nb+j+2), &
                  & v(Nb+i,Nb+j-3:Nb+j+2), e(Nb+i,Nb+j-3:Nb+j+2), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 5) then
          call reconstruction_cons_MP5( &
                  & rho(Nb+i,Nb+j-3:Nb+j+2), u(Nb+i,Nb+j-3:Nb+j+2), &
                  & v(Nb+i,Nb+j-3:Nb+j+2), e(Nb+i,Nb+j-3:Nb+j+2), &
                  & S_s(Ns/2-2:Ns/2+3), SA, &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        else if(con_sch == 6) then
          call reconstruction_char_MP5( &
                  & rho(Nb+i,Nb+j-3:Nb+j+2), u(Nb+i,Nb+j-3:Nb+j+2), &
                  & v(Nb+i,Nb+j-3:Nb+j+2), e(Nb+i,Nb+j-3:Nb+j+2), &
                  & S_s(Ns/2-2:Ns/2+3), SA, &
                  & jxS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & jyS(i:i+1,j+1)/S_s(Ns/2:Ns/2+1), &
                  & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
        end if
        !================================================================
        !
        !------------------calculating convection flux---------------
        ! numerical flux scheme
        !=======================Roe-type FDS=======================
        if(ma_sch == 1) then
          call Roe_FDS(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                            & jxSA, jySA, SA, Fj(:,i,j))
        !=========================RoeMAS============================
        else if(ma_sch == 2) then
          call RoeMAS(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
                    & jxSA, jySA, SA, p(Nb+i,Nb+j), p(Nb+i+1,Nb+j), &
                    & p(Nb+i-1,Nb+j), p(Nb+i,Nb+j+1), &
                    & p(Nb+i+1,Nb+j+1), p(Nb+i-1,Nb+j+1), Fj(:,i,j))
        end if
        !===========================================================
        ! Qiita の記事ではオイラー方程式を解くため、粘性項は計算しない
        !-----------------calculating diffusion flux----------------
        !flag_i = 0.0_8
        !flag_j = 1.0_8
        !call construct_numflux_diff(Re, Pr, rhoL, uL, vL, eL, &
        !    & rhoR, uR, vR, eR, ixSA, iySA, jxSA, jySA, SA, &
        !    & u(Nb+i,Nb+j-2:Nb+j+1), v(Nb+i,Nb+j-2:Nb+j+1), &
        !    & temp(Nb+i,Nb+j-2:Nb+j+1), &
        !    & transpose(u(Nb+i-2:Nb+i+2,Nb+j-1:Nb+j)), &
        !    & transpose(v(Nb+i-2:Nb+i+2,Nb+j-1:Nb+j)), &
        !    & transpose(temp(Nb+i-2:Nb+i+2,Nb+j-1:Nb+j)), &
        !    & flag_i, flag_j, Fd)
        !Fj(:,i,j) = Fj(:,i,j) - Fd
      end do
    end do
    !$omp end do
    !$omp end parallel
    !
    !------------------------calculating RHS------------------------
    ret     = Fi(:,1:Ni-2*Nb,1:Nj-2*Nb) - Fi(:,2:Ni-2*Nb+1,1:Nj-2*Nb) &
              & + Fj(:,1:Ni-2*Nb,1:Nj-2*Nb) - Fj(:,1:Ni-2*Nb,2:Nj-2*Nb+1)
  end subroutine RHS
  !======================================================================
  !------------------reconstruction by basic variables-------------------
  !======================================================================
  subroutine reconstruction_basic_MUSCL3(rho_s, u_s, v_s, e_s, &
                          & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
    implicit none
    real(8), intent(in)  :: rho_s(4), u_s(4), v_s(4), e_s(4)
    real(8), intent(out) :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    call MUSCL3(rho_s(1), rho_s(2), rho_s(3), rho_s(4), &
                              & rhoL, rhoR)
    call MUSCL3(u_s(1), u_s(2), u_s(3), u_s(4), &
                              & uL, uR)
    call MUSCL3(v_s(1), v_s(2), v_s(3), v_s(4), &
                              & vL, vR)
    call MUSCL3(e_s(1), e_s(2), e_s(3), e_s(4), &
                              & eL, eR)
  end subroutine reconstruction_basic_MUSCL3
  subroutine reconstruction_basic_MP5(rho_s, u_s, v_s, e_s, &
                          & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
    implicit none
    real(8), intent(in)  :: rho_s(6), u_s(6), v_s(6), e_s(6)
    real(8), intent(out) :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    call MP5(rho_s(1), rho_s(2), rho_s(3), rho_s(4), rho_s(5), rho_s(6),&
                          & rhoL, rhoR)
    call MP5(u_s(1), u_s(2), u_s(3), u_s(4), u_s(5), u_s(6),&
                          & uL, uR)
    call MP5(v_s(1), v_s(2), v_s(3), v_s(4), v_s(5), v_s(6),&
                          & vL, vR)
    call MP5(e_s(1), e_s(2), e_s(3), e_s(4), e_s(5), e_s(6),&
                          & eL, eR)
  end subroutine reconstruction_basic_MP5
  !======================================================================
  !---------------reconstruction by conservative variables---------------
  !======================================================================
  subroutine reconstruction_cons_MUSCL3(rho_s, u_s, v_s, e_s, S_s, SA, &
                          & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
    implicit none
    real(8), intent(in)  :: rho_s(4), u_s(4), v_s(4), e_s(4), S_s(4)
    real(8), intent(in)  :: SA
    real(8), intent(out) :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8)              :: Q(4,4), QL(4), QR(4)
    integer              :: i
    do i = 1, 4
      call calc_conservative(rho_s(i), u_s(i), v_s(i), e_s(i), S_s(i), Q(:,i))
    end do
    do i = 1, 4
      call MUSCL3(Q(i,1), Q(i,2), Q(i,3), Q(i,4), QL(i), QR(i))
    end do
    call calc_basic(QL, SA, rhoL, uL, vL, eL)
    call calc_basic(QR, SA, rhoR, uR, vR, eR)
  end subroutine reconstruction_cons_MUSCL3
  subroutine reconstruction_cons_MP5(rho_s, u_s, v_s, e_s, S_s, SA, &
                          & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
  implicit none
  real(8), intent(in)  :: rho_s(6), u_s(6), v_s(6), e_s(6), S_s(6)
  real(8), intent(in)  :: SA
  real(8), intent(out) :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
  real(8)              :: Q(4,6), QL(4), QR(4)
  integer              :: i
  do i = 1, 6
  call calc_conservative(rho_s(i), u_s(i), v_s(i), e_s(i), S_s(i), Q(:,i))
  end do
  do i = 1, 4
  call MP5(Q(i,1), Q(i,2), Q(i,3), Q(i,4), Q(i,5), Q(i,6), QL(i), QR(i))
  end do
  call calc_basic(QL, SA, rhoL, uL, vL, eL)
  call calc_basic(QR, SA, rhoR, uR, vR, eR)
  end subroutine reconstruction_cons_MP5
  !======================================================================
  !-------------reconstruction by characteristic variables---------------
  !======================================================================
  subroutine reconstruction_char_MUSCL3(rho_s, u_s, v_s, e_s, S_s, SA, &
                          & ix_s, iy_s, &
                          & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
    implicit none
    integer, parameter   :: Nv = 4, Ns = 4
    real(8), intent(in)  :: rho_s(Ns), u_s(Ns), v_s(Ns), e_s(Ns), S_s(Ns)
    real(8), intent(in)  :: SA, ix_s(2), iy_s(2)
    real(8), intent(out) :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8)              :: Q(Nv,Ns), LAM(Nv), R(Nv,Nv), Rinv(Nv,Nv)
    real(8)              :: W(Ns), Wc(Nv), tmp, Qc(Nv)
    integer              :: i, j, l
    do i = 1, Ns
      call calc_conservative(rho_s(i), u_s(i), v_s(i), e_s(i), S_s(i), Q(:,i))
    end do
    ! QL
    call calc_eigen(rho_s(2), u_s(2), v_s(2), e_s(2), ix_s(1), iy_s(1), &
                          & LAM, R, Rinv)
    do l = 1, Nv ! index for char-variables
      W = 0.0_8
      do j = 1, Ns ! index for stencil-size
        do i = 1, Nv ! index for cons-variables
          W(j) = W(j) + Rinv(l,i)*Q(i,j)
        end do
      end do
      call MUSCL3(W(1), W(2), W(3), W(4), Wc(l), tmp)
    end do
    Qc = 0.0_8
    do j = 1, Nv
      do i = 1, Nv
        Qc(i) = Qc(i) + R(i,j)*Wc(j)
      end do
    end do
    call calc_basic(Qc, SA, rhoL, uL, vL, eL)
    ! QR
    call calc_eigen(rho_s(3), u_s(3), v_s(3), e_s(3), ix_s(2), iy_s(2), &
                          & LAM, R, Rinv)
    do l = 1, Nv ! index for char-variables
      W = 0.0_8
      do j = 1, Ns ! index for stencil-size
        do i = 1, Nv ! index for cons-variables
          W(j) = W(j) + Rinv(l,i)*Q(i,j)
        end do
      end do
      call MUSCL3(W(1), W(2), W(3), W(4), tmp, Wc(l))
    end do
    Qc = 0.0_8
    do j = 1, Nv
      do i = 1, Nv
        Qc(i) = Qc(i) + R(i,j)*Wc(j)
      end do
    end do
    call calc_basic(Qc, SA, rhoR, uR, vR, eR)
  end subroutine reconstruction_char_MUSCL3
  subroutine reconstruction_char_MP5(rho_s, u_s, v_s, e_s, S_s, SA, &
    & ix_s, iy_s, &
    & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
    implicit none
    integer, parameter   :: Nv = 4, Ns = 6
    real(8), intent(in)  :: rho_s(Ns), u_s(Ns), v_s(Ns), e_s(Ns), S_s(Ns)
    real(8), intent(in)  :: SA, ix_s(2), iy_s(2)
    real(8), intent(out) :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8)              :: Q(Nv,Ns), LAM(Nv), R(Nv,Nv), Rinv(Nv,Nv)
    real(8)              :: W(Ns), Wc(Nv), Qc(Nv)
    integer              :: i, j, l
    do i = 1, Ns
      call calc_conservative(rho_s(i), u_s(i), v_s(i), e_s(i), S_s(i), Q(:,i))
    end do
    ! QL
    call calc_eigen(rho_s(2), u_s(2), v_s(2), e_s(2), ix_s(1), iy_s(1), &
        & LAM, R, Rinv)
    do l = 1, Nv ! index for char-variables
      W = 0.0_8
      do j = 1, Ns ! index for stencil-size
        do i = 1, Nv ! index for cons-variables
          W(j) = W(j) + Rinv(l,i)*Q(i,j)
        end do
      end do
      call MP5_sub(W(1), W(2), W(3), W(4), W(5), Wc(l))
    end do
    Qc = 0.0_8
    do j = 1, Nv
      do i = 1, Nv
        Qc(i) = Qc(i) + R(i,j)*Wc(j)
      end do
    end do
    call calc_basic(Qc, SA, rhoL, uL, vL, eL)
    ! QR
    call calc_eigen(rho_s(3), u_s(3), v_s(3), e_s(3), ix_s(2), iy_s(2), &
        & LAM, R, Rinv)
    do l = 1, Nv ! index for char-variables
      W = 0.0_8
      do j = 1, Ns ! index for stencil-size
        do i = 1, Nv ! index for cons-variables
          W(j) = W(j) + Rinv(l,i)*Q(i,j)
        end do
      end do
      call MP5_sub(W(6), W(5), W(4), W(3), W(2), Wc(l))
    end do
    Qc = 0.0_8
    do j = 1, Nv
      do i = 1, Nv
        Qc(i) = Qc(i) + R(i,j)*Wc(j)
      end do
    end do
    call calc_basic(Qc, SA, rhoR, uR, vR, eR)
  end subroutine reconstruction_char_MP5
  !======================================================================
  !-----------------------reconstruction by flux-------------------------
  !======================================================================
  subroutine reconstruction_flux_MUSCL3(rho_s, u_s, v_s, e_s, S_s, SA, &
    & ix_s, iy_s, &
    & rhoL, uL, vL, eL, rhoR, uR, vR, eR)
    implicit none
    integer, parameter   :: Nv = 4, Ns = 4
    real(8), intent(in)  :: rho_s(Ns), u_s(Ns), v_s(Ns), e_s(Ns), S_s(Ns)
    real(8), intent(in)  :: SA, ix_s(2), iy_s(2)
    real(8), intent(out) :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8)              :: Q(Nv,Ns), LAM(Nv), R(Nv,Nv), Rinv(Nv,Nv)
    real(8)              :: F(Ns), Fc(Nv), tmp, Qc(Nv)
    integer              :: i, j, l, m
    do i = 1, Ns
      call calc_conservative(rho_s(i), u_s(i), v_s(i), e_s(i), S_s(i), Q(:,i))
    end do
    ! QL
    call calc_eigen(rho_s(2), u_s(2), v_s(2), e_s(2), ix_s(1), iy_s(1), &
        & LAM, R, Rinv)
    do l = 1, Nv ! index for flux components
      F = 0.0_8
      do m = 1, Ns ! index for stencil-size
        do j = 1, Nv ! index for cons-variables
          do i = 1, Nv
            F(m) = F(m) + R(l,i)*LAM(i)*Rinv(i,j)*Q(j,m)
          end do
        end do
      end do
      call MUSCL3(F(1), F(2), F(3), F(4), Fc(l), tmp)
    end do
    Qc = 0.0_8
    do l = 1, Nv
      do j = 1, Nv
        do i = 1, Nv
          Qc(i) = Qc(i) + R(i,j)/LAM(j)*Rinv(j,l)*Fc(l)
        end do
      end do
    end do
    call calc_basic(Qc, SA, rhoL, uL, vL, eL)
    ! QR
    call calc_eigen(rho_s(3), u_s(3), v_s(3), e_s(3), ix_s(2), iy_s(2), &
        & LAM, R, Rinv)
    do l = 1, Nv ! index for char-variables
      F = 0.0_8
      do m = 1, Ns ! index for stencil-size
        do j = 1, Nv ! index for cons-variables
          do i = 1, Nv
            F(m) = F(m) + R(l,i)*LAM(i)*Rinv(i,j)*Q(j,m)
          end do
        end do
      end do
      call MUSCL3(F(1), F(2), F(3), F(4), tmp, Fc(l))
    end do
    Qc = 0.0_8
    do l = 1, Nv
      do j = 1, Nv
        do i = 1, Nv
          Qc(i) = Qc(i) + R(i,j)/LAM(j)*Rinv(j,l)*Fc(l)
        end do
      end do
    end do
    call calc_basic(Qc, SA, rhoR, uR, vR, eR)
  end subroutine reconstruction_flux_MUSCL3
  !======================================================================
  !---------------calculation of numerical diffusion flux----------------
  !======================================================================
  subroutine construct_numflux_diff(Re, Pr, rhoL, uL, vL, eL, &
              & rhoR, uR, vR, eR, ixSA, iySA, jxSA, jySA, SA, &
              & u_forigrad, v_forigrad, temp_forigrad, &
              & u_forjgrad, v_forjgrad, temp_forjgrad, &
              & flag_i, flag_j, Fd)
    implicit none
    real(8), intent(in)  :: Re, Pr
    real(8), intent(in)  :: rhoL, uL, vL, eL, rhoR, uR, vR, eR
    real(8), intent(in)  :: ixSA, iySA, jxSA, jySA, SA
    real(8), intent(in)  :: u_forigrad(4), v_forigrad(4), temp_forigrad(4)
    real(8), intent(in)  :: u_forjgrad(2,5),v_forjgrad(2,5),temp_forjgrad(2,5)
    real(8), intent(in)  :: flag_i, flag_j
    real(8), intent(inout) :: Fd(4)
    real(8)              :: ui, vi, Ti, qijm, qijp, uj, vj, Tj
    real(8)              :: rhoA, uA, vA, eA
    ! i-gradient for stress tensor
        ! calculating by 3rd order central difference
        ! evaluated at (i+0.5, j) from stencil [(i-1,j),...,(i+2,j)]
    call central_diff3(u_forigrad(1), u_forigrad(2), &
                      & u_forigrad(3), u_forigrad(4), ui)
    call central_diff3(v_forigrad(1), v_forigrad(2), &
                      & v_forigrad(3), v_forigrad(4), vi)
    call central_diff3(temp_forigrad(1), temp_forigrad(2), &
                      & temp_forigrad(3), temp_forigrad(4), Ti)
    ! j-gradient for stress tensor
    ! calculating by average of 4th order central difference
    ! evaluated at (i+0.5, j) by averaging (i,j) and (i+1,j),
    ! which are calculated from stencil [(i,j-2),...,(i,j+2)]
    ! and [(i+1,j-2),...,(i+1,j+2)], respectively
    call central_diff4(u_forjgrad(1,1), u_forjgrad(1,2), &
        & u_forjgrad(1,4), u_forjgrad(1,5), qijm)
    call central_diff4(u_forjgrad(2,1), u_forjgrad(2,2), &
        & u_forjgrad(2,4), u_forjgrad(2,5), qijp)
    uj     = 0.5_8*(qijm + qijp)
    call central_diff4(v_forjgrad(1,1), v_forjgrad(1,2), &
        & v_forjgrad(1,4), v_forjgrad(1,5), qijm)
    call central_diff4(v_forjgrad(2,1), v_forjgrad(2,2), &
        & v_forjgrad(2,4), v_forjgrad(2,5), qijp)
    vj     = 0.5_8*(qijm + qijp)
    call central_diff4(temp_forjgrad(1,1), temp_forjgrad(1,2), &
        & temp_forjgrad(1,4), temp_forjgrad(1,5), qijm)
    call central_diff4(temp_forjgrad(2,1), temp_forjgrad(2,2), &
        & temp_forjgrad(2,4), temp_forjgrad(2,5), qijp)
    Tj     = 0.5_8*(qijm + qijp)
    ! calculating diffusion flux
    call Roe_average(rhoL, uL, vL, eL, rhoR, uR, vR, eR, &
        & rhoA, uA, vA, eA)
    call calc_flux_diff(Re, Pr, uA, vA, &
          & ixSA*flag_i + jxSA*flag_j, iySA*flag_i + jySA*flag_j, &
          & ixSA/SA, iySA/SA, jxSA/SA, jySA/SA, &
          & ui*flag_i + uj*flag_j, uj*flag_i + ui*flag_j, &
          & vi*flag_i + vj*flag_j, vj*flag_i + vi*flag_j, &
          & Ti*flag_i + Tj*flag_j, Tj*flag_i + Ti*flag_j, Fd)
  end subroutine construct_numflux_diff
end module proceed