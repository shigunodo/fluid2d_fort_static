! include modules
include 'parameters_io.f90'
include 'fundamentals.f90'
include 'eos.f90'
include 'bc.f90'
include 'flux_scheme.f90'
include 'proceed.f90'
!
!============================MAIN PROGRAM==============================
program main
  !$ use omp_lib
  use parameters_io
  use proceed
  implicit none
  !--------------declaration of fluid basic variables------------------
  real(8) :: rho(Ni, Nj), u(Ni, Nj), v(Ni, Nj), e(Ni, Nj)
  !--------------declaration of time variables-------------------------
  real(8) :: t = 0.0_8 ! real time
  integer :: tstep ! time steps
  integer :: iter  ! iteration in the subtroutine proceed_time
  real(8) :: cpu_begin ! CPU time at begin
  !--------------starting the mesurement of CPU time-------------------
  call cpu_time(cpu_begin)
  !$ cpu_begin = omp_get_wtime()
  !-----------setting coordinate and initialization--------------------
  call set_coordinate()
  call input_initial_basic(rho, u, v, e)
  call output_settings()
  !-----------------------time integration-----------------------------
  ! proceeding time by Tout and outputting
  do tstep = 1, Nout
    ! time integration
    call proceed_time(Tout, tstep, t, iter, rho, u, v, e)
    ! output
    call output_basic(tstep, t, iter, rho, u, v, e, cpu_begin)
  end do
  !=============================END====================================
  write(*, *) 'Main program ended.'
end program main