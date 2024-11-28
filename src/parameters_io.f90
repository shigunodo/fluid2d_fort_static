module parameters_io
  !$ use omp_lib
  implicit none
  !=========================COMMON PARAMETERS============================
  !------------------parameters for coordinate---------------------------
  ! coordinate is given as (Ni, Nj) arrays
  ! but Nb columns from ends are dummy/bc for calculating diffs
  ! so evaluation of eqs is conducted on (Ni-2*Nb, Nj-2*Nb) arrays
  !----------------------------------------------------------------------
  ! total grid number and dummy grids for boundary conditions
  ! Nb >= 3 for MUSCL3, Nb >= 4 for MP5
  integer, parameter :: Ni = 408, Nj = 408, Nb = 4
  ! coordinate: (Ni, Nj) array
  real(8), protected, save :: x(Ni, Nj), y(Ni, Nj)
  ! dx for CFL condition:evaluated on the points (i, j)
  ! metrix              :evaluated on the points (i, j)
  ! inverse of Jacobian :evaluated on the points (i+0.5, j+0.5)
  real(8), protected, save :: dx(Ni-2*Nb, Nj-2*Nb)
  real(8), protected, save :: ixS(Ni-2*Nb+2, Nj-2*Nb+2), iyS(Ni-2*Nb+2, Nj-2*Nb+2), &
                & jxS(Ni-2*Nb+2, Nj-2*Nb+2), jyS(Ni-2*Nb+2, Nj-2*Nb+2)
  real(8), protected, save :: S(Ni-1, Nj-1)
  !----------------------dimensionless numbers---------------------------
  ! details are refffered to README.pdf
  !Alp = p_inf / (rho_inf * U_inf * U_inf)
  real(8), parameter :: Re = 1.0d+3, Pr=1.0_8
  !real(8), parameter :: Alp = 7.936507936507936_8
  !----------------------directory for io--------------------------------
  !!!!!!!!!ワーキングディレクトリを設定してください!!!!!!!!!!!!
  character(len=256), parameter, private :: dir = &
          & '../data/'
  !----------------------parameters for time integration-----------------
  real(8), parameter :: Tmax = 3.0_8 ! integration time
  integer, parameter :: Nout = 100     ! output times
  real(8), parameter :: Tout = Tmax / real(Nout, 8) ! output interval
  !-----------------------parameters for boundary condition--------------
  ! for cyinder bc
  integer, protected :: Nc1, Nc2
  !-----------------------scientific parameters--------------------------
  real(8), parameter :: pi = 3.14159265358979323846264338327950288_8
  !
  !
  !
  !
  !==============================SUBROUTINES=============================
  ! subroutines for initializing coordinate and I/O from/into files
  ! private setting for subroutines
  private :: input_coordinate, set_dx
  !======================================================================
  contains
  !====================================================================
  !------------------outputting basic variables------------------------
  ! outputting (rho, u, v, e) into one file
  !====================================================================
  subroutine output_basic(tstep, t, iter, rho, u, v, e, cpu_begin)
    implicit none
    character(len=256)  :: name, line
    integer             :: idf, hrs, mins, secs, hrs_r, mins_r, secs_r
    real(8)             :: cpu_now, time_r
    integer, intent(in) :: tstep, iter
    real(8), intent(in) :: t, cpu_begin
    real(8), intent(in) :: rho(:,:), u(:,:), v(:,:), e(:,:)
    ! make str for the name of the output file
    write(name, '("b",i7.7,".dat")') tstep
    name = trim(dir) // trim(name)
    ! open file
    open(newunit=idf, file=name,form='formatted')
    ! write into the file
    write(idf,*) rho, u, v, e
    ! close the file
    close(idf)
    ! mesure the CPU time
    call cpu_time(cpu_now)
    !$ cpu_now = omp_get_wtime()
    cpu_now = cpu_now - cpu_begin
    ! convert it into hours and minutes
    call time_to_hms(cpu_now, hrs, mins, secs)
    ! rest time estimate
    time_r = cpu_now/real(tstep,8)*real(Nout-tstep,8)
    call time_to_hms(time_r, hrs_r, mins_r, secs_r)
    ! write the status
    write(line, '("elapsed: ",i3," h ",i2.2," m ",i2.2, &
    & " s | tstep = ",i5," | t = ",f10.4," | iteration = ",i7, &
    & " | rest: ",i3," h ",i2.2," m ",i2.2," s")') &
    & hrs, mins, secs, tstep, t, iter, hrs_r, mins_r, secs_r
    ! output status to standard
    write(*,*) trim(line)
    ! open settings.tda and writing status
    name = trim(dir) // "settings.dat"
    open(newunit=idf, file=name, form='formatted', &
            & action='write', status="old", position='append')
    write(idf,*) trim(line)
    close(idf)
  end subroutine output_basic
  !----------------------convert time into h/mm/ss---------------------
  subroutine time_to_hms(time, hrs, mins, secs)
    implicit none
    real(8), intent(in)  :: time
    integer, intent(out) :: hrs, mins, secs
    real(8)              :: tmp
    hrs    = int(time/3600.0_8)
    tmp = time - 3600.0_8 * real(hrs,8)
    mins    = int(tmp/60.0_8)
    tmp = tmp - 60.0_8 * real(mins,8)
    secs    = int(tmp)
  end subroutine
  !====================================================================
  !--------------------------output settings---------------------------
  !====================================================================
  subroutine output_settings()
    implicit none
    integer            :: idf
    character(len=256) :: name
    ! make str for the name of the output file
    name = trim(dir) // "settings.dat"
    ! open the file
    open(newunit=idf, file=name,form='formatted')
    ! output
    write(idf,*) "Ni= ", Ni
    write(idf,*) "Nj= ", Nj
    write(idf,*) "Nb= ", Nb
    write(idf,*) "Tmax= ", Tmax
    write(idf,*) "Nout= ", Nout
    write(idf,*) "Re= ", Re
    write(idf,*) "Pr= ", Pr
    write(idf,*) ""
    ! close file
    close(idf)
  end subroutine output_settings
  !====================================================================
  !----------------------input the coordinate--------------------------
  !====================================================================
  subroutine input_coordinate()
    implicit none
    integer            :: idf
    character(len=256) :: name
    ! make str for the name of the output file
    name = trim(dir) // "coordinate.dat"
    ! open the file
    open(newunit=idf, file=name,form='formatted')
    ! read the coordinate
    read(idf,*) x, y
    ! close the file
    close(idf)
  end subroutine input_coordinate
  !====================================================================
  !--------------------input initial basic variables-------------------
  !====================================================================
  subroutine input_initial_basic(rho, u, v, e)
    implicit none
    integer                :: idf
    character(len=256)     :: name
    real(8), intent(inout) :: rho(:,:), u(:,:), v(:,:), e(:,:)
    ! make str for the name of the output file
    name = trim(dir) // "b0000000.dat"
    ! open the file
    open(newunit=idf, file=name,form='formatted')
    ! read the coordinate
    read(idf,*) rho, u, v, e
    ! close the file
    close(idf)
  end subroutine input_initial_basic
  !----------------------for DEBUG--------------------------------------
  subroutine output_coordinate()
    implicit none
    integer :: idf, i, j
    real(8) :: r(Ni, Nj), th(Ni, Nj)
    r = sqrt(x*x+y*y)
    th = acos(x/r)
    open(newunit=idf, file="coordinate_r_test.txt")
    do j = 1, Nj
      write(idf,*) (r(i, j), i=1,Ni)
    end do
    close(idf)
    open(newunit=idf, file="coordinate_th_test.txt")
    do j = 1, Nj
      write(idf,*) (th(i,j), i=1,Ni)
    end do
    close(idf)
  end subroutine output_coordinate
  !
  !
  !
  !
  !====================================================================
  !---------------------setting the coordinate-------------------------
  ! inputting the coordinate
  ! calculating metrix and Jacobian
  ! calculating dx for CFL condition
  !====================================================================
  subroutine set_coordinate()
    implicit none
    integer :: i, j
    ! input the coordinate
    call input_coordinate()
    ! calculate the metrix
    ! the values are evaluated at the points (i, j)
    ! returning (Ni-2*Nb+2,Nj-2*Nb+2) arrays
    ! using 2nd order central difference
    ixS =  0.5_8 * (y(Nb:Ni-Nb+1, Nb+1:Nj-Nb+2) &
                & - y(Nb:Ni-Nb+1, Nb-1:Nj-Nb))
    iyS = -0.5_8 * (x(Nb:Ni-Nb+1, Nb+1:Nj-Nb+2) &
                & - x(Nb:Ni-Nb+1, Nb-1:Nj-Nb))
    jxS = -0.5_8 * (y(Nb+1:Ni-Nb+2, Nb:Nj-Nb+1) &
                & - y(Nb-1:Ni-Nb, Nb:Nj-Nb+1))
    jyS =  0.5_8 * (x(Nb+1:Ni-Nb+2, Nb:Nj-Nb+1) &
                & - x(Nb-1:Ni-Nb, Nb:Nj-Nb+1))
    ! calculatw the inverse of Jacobian
    ! the values are evaluated ar the points (i+0.5, j+0.5)
    ! returning (Ni-1, Nj-1) array
    S = (x(2:Ni, 2:Nj) - x(1:Ni-1, 1:Nj-1)) &
    & * (y(1:Ni-1, 2:Nj) - y(2:Ni, 1:Nj-1))
    S = S - (y(2:Ni, 2:Nj) - y(1:Ni-1, 1:Nj-1)) &
    & * (x(1:Ni-1, 2:Nj) - x(2:Ni, 1:Nj-1))
    S = S * 0.5_8
    ! calculating dx for CFL condition
    ! evaluated at the points (i, j)
    ! returning (Ni-2*Nb, Nj-2*Nb) array
    call set_dx()
    ! for cylinder boundary condition
    ! j = Nc1, Nc2 are the points where theta=pi/2,3pi/2
    Nc1 = 0
    do j = 1, int(Nj/2)
      if( x(Ni,j)>0.0_8 ) then
        Nc1 = Nc1+1
      end if
    end do
    Nc2 = Nj+1
    do j = Nj, int(Nj/2), -1
      if( x(Ni,j)>0.0_8 ) then
        Nc2 = Nc2-1
      end if
    end do
  end subroutine set_coordinate
  !======================================================================
  !---------------calculating dx for CFL condition-----------------------
  ! calculate the minimum of the distances between adjacent 4 grids
  ! returning (Ni-2*Nb, Nj-2*Nb) array
  !======================================================================
  subroutine set_dx()
    implicit none
    real(8) :: dx1, dy1, dd
    integer :: i, j

    do j = 1, Nj - 2*Nb
      do i = 1, Ni - 2*Nb
        dx1 = x(Nb+i+1, Nb+j) - x(Nb+i, Nb+j)
        dy1 = y(Nb+i+1, Nb+j) - y(Nb+i, Nb+j)
        dx(i,j) = dx1 * dx1 + dy1 * dy1
        dx1 = x(Nb+i, Nb+j+1) - x(Nb+i, Nb+j)
        dy1 = y(Nb+i, Nb+j+1) - y(Nb+i, Nb+j)
        dd  = dx1 * dx1 + dy1 * dy1
        dx(i,j) = min(dx(i,j), dd)
        dx1 = x(Nb+i-1, Nb+j) - x(Nb+i, Nb+j)
        dy1 = y(Nb+i-1, Nb+j) - y(Nb+i, Nb+j)
        dd  = dx1 * dx1 + dy1 * dy1
        dx(i,j) = min(dx(i,j), dd)
        dx1 = x(Nb+i, Nb+j-1) - x(Nb+i, Nb+j)
        dy1 = y(Nb+i, Nb+j-1) - y(Nb+i, Nb+j)
        dd  = dx1 * dx1 + dy1 * dy1
        dx(i,j) = min(dx(i,j), dd)
        dx(i,j) = sqrt(dx(i,j))
      end do
    end do
  end subroutine set_dx
end module parameters_io