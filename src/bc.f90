module bc
  use parameters_io
  use ideal_eos
  implicit none
  
  contains
  !=====================================================================
  !-------------------bc for the flow around cylinder-------------------
  !=====================================================================
  subroutine bc_periodic(rho, u, v, e)
    implicit none
    real(8), intent(inout)  :: rho(Ni,Nj), u(Ni,Nj), v(Ni,Nj), e(Ni,Nj)
    integer                 :: i, j
    real(8)                 :: p, T

    ! periodic in i-direction
    rho(1:Nb,:)       = rho(Ni-Nb-Nb+1:Ni-Nb,:)
    u(1:Nb,:)         = u(Ni-Nb-Nb+1:Ni-Nb,:)
    v(1:Nb,:)         = v(Ni-Nb-Nb+1:Ni-Nb,:)
    e(1:Nb,:)         = e(Ni-Nb-Nb+1:Ni-Nb,:)
    rho(Ni-Nb+1:Ni,:) = rho(Nb+1:Nb+Nb,:)
    u(Ni-Nb+1:Ni,:)   = u(Nb+1:Nb+Nb,:)
    v(Ni-Nb+1:Ni,:)   = v(Nb+1:Nb+Nb,:)
    e(Ni-Nb+1:Ni,:)   = e(Nb+1:Nb+Nb,:)

    ! No updates for j-boundaries since Dirichlet

  end subroutine bc_periodic
end module bc