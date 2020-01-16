module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => rm1d_init_one_grid

    call set_coordinate_system("Cartesian")
    call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    ! iprob==1 Sod test- rarefaction and shock
    if (iprob==1) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 1.0d0
        elsewhere
           w(ix^S,rho_)   = 1.25d-1
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 0.1d0
        end where
    ! iprob==2 Lax problem
    else if (iprob==2) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 4.45d-1
           w(ix^S,mom(1)) = 6.98d-1
           w(ix^S,e_)     = 3.53d0
        elsewhere
           w(ix^S,rho_)   = 5.0d-1
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 5.71d-1
        end where
    ! iprob==3  weaker version of the 123 problem
    else if (iprob==3) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = -0.5d0
           w(ix^S,e_)     = 1.0d0
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.5d0
           w(ix^S,e_)     = 1.0d0
        end where
    ! iprob==4 blast wave- right shock, contact and left rarefaction
    else if (iprob==4) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 1.0d3
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 0.1d-1
        end where
    ! iprob==5 reverse blast wave
    else if (iprob==5) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 1.0d-2
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 1.0d2
        end where
    ! iprob==6 colliding shocks
    else if (iprob==6) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 5.9924d0
           w(ix^S,mom(1)) = 19.5975d0
           w(ix^S,e_)     = 460.894d0
        elsewhere
           w(ix^S,rho_)   = 5.99242d0
           w(ix^S,mom(1)) = -6.19633d0
           w(ix^S,e_)     = 46.095d0
        end where
    else
        call mpistop("iprob not available!")
    end if

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine rm1d_init_one_grid

end module mod_usr
