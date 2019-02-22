!=============================================================================!
!=============================================================================!
!                        Functions for HW6 Code                               !
!=============================================================================!
!=============================================================================!
!------------------------------------------------------------------------------
! 
!------------------------------------------------------------------------------
module hw6_funcs
  implicit none
  
  public :: lorentz
  public :: magfield, elecfield
  
contains
!------------------------------------------------------------------------------
  subroutine lorentz(x,v,t,dxdt,dvdt)
    implicit none
    !Passed
    real, intent(in), dimension(1:3) :: x,v !Position and velocity 
    real, intent(in) :: t                   !Time
    real, intent(out), dimension(1:3)  :: dxdt,dvdt !Time derivatives
    !Local
    real, dimension(1:3) :: B0   !Vector magnetic field
    real, dimension(1:3) :: E0   !Vector magnetic field
   
    !Get electric and magnetic fields at x and t
    E0=elecfield(x,t)
    B0=magfield(x,t)

    dxdt(1:3)=v(1:3)
    dvdt(1)= E0(1) + v(2)*B0(3) -v(3)*B0(2) 
    dvdt(2)= E0(2) + v(3)*B0(1) -v(1)*B0(3) 
    dvdt(3)= E0(3) + v(1)*B0(2) -v(2)*B0(1) 
    
  end subroutine lorentz
!------------------------------------------------------------------------------
  function magfield(x,t)  result(B)
    !  MAGNETIC FIELD Geometry choice (geom)
    !    1 Straight Field
    implicit none
    !Passed
    real, dimension(1:3) :: x   !Position
    real :: t                   !Time
    real, dimension(1:3) :: B   !Vector magnetic field
    !Local
    integer, parameter :: bgeom=1 !Geometry

    select case(bgeom)
    case(1) !Straight Magnetic Field
       B(1)=0.
       B(2)=0. 
       B(3)=1.0
    case default
       write(*,*)'ERR: Magnetic Field geometry choice bgeom unrecognized'
       stop
    end select
    
  end function magfield

!------------------------------------------------------------------------------
  function elecfield(x,t)  result(E)
    !  ELECTRIC FIELD Geometry choice (geom)
    !    1 Straight Field
    implicit none
    !Passed
    real, dimension(1:3) :: x   !Position
    real :: t                   !Time
    real, dimension(1:3) :: E   !Vector magnetic field
    !Local
    integer, parameter :: egeom=1 !Geometry

    select case(egeom)
    case(1) !Straight Electric Field
       E(1)=0.
       E(2)=0.1 
       E(3)=0.
    case default
       write(*,*)'ERR: Electric Field geometry choice egeom unrecognized'
       stop
    end select
    
  end function elecfield
!------------------------------------------------------------------------------
end module hw6_funcs
