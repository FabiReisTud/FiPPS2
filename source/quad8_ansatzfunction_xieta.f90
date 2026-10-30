!
!  This program developed in FORTRAN is a Finite Element solver for linear-
!  static analyses as well as linearized stability analyses. It is inherently
!  coupled to the open-source panel method APAME for providing fluid-structure-
!  interaction capabilites.
!
!  Copyright (C) 2024 TUD Dresden University of Technology
!
!  This file is part of FiPPS².
!
!  FiPPS² is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  FiPPS² is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
! =================================================================================================
!
!> @brief
!> Calculates all necessary values for numerical integration of a quad8 finite
!> shell element
!
!> @details
!
!> @author   Martin Rädel, TU Dresden, Diplomarbeit, 08.12.2009
!
!> @author   Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 17.10.2012
!
!> $Id: quad8_ansatzfunction_xieta.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)

!
! Use
!

!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
!
! Data types
!
  double precision, intent(in)                                  :: xi           !< xi coordinate inside the element (-1 <= xi <= 1)
  double precision, intent(in)                                  :: eta          !< eta coordinate inside the element (-1 <= eta <= 1)

  double precision, dimension(8), optional, intent(out)         :: Ni           !< value of the shape function for each node
  double precision, dimension(8), optional, intent(out)         :: dNidxi       !< value of the first derivation of the shape function which respect to xi for each node
  double precision, dimension(8), optional, intent(out)         :: dNideta      !< value of the first derivation of the shape function which respect to eta for each node
  double precision, dimension(8), optional, intent(out)         :: d2Nidxi2     !< value of the second derivation of the shape function which respect to xi for each node
  double precision, dimension(8), optional, intent(out)         :: d2Nidxideta  !< value of the second derivation of the shape function which respect to xi and eta for each node
  double precision, dimension(8), optional, intent(out)         :: d2Nidetadxi  !< value of the second derivation of the shape function which respect to eta and ci for each node
  double precision, dimension(8), optional, intent(out)         :: d2Nideta2    !< value of the second derivation of the shape function which respect to eta for each node

  integer                                                       :: err_code=0
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation
!
  ! calculate shape-function values with respect to natural coordinates

  if (present(Ni) .eq. .true.) then
    ! corner nodes

    Ni(1) = 0.25D0*(xi**2-xi**2*eta+xi*eta-xi*eta**2+eta**2-1.D0)
    Ni(2) = 0.25D0*(xi**2-xi**2*eta-xi*eta+xi*eta**2+eta**2-1.D0)
    Ni(3) = 0.25D0*(xi**2+xi**2*eta+xi*eta+xi*eta**2+eta**2-1.D0)
    Ni(4) = 0.25D0*(xi**2+xi**2*eta-xi*eta-xi*eta**2+eta**2-1.D0)

    ! mid-side nodes

    Ni(5) = 0.5D0*( xi**2*(-1.D0+eta)-eta+1.D0)
    Ni(6) = 0.5D0*(eta**2*(-1.D0- xi)+ xi+1.D0)
    Ni(7) = 0.5D0*( xi**2*(-1.D0-eta)+eta+1.D0)
    Ni(8) = 0.5D0*(eta**2*(-1.D0+ xi)- xi+1.D0)

  end if


  ! calculate shape-function derivatives with respect to natural coordinates

  if (present(dNidxi) .eq. .true.) then

    dNidxi(1) = 0.25D0*(2.D0*xi*(1.D0-eta)+eta-eta**2)
    dNidxi(2) = 0.25D0*(2.D0*xi*(1.D0-eta)-eta+eta**2)
    dNidxi(3) = 0.25D0*(2.D0*xi*(1.D0+eta)+eta+eta**2)
    dNidxi(4) = 0.25D0*(2.D0*xi*(1.D0+eta)-eta-eta**2)

    dNidxi(5) = xi*( eta-1.D0)
    dNidxi(6) = 0.5D0*( 1.D0-eta**2)
    dNidxi(7) = xi*(-eta-1.D0)
    dNidxi(8) = 0.5D0*(-1.D0+eta**2)

  end if

  if (present(dNideta) .eq. .true.) then

    dNideta(1) = 0.25D0*(-xi**2+xi+2.D0*eta*(1.D0-xi))
    dNideta(2) = 0.25D0*(-xi**2-xi+2.D0*eta*(1.D0+xi))
    dNideta(3) = 0.25D0*( xi**2+xi+2.D0*eta*(1.D0+xi))
    dNideta(4) = 0.25D0*( xi**2-xi+2.D0*eta*(1.D0-xi))

    dNideta(5) = 0.5D0*(-1.D0+xi**2)
    dNideta(6) = eta*(-xi-1.D0)
    dNideta(7) = 0.5D0*( 1.D0-xi**2)
    dNideta(8) = eta*( xi-1.D0)

  end if


  ! calculate shape-function second derivatives with respect to natural coordinates

  if (present(d2Nidxi2) .eq. .true.) then

    d2Nidxi2(1) = 0.5D0*(1.D0-eta)
    d2Nidxi2(2) = 0.5D0*(1.D0-eta)
    d2Nidxi2(3) = 0.5D0*(1.D0+eta)
    d2Nidxi2(4) = 0.5D0*(1.D0+eta)

    d2Nidxi2(5) =  eta-1.D0
    d2Nidxi2(6) = 0.0D0
    d2Nidxi2(7) = -eta-1.D0
    d2Nidxi2(8) = 0.0D0

  end if

  if (present(d2Nidxideta) .eq. .true.) then

    d2Nidxideta(1) = 0.25D0*(-2.D0*xi+1.D0-2.D0*eta)
    d2Nidxideta(2) = 0.25D0*(-2.D0*xi-1.D0+2.D0*eta)
    d2Nidxideta(3) = 0.25D0*(+2.D0*xi+1.D0+2.D0*eta)
    d2Nidxideta(4) = 0.25D0*(+2.D0*xi-1.D0-2.D0*eta)

    d2Nidxideta(5) =  xi
    d2Nidxideta(6) = -eta
    d2Nidxideta(7) = -xi
    d2Nidxideta(8) =  eta

  end if

  if (present(d2Nidetadxi) .eq. .true.) then

    d2Nidetadxi(1) = 0.25D0*(-2.D0*xi+1.D0-2.D0*eta)
    d2Nidetadxi(2) = 0.25D0*(-2.D0*xi-1.D0+2.D0*eta)
    d2Nidetadxi(3) = 0.25D0*( 2.D0*xi+1.D0+2.D0*eta)
    d2Nidetadxi(4) = 0.25D0*( 2.D0*xi-1.D0-2.D0*eta)

    d2Nidetadxi(5) =  xi
    d2Nidetadxi(6) = -eta
    d2Nidetadxi(7) = -xi
    d2Nidetadxi(8) =  eta

  end if

  if (present(d2Nideta2) .eq. .true.) then

    d2Nideta2(1) = 0.5D0*(1.D0-xi)
    d2Nideta2(2) = 0.5D0*(1.D0+xi)
    d2Nideta2(3) = 0.5D0*(1.D0+xi)
    d2Nideta2(4) = 0.5D0*(1.D0-xi)

    d2Nideta2(5) = 0.D0
    d2Nideta2(6) = -xi-1.D0
    d2Nideta2(7) = 0.D0
    d2Nideta2(8) =  xi-1.D0

  end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_ansatzfunction_xieta'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

end subroutine quad8_ansatzfunction_xieta
