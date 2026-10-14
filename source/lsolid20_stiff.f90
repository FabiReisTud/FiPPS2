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
!> Calculates element stiffness matrix for layered 20-node solid element Lsolid20
!
!> @details
!> Calculation of element stiffness matrix for layered 20-node solid element
!> Lsolid20; integration is done layerwise, zeta = [-1,1] is substituted for
!> each layer to run over thickness of layer such that the same integration point
!> values can be used in every layer; by this, the layer thicknesses are taken into
!> account as relative to the total laminat thickness which therefore needs not to be
!> equal to the actual thickness of the element as layer thicknesses are scaled to it;
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_stiff.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_stiff(fesim,nlay,Clay,lth,tth,node_coords,nip,nop,cid,lsolid20Num,elemNum,lsolid20_K)
! =================================================================================================
! use
!
 use integration_schemes
 use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation), intent(in)                   :: fesim
integer, intent(in)                               :: nlay !< Total number of layers in current element
integer, intent(in)                               :: nip !< Number of integration points in plane per direction
integer, intent(in)                               :: nop(nlay) !< (nlay)-Array with number of integration points over layer thickness for each layer
integer, intent(in)                               :: cid !< ID of coordinate system for orientation of element coosy
double precision, intent(in)                      :: Clay(6,6,nlay) !< (6,6,nlay)-Array which contains the material stiffness matrix in the element coordinate system for each layer
double precision, intent(in)                      :: lth(nlay) !< (nlay)-Array which contains the thickness for each layer
double precision, intent(in)                      :: tth !< Total thickness of laminat (sum of all values in lth)
double precision, intent(in)                      :: node_coords(20,3) !< (20,3)-Array containing global coordinates of element nodes
integer, intent(in)                               :: lsolid20Num           !< number of the lsolid20 element (only used for error messages)
integer, intent(in)                               :: elemNum            !< element number of the lsolid20 element (only used for error messages)
!
! Output
!
double precision, intent(out), dimension(60,60)   :: lsolid20_K !< 60x60-Element stiffness matrix (no values for rotational dofs)
!
! Internal
!
integer                                           :: err_code=0, N, ii
integer                                           :: jj, kk, js, ks
integer                                           :: mm, nn
integer                                           :: lay
double precision                                  :: fac
double precision                                  :: Bxj, Byj, Bzj, Bxk, Byk, Bzk
double precision, allocatable, dimension(:)       :: xi_ip, eta_ip, zeta_ip
double precision, allocatable, dimension(:)       :: w_xi, w_eta, w_zeta
double precision, dimension(20)                   :: dNidx, dNidy, dNidz
double precision, dimension(3,3)                  :: Kpart
double precision                                  :: weight,detJac,J(3,3),TC(6,6)
double precision                                  :: zeta_elem, lth_sum
! double precision                                  :: first_detJac
! logical                                           :: sign_switched
!
! =================================================================================================
!
! Initialisation
!
 lsolid20_K = 0.d0
 lth_sum = 0.d0
!  sign_switched = .false.
!
! =================================================================================================
!
! Calculation
!
 do lay = 1,nlay

   ! Total number of integration points in current layer
   N = nip*nip*nop(lay)

   ! allocate layer specific arrays for integration point values
   allocate(xi_ip(N),eta_ip(N),zeta_ip(N),w_xi(N),w_eta(N),w_zeta(N))

   ! get integration points and corresponding weights
   call integration_points_3d(     nip, "Gauss",   xi_ip, w_xi, &
                                   nip, "Gauss",  eta_ip, w_eta, &
                              nop(lay), "Simps", zeta_ip, w_zeta)

   ! calculate summed thickness until current layer
   lth_sum = lth_sum + lth(lay)

   ! Loop over integration points
   do ii = 1,N

     weight = w_xi(ii)*w_eta(ii)*w_zeta(ii)

     ! transform layerwise zeta-value of integration point zeta_ip(ii)
     ! to elementwise zeta_value zeta_elem
     ! (see eq.(22) in Panda, S.; Natarajan, R.: Analysis of laminated composite shell
     !  structures by finite element method. Computers & Structures 14(3-4). 1981. pp. 225-230)
     zeta_elem = -1.d0 + (2.d0*lth_sum - lth(lay)*(1.d0-zeta_ip(ii)))/tth

     ! Get derivates of shape functions with respect to global coordinates
     call hexa20_shapefunctions_x_y_z(xi_ip(ii),eta_ip(ii),zeta_elem,node_coords,dNidx,dNidy,dNidz,detJac,J)

!      if (ii .eq. 1) then
!        first_detJac = detJac
!      end if

!       ! Check i sign of Jacobian has changed
!      if (detJac*first_detJac .lt. 0.d0) then
!        sign_switched = .true.
!      end if

     if (detJac < 0.d0) then
       write(*,*) 'for lsolid20 number:                  ', lsolid20Num
       write(*,*) 'element number:                       ', elemNum
       write(*,*) 'negative determinant of the jacobian: ', detJac
       err_code = 2
       goto 9999
     end if

!      if (sign_switched .eq. .true.) then
!        write(*,*) 'for lsolid20 number:                  ', lsolid20Num
!        write(*,*) 'element number:                       ', elemNum
!        write(*,*) 'sign switch in determinant of the jacobian'
!      end if

     ! Transform material stiffness matrix from element
     ! to global coordinate system
     call transform_eg_cmat(fesim, Clay(:,:,lay), J, cid, TC)
     ! This loop adds the stiffness contributions of all nodes
     do kk = 1,20   ! node associated with row block
       do jj = 1,kk ! node associated with column block, only lower third of stiffness computed and later mirrored
         ! Ableitung der Interpolationsfunkten des jj-ten Knoten
         Bxj = dNidx(jj)
         Byj = dNidy(jj)
         Bzj = dNidz(jj)
         ! Ableitung der Interpolationsfunkten des kk-ten Knoten
         Bxk = dNidx(kk)
         Byk = dNidy(kk)
         Bzk = dNidz(kk)
         ! Anteile des Knotens
         Kpart( 1, 1) = Bxj*(Bxk*TC(1, 1) + Byk*TC(6, 1) + Bzk*TC(5, 1)) + Byj*(Bxk*TC(1, 6) + Byk*TC(6, 6) + Bzk*TC(5, 6)) + Bzj*(Bxk*TC(1, 5) + Byk*TC(6, 5) + Bzk*TC(5, 5))
         Kpart( 2, 1) = Bxj*(Bxk*TC(6, 1) + Byk*TC(2, 1) + Bzk*TC(4, 1)) + Byj*(Bxk*TC(6, 6) + Byk*TC(2, 6) + Bzk*TC(4, 6)) + Bzj*(Bxk*TC(6, 5) + Byk*TC(2, 5) + Bzk*TC(4, 5))
         Kpart( 3, 1) = Bxj*(Bxk*TC(5, 1) + Byk*TC(4, 1) + Bzk*TC(3, 1)) + Byj*(Bxk*TC(5, 6) + Byk*TC(4, 6) + Bzk*TC(3, 6)) + Bzj*(Bxk*TC(5, 5) + Byk*TC(4, 5) + Bzk*TC(3, 5))
         Kpart( 1, 2) = Bxj*(Bxk*TC(1, 6) + Byk*TC(6, 6) + Bzk*TC(5, 6)) + Byj*(Bxk*TC(1, 2) + Byk*TC(6, 2) + Bzk*TC(5, 2)) + Bzj*(Bxk*TC(1, 4) + Byk*TC(6, 4) + Bzk*TC(5, 4))
         Kpart( 2, 2) = Bxj*(Bxk*TC(6, 6) + Byk*TC(2, 6) + Bzk*TC(4, 6)) + Byj*(Bxk*TC(6, 2) + Byk*TC(2, 2) + Bzk*TC(4, 2)) + Bzj*(Bxk*TC(6, 4) + Byk*TC(2, 4) + Bzk*TC(4, 4))
         Kpart( 3, 2) = Bxj*(Bxk*TC(5, 6) + Byk*TC(4, 6) + Bzk*TC(3, 6)) + Byj*(Bxk*TC(5, 2) + Byk*TC(4, 2) + Bzk*TC(3, 2)) + Bzj*(Bxk*TC(5, 4) + Byk*TC(4, 4) + Bzk*TC(3, 4))
         Kpart( 1, 3) = Bxj*(Bxk*TC(1, 5) + Byk*TC(6, 5) + Bzk*TC(5, 5)) + Byj*(Bxk*TC(1, 4) + Byk*TC(6, 4) + Bzk*TC(5, 4)) + Bzj*(Bxk*TC(1, 3) + Byk*TC(6, 3) + Bzk*TC(5, 3))
         Kpart( 2, 3) = Bxj*(Bxk*TC(6, 5) + Byk*TC(2, 5) + Bzk*TC(4, 5)) + Byj*(Bxk*TC(6, 4) + Byk*TC(2, 4) + Bzk*TC(4, 4)) + Bzj*(Bxk*TC(6, 3) + Byk*TC(2, 3) + Bzk*TC(4, 3))
         Kpart( 3, 3) = Bxj*(Bxk*TC(5, 5) + Byk*TC(4, 5) + Bzk*TC(3, 5)) + Byj*(Bxk*TC(5, 4) + Byk*TC(4, 4) + Bzk*TC(3, 4)) + Bzj*(Bxk*TC(5, 3) + Byk*TC(4, 3) + Bzk*TC(3, 3))
         ! Factor in Abhängigkeit von dem Knotengewicht und dem Anteil der Laminatdicke zur Gesamtdicke
         fac = weight*detJac*lth(lay)/tth

         ks = (kk-1)*3  ! Beginning row of 3x3 block
         js = (jj-1)*3  ! Beginning column of 3x3 block

         do nn = 1,3
           do mm = 1,3
             lsolid20_K(ks+mm,js+nn) = lsolid20_K(ks+mm,js+nn) + Kpart(mm,nn)*fac   ! Hinzufügen zur Gesamtsteifigkeit des Elements unterhalb der Diagonalen

             if (jj .NE. kk) then
               lsolid20_K(js+mm,ks+nn) = lsolid20_K(js+mm,ks+nn) + Kpart(nn,mm)*fac ! Hinzufügen zur Gesamtsteifigkeit des Elements oberhalb der Diagonalen
             end if
           end do
         end do

       end do
     end do

   end do

   ! deallocate layer specific arrays for integration point values
   deallocate(xi_ip,eta_ip,zeta_ip,w_xi,w_eta,w_zeta)

 end do

!     ! If sign of all Jacobians is negative, multiply with -1.d0
!    lsolid20_K = SIGN(1.d0,first_detJac) * lsolid20_K
!
! =================================================================================================
!
! Error handling
!
9999 continue
!
if (err_code /= 0) then
!
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_stiff'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!
end if
!
return
!
end subroutine lsolid20_stiff
