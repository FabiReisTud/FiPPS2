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
!> control subroutine for computing element conduction matrix for quadratic
!> 20-node solid element and assemble it to the total
!> stiffness matrix
!
!> @details
!> subroutine calculates the element conduction matrix for quadratic 20-node
!> solid element and assembles it to the total stiffness matrix
!
!> @author Fabian Reisener, TU Dresden, WM 2026
!
! =================================================================================================
subroutine sol3_tsolid20_cond(fesim,node_coords,node_ids,lsolid20Num,elemNum,k,tsolid20_C)
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
! Interfaces
! =================================================================================================
!
! Data types
!
! Input
type(fe_simulation), intent(in)                   :: fesim
double precision, intent(in)                      :: node_coords(20,3)  !< (20,3)-Array containing global coordinates of element nodes
integer, intent(in)                               :: node_ids(20)       !< 20-Array containing all node-ids for eval of IP-temperature
integer, intent(in)                               :: lsolid20Num        !< number of the lsolid20 element (only used for error messages)
integer, intent(in)                               :: elemNum            !< element number of the lsolid20 element (only used for error messages)
type(thermal_curve_type), intent(in)              :: k                  !< thermal curve of element
!
! Output
double precision, intent(out), dimension(20,20)   :: tsolid20_C !< 20x20-Element conduction matrix
!
! Internal
integer                                           :: err_code=0, ii, jj, kk, N
double precision, allocatable, dimension(:)       :: xi_ip, eta_ip, zeta_ip
double precision, allocatable, dimension(:)       :: w_xi, w_eta, w_zeta
double precision, dimension(20)                   :: dNidx, dNidy, dNidz, J
double precision, dimension(20)                   :: Ni,dNidxi, dNideta, dNidzeta
double precision                                  :: weight, detJac


double precision                                  :: c, T_ip
!
! =================================================================================================
!
! Initialisation
!
tsolid20_C = 0.d0
!
! =================================================================================================
!
! Calculation
!
! Total number of integration points in Tsolid20-Element
N = 8

! allocate layer specific arrays for integration point values
allocate(xi_ip(N),eta_ip(N),zeta_ip(N),w_xi(N),w_eta(N),w_zeta(N))

! set integration points and corresponding weights
call integration_points_3d( 2, "Gauss",   xi_ip, w_xi, &
                            2, "Gauss",  eta_ip, w_eta, &
                            2, "Gauss", zeta_ip, w_zeta)

! Loop over integration points
do ii=1,N
    ! Get the weight of the integration point
    weight = w_xi(ii)*w_eta(ii)*w_zeta(ii)

    ! Temperatur an IP über Quadratur holen, falls ein Temperaturfeld bereits berechnet wurde
    if (allocated(fesim%thermal%temperature)) then
        T_ip = 0.d0
        call hexa20_shapefunctions_xi_eta_zeta(xi_ip(ii), eta_ip(ii), zeta_ip(ii), Ni, dNidxi, dNideta, dNidzeta)
        do jj=1,20
            if (node_ids(jj) >= 1 .and. node_ids(jj) <= size(fesim%thermal%temperature)) then
                T_ip = T_ip + Ni(jj) * fesim%thermal%temperature(node_ids(jj)) ! Knotenanteil der IP-Temperatur hinzufügen
            end if
        end do
    else
        T_ip = fesim%thermal%tref ! IP-Temperatur wird ansonsten auf Referenztemperatur der Simulation gesetzt
    end if
    ! Wärmeleitungskoeffizient über Temperatur am IP holen
    c = thermal_curve_eval(k,T_ip)

    ! This adds the conduction contributions of all nodes
    ! Get derivates of shape functions with respect to global coordinates
    call hexa20_shapefunctions_x_y_z(xi_ip(ii),eta_ip(ii),zeta_ip(ii),node_coords,dNidx,dNidy,dNidz,detJac,J)
    ! check if Element is malformed
    if (detJac < 0.d0) then
       write(*,*) 'for lsolid20 number:                  ', lsolid20Num
       write(*,*) 'element number:                       ', elemNum
       write(*,*) 'negative determinant of the jacobian: ', detJac
       err_code = 2
       goto 9999
    end if

    ! Loop over all nodes
    do kk = 1,20   ! node associated with row block
        do jj = 1,kk ! node associated with column block, only until kk, as upper third of conduction matrix is later simply mirrored
            ! Hinzufügen des Knotenanteils zur Leitungsmatrix
            tsolid20_C(jj,kk) = tsolid20_C(jj,kk) + c*weight*dNidx(jj)*dNidx(kk) + c*weight*dNidy(jj)*dNidy(kk) + c*weight*dNidz(jj)*dNidz(kk)
            if (abs(tsolid20_C(jj,kk))<1e-10) then
                tsolid20_C(jj,kk) = 0.d0 ! Remove veeeery small values (<1e-10) and
            end if
        end do
    end do
end do

! mirror matrix
do kk = 1, 20
    do jj = 1, kk-1
        tsolid20_C(kk,jj) = tsolid20_C(jj,kk)
    end do
end do

! cleanup
deallocate(xi_ip,eta_ip,zeta_ip,w_xi,w_eta,w_zeta)
! =================================================================================================
!
! Error handling
!
9999 continue
!
if (err_code /= 0) then
!
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'tsolid20_cond'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!
end if
!
return
!
end subroutine
