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
!> \file sol3_set_thermal_dof_map.f90
! =================================================================================================
!
!>  \brief      subroutine for setting the thermal_dof_map
!
!>  \details    subroutine for setting the thermal_dof_map
!>              1. temperatureBC
!>              2. fluxBC
!>              3. convectionBC
!
!>  \author     Fabian Reisener, WM, 2026
!
! =================================================================================================

subroutine sol3_set_thermal_dof_map(fesim)
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscksp.h"
  use petscksp
  use fesimulation_typen
!
! =================================================================================================
!
! Data types
!
  implicit none
! Input
  type(fe_simulation),intent(inout) :: fesim
! Output
! Internal other
  integer                           :: ii,jj,ll,nn
  integer                           :: err_code = 0
  integer,dimension(8)              :: node_ids
!
! =================================================================================================
!
! Initialisation
!
!
! =================================================================================================
!
! Calculation
 ! 11. Loop over all nodes, assign 0 to their row (NID) fesim%internal%thermal_dof_map
 nn = size(fesim%knoten%nodes)
 allocate(fesim%internals%thermal_dof_map(nn,4))
 do ii=1,nn
    fesim%internals%thermal_dof_map(ii,1) = fesim%knoten%nodes(ii)%nid
    fesim%internals%thermal_dof_map(ii,2) = 0.d0 ! No BC has been assigned yet
 end do

 ! 12. Loop over all thermal boundary conditions, assign values to the row in thermal_dof_map if assigned
 ! 121. Loop over temperature BCs
 do ii=1,size(fesim%thermal%tbcs)
    ll = findloc(fesim%internals%thermal_dof_map(:,1),dble(fesim%thermal%tbcs(ii)%nid),dim=1) ! Find index of nid in thermal_dof_map
    if (ll == 0) then
        write(*,*) 'NID ',fesim%thermal%tbcs(ii)%nid,' from TemperatureBC could not be found in thermal_dof_map!'
        err_code = 2
        goto 9999
    end if

    fesim%internals%thermal_dof_map(ll,2) = 1.                                          ! Set the type to 1 to signal temperature bc
    fesim%internals%thermal_dof_map(ll,3) = fesim%thermal%tbcs(ii)%temperature          ! Set the temperature value at the node
 end do

 ! 122. Loop over flux BCs
 do ii=1,size(fesim%thermal%fluxes)
    node_ids(1) = fesim%thermal%fluxes(ii)%nid1
    node_ids(2) = fesim%thermal%fluxes(ii)%nid2
    node_ids(3) = fesim%thermal%fluxes(ii)%nid3
    node_ids(4) = fesim%thermal%fluxes(ii)%nid4
    node_ids(5) = fesim%thermal%fluxes(ii)%nid5
    node_ids(6) = fesim%thermal%fluxes(ii)%nid6
    node_ids(7) = fesim%thermal%fluxes(ii)%nid7
    node_ids(8) = fesim%thermal%fluxes(ii)%nid8

    do jj=1,8
        ll = findloc(fesim%internals%thermal_dof_map(:,1),dble(node_ids(jj)),dim=1)   ! Find index of nid in thermal_dof_map
        if (ll == 0) then
            write(*,*) 'NID ',node_ids(jj),' from fluxBc could not be found in thermal_dof_map!'
            err_code = 2
            goto 9999
        end if

        if (fesim%internals%thermal_dof_map(ll,2) .NE. 0.) then
            write(*,*) 'Cannot assign fluxBC to Node with ID ',node_ids(jj)
            if (fesim%internals%thermal_dof_map(ll,2) .EQ. 1.) then
                write(*,*) 'Since it has already been assigned a temperatureBC'
            else
                write(*,*) 'Since it has already been assigned an invalid BC-type-Nr.: ',fesim%internals%thermal_dof_map(ll,2)
            endif
            err_code = 2
            goto 9999
        end if
        fesim%internals%thermal_dof_map(ll,2) = 2.                              ! Set the type to 1 to signal temperature bc
        fesim%internals%thermal_dof_map(ll,3) = fesim%thermal%fluxes(ii)%q      ! Set the temperature value at the node
    end do
 end do

 ! 123. Loop over convection BCs
 do ii=1,size(fesim%thermal%convections)
    node_ids(1) = fesim%thermal%convections(ii)%nid1
    node_ids(2) = fesim%thermal%convections(ii)%nid2
    node_ids(3) = fesim%thermal%convections(ii)%nid3
    node_ids(4) = fesim%thermal%convections(ii)%nid4
    node_ids(5) = fesim%thermal%convections(ii)%nid5
    node_ids(6) = fesim%thermal%convections(ii)%nid6
    node_ids(7) = fesim%thermal%convections(ii)%nid7
    node_ids(8) = fesim%thermal%convections(ii)%nid8

    do jj=1,8
        ll = findloc(fesim%internals%thermal_dof_map(:,1),dble(node_ids(jj)),dim=1)   ! Find index of nid in thermal_dof_map
        if (ll == 0) then
            write(*,*) 'NID ',node_ids(jj),' from convectionBc could not be found in thermal_dof_map!'
            err_code = 2
            goto 9999
        end if
        if (fesim%internals%thermal_dof_map(ll,2) .NE. 0.) then
            write(*,*) 'Cannot assign convectionBc to Node with ID ',node_ids(jj)
            if (fesim%internals%thermal_dof_map(ll,2) .EQ. 1.) then
                write(*,*) 'Since it has already been assigned a temperatureBC'
            elseif (fesim%internals%thermal_dof_map(ll,2) .EQ. 2.) then
                write(*,*) 'Since it has already been assigned a fluxBC'
            else
                write(*,*) 'Since it has already been assigned an invalid BC-type-Nr.: ',fesim%internals%thermal_dof_map(ll,2)
            end if
            err_code = 2
            goto 9999
        end if
        fesim%internals%thermal_dof_map(ll,2) = 3.                              ! Set the type to 1 to signal temperature bc
        fesim%internals%thermal_dof_map(ll,3) = fesim%thermal%convections(ii)%h ! Set the temperature value at the node
        fesim%internals%thermal_dof_map(ll,4) = fesim%thermal%convections(ii)%T_amb ! Set the temperature value at the node
    end do
 end do
 ! 124. Loop over MPCs
 do ii=1,size(fesim%randbedingungen%mpcs)
    ll = findloc(fesim%internals%thermal_dof_map(:,1),dble(fesim%randbedingungen%mpcs(ii)%dependend%nid),dim=1)   ! Find index of nid in thermal_dof_map
    if (ll == 0) then
        write(*,*) 'NID ',fesim%randbedingungen%mpcs(ii)%dependend%nid,' from MPCID: ',fesim%randbedingungen%mpcs(ii)%sid,' could not be found in thermal_dof_map!'
        err_code = 2
        goto 9999
    end if

    if (fesim%internals%thermal_dof_map(ll,2) .NE. 0.) then
        write(*,*) 'Cannot assign MPCID: ',fesim%randbedingungen%mpcs(ii)%sid,' to dependent Node with ID ',fesim%randbedingungen%mpcs(ii)%dependend%nid
        if (fesim%internals%thermal_dof_map(ll,2) .EQ. 1.) then
            write(*,*) 'Since it has already been assigned a temperatureBC'
        elseif (fesim%internals%thermal_dof_map(ll,2) .EQ. 2.) then
            write(*,*) 'Since it has already been assigned a fluxBC'
        elseif (fesim%internals%thermal_dof_map(ll,2) .EQ. 3.) then
            write(*,*) 'Since it has already been assigned a convectionBC'
        else
            write(*,*) 'Since it has already been assigned an invalid BC-type-Nr.: ',fesim%internals%thermal_dof_map(ll,2)
        end if
        err_code = 2
        goto 9999
    end if
    fesim%internals%thermal_dof_map(ll,2) = 4.
    fesim%internals%thermal_dof_map(ll,3) = fesim%randbedingungen%mpcs(ii)%sid
 end do
 ! 13. Loop over thermal_dof_map and create the permutation matrix PermutMat for the reordering of all matrices

! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol3_set_thermal_dof_map'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return
end subroutine sol3_set_thermal_dof_map
