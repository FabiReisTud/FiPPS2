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
!> Subroutine computes or assignes the temperatures regarding to the specified
!> load case to each node.
!
!> @details
!> Subroutine computes or assignes the temperatures regarding to the specified
!> load case to each node.
!
!> @author Thomas Dylla, TU Dresden, Diplomarbeit, 31.01.2014
!
!> @author Florian Dexl, 2015
!
!> @author Andreas Hauffe, (Quad8), 16.05.2017
!
!> $Id: get_node_temperatures.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine get_node_temperatures(fesim,lcID,nodal_temperatures)
! =================================================================================================
! use
!
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
type(fe_simulation)                    :: fesim
integer, intent(in)                    :: lcID !< Load case-ID
!
! Output
double precision, dimension(fesim%num_nodes), intent(out) :: nodal_temperatures !< Temperatures on each node for this load case
!
! Internal
!
integer                                :: ii, jj, j1, j2, nid
integer, dimension(:),allocatable      :: node_ids
logical, dimension(:),allocatable      :: temp_flag
double precision                       :: T_ref

integer                                :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
! nodal_temperatures = 0.D0 ! FR: this is redundant, since nodal_temperature is already initialized in "temperature_loads"
 T_ref = 0.d0
 allocate(temp_flag(fesim%num_nodes))
 temp_flag = .FALSE.
!
! =================================================================================================
!
! Calculation
!
 ! Find reference temperature for current load case
 ! If no reference temperature is given, it is set to 0
 do ii = 1,size(fesim%lasten%nodetemps,1)
   do jj = 1,size(fesim%lasten%loads,1)
     if ((fesim%lasten%loads(jj)%lcid .eq. lcID) .and. (fesim%lasten%loads(jj)%lidi .eq. fesim%lasten%nodetemps(ii)%lid) .and. (fesim%lasten%nodetemps(ii)%nid .eq. 0)) then
       T_ref = fesim%lasten%nodetemps(ii)%temp ! Tref is set from the first entry(nid==0) of the nodetemps-array.
     end if
   end do
 end do

 ! Loop over each temperature load to get all node temperatures assigned to this load case;
 ! Set temperature at node to (given temperature - reference temperature);
 ! Temperature = 0 is set to all nodes where no explicit temperature was assigned
 do ii=1,size(fesim%lasten%nodetemps,1)
   do jj = 1,size(fesim%lasten%loads,1)
     if((fesim%lasten%loads(jj)%lcid .eq. lcID) .and. (fesim%lasten%loads(jj)%lidi .eq. fesim%lasten%nodetemps(ii)%lid) .and. (fesim%lasten%nodetemps(ii)%nid .ne. 0)) then
       ! Assign temperature to node
       nid = fesim%lasten%nodetemps(ii)%nid
       nodal_temperatures(nid) = fesim%lasten%nodetemps(ii)%temp - T_ref
       temp_flag(nid) = .TRUE.
     end if
   end do
 end do

! Compute temperatures on elements (only needed in case of elements with mid nodes)

! quad8 (mid nodes without explicit temperature assignment get the average temperature of the neighbour nodes)
 if(allocated(fesim%elemente%quad8s)) then
   allocate(node_ids(8))
   ! Loop over each element
   do ii=1,size(fesim%elemente%quad8s,1)
     ! Get node IDs
     node_ids = fesim%elemente%quad8s(ii)%nids
     ! Compute mid node temperatures as mean value of the neighbour nodes, if no
     ! node temperature was defined before
     if (.NOT. temp_flag(node_ids(5))) then
       nodal_temperatures(node_ids(5)) = (nodal_temperatures(node_ids(1)) + nodal_temperatures(node_ids(2)))/2.d0
     end if
     if (.NOT. temp_flag(node_ids(6))) then
       nodal_temperatures(node_ids(6)) = (nodal_temperatures(node_ids(2)) + nodal_temperatures(node_ids(3)))/2.d0
     end if
     if (.NOT. temp_flag(node_ids(7))) then
       nodal_temperatures(node_ids(7)) = (nodal_temperatures(node_ids(3)) + nodal_temperatures(node_ids(4)))/2.d0
     end if
     if (.NOT. temp_flag(node_ids(8))) then
       nodal_temperatures(node_ids(8)) = (nodal_temperatures(node_ids(4)) + nodal_temperatures(node_ids(1)))/2.d0
     end if
   end do
   deallocate(node_ids)
 end if

! Lsolid20 (mid nodes without explicit temperature assignment get the average temperature of the neighbour nodes)
 if(allocated(fesim%elemente%lsolid20s)) then
   allocate(node_ids(20))
   ! Loop over each element
   do ii=1,size(fesim%elemente%lsolid20s,1)
     ! Get node IDs
     node_ids = fesim%elemente%lsolid20s(ii)%nids
     ! Compute mid node temperatures as mean value of the neighbour nodes, if no
     ! node temperature was defined before
     do jj=9,20
       if (jj .LE. 16) then
         ! Midside nodes in plane
         j1 = jj-8
         j2 = jj-7
         if (MOD((j2-1),4) == 0) j2 = j2-4
       else
         ! Midside nodes in thickness direction
         j1 = jj-16
         j2 = jj-12
       end if
       if (.NOT. temp_flag(node_ids(jj))) then
         nodal_temperatures(node_ids(jj)) = (nodal_temperatures(node_ids(j1)) + nodal_temperatures(node_ids(j2)))/2.d0
       end if
     end do
   end do
   deallocate(node_ids)
 end if

 deallocate(temp_flag)

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'get_node_temperatures'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return

end subroutine get_node_temperatures
