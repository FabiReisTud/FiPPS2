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
!> Adds temperature loadvector to the total loadvector
!
!> @details
!> Temperature load vector of each element is calculated by the corresponding
!> subroutine and then added to the total loadvector
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: temperature_loads.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine temperature_loads(fesim, lcID, pos, lvec)
! =================================================================================================
!   Input:      hh:   Index of current load
!               pos:  Position of loadvector in loadvector array
!
!   In-/Output: lvec: Total loadvector for current load
!
! =================================================================================================
! use
!
use konstanten
use fesimulation_typen
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
! Input
!
 type(fe_simulation)                             :: fesim
 integer, intent(in)                             :: lcID !< Loadcase-ID
 integer, intent(in)                             :: pos !< Position of loadvector in loadvector array
!
! Input + Output
 double precision, intent(inout)                 :: lvec(1:fesim%num_dof,1:fesim%num_loadcases)
!
! Internal
!
 double precision                                :: nodal_temperatures(fesim%num_nodes)
 double precision                                :: beam2_temperatures(size(fesim%elemente%beam2s,1))
 double precision                                :: quad8_temperatures(size(fesim%elemente%quad8s,1))
 double precision                                :: lsolid20_temperatures(size(fesim%elemente%lsolid20s,1))
 double precision, dimension(:), allocatable     :: elementtemperatures_vtk
 integer                                         :: eid_offset
 integer                                         :: ii
!
! =================================================================================================
!
! Initialisation
!
! =================================================================================================
!
! Calculation
!
  nodal_temperatures(:)    = 0.d0
  beam2_temperatures(:)    = 0.d0
  quad8_temperatures(:)    = 0.d0
  lsolid20_temperatures(:) = 0.d0

  ! Get nodal temperatures from the thermal simulation
  if (fesim%thermal%enabled) then
    nodal_temperatures = fesim%thermal%temperature
  end if

  ! Get nodal temperature at every node
  if (fesim%is_temperature .eqv. .true.) then
    call get_node_temperatures(fesim, lcID, nodal_temperatures)
    if (fesim%ausgabe%outputVTK .eqv. .true.) then
      ! Write nodal temperature to vtk
      call vtkoutpointscalar(nodal_temperatures, fesim%num_nodes, 'Temperatur')
    endif
  end if
!
  ! Get elemental temperatures at beam2 elements
  if (fesim%is_beam2temp .eqv. .true.) call get_beam2_temperatures(fesim, lcID, beam2_temperatures)
!
  ! Get elemental temperatures at quad8 elements
  if (fesim%is_quad8temp .eqv. .true.) call get_quad8_temperatures(fesim, lcID, quad8_temperatures)
!
  ! Get elemental temperatures at lsolid20 elements
  if (fesim%is_lsolid20temp .eqv. .true.) call get_lsolid20_temperatures(fesim, lcID, lsolid20_temperatures)
!
  ! Write elemental temperatures to vtk file
  if ((fesim%is_beam2temp .eqv. .true.) .or. (fesim%is_quad8temp .eqv. .true.) .or. (fesim%is_lsolid20temp .eqv. .true.)) then
    if (fesim%ausgabe%outputVTK .eqv. .true.) then
      allocate(elementtemperatures_vtk(fesim%num_elements))
      elementtemperatures_vtk(1:size(fesim%elemente%beam2s,1)) = beam2_temperatures(:)
      eid_offset = size(fesim%elemente%beam2s,1)
      elementtemperatures_vtk(1+eid_offset:eid_offset+size(fesim%elemente%quad8s,1)) = quad8_temperatures(:)
      eid_offset = eid_offset + size(fesim%elemente%quad8s,1)
      elementtemperatures_vtk(1+eid_offset:eid_offset+size(fesim%elemente%lsolid20s,1)) = lsolid20_temperatures(:)
      ! Write elemental temperature to vtk
      call vtkoutcelldata(elementtemperatures_vtk, fesim%num_elements, 'Temperatur')
      deallocate(elementtemperatures_vtk)
    end if
  end if

  call beam2_temploadvec_control(fesim, lvec, pos, nodal_temperatures, beam2_temperatures)

  call quad8_temploadvec_control(fesim, lvec, pos, nodal_temperatures, quad8_temperatures)

  call lsolid20_temploadvec_control(fesim, lvec, pos, nodal_temperatures, lsolid20_temperatures)

end subroutine temperature_loads
