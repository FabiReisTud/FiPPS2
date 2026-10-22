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
!> \file sol3_assemble.f90
! =================================================================================================
!
!>  \brief      subroutine for assembling all matrices and vectors of the thermal problem
!
!>  \details    subroutine for assembling all matrices and vectors of the thermal problem
!>              1. the temperatur vectors
!>              11. Tf = free temperature dofs
!>              12. Tc = constrained temperature dofs
!>              2. the load vectors (if convection, based on prev. temp)
!>              21. Qf = free flux dofs (between elements)
!>              22. Qc = constrained dofs (at boundaries, adiabitc=0 if not defined otherwise)
!>              3. the conduction matrices
!>              31. Cff = coupling between free and free temperature dofs
!>              31. Cfc = coupling between constrained and free temperature dofs
!
!>  \author     Fabian Reisener, WM, 2026
!
! =================================================================================================

subroutine sol3_assemble(fesim,Tf,Tc,Qf,Qc,Cff,Cfc)
!
! =================================================================================================
!
! Include
!
  use petscmat
  use petscvec
  use petsc
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
  Mat, intent(out)                  :: Cff, Cfc
  Vec, intent(out)                  :: Tf,Tc,Qf,Qc
! Internal Vectors/Matrices
  Mat                               :: C
  Vec                               :: T,Q
  Mat                               :: PermutMat
! Internal other
  integer                           :: ii,jj,ll,nn
  integer                           :: err_code
  integer,dimension(:)              :: node_ids
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
 end do

 ! 12. Loop over all thermal boundary conditions, assign values to the row in thermal_dof_map if assigned
 ! 121. Loop over temperature BCs
 do ii=1,size(fesim%thermal%tbcs)
    ll = findloc(fesim%internals%thermal_dof_map(:,1),fesim%thermal%tbcs(ii)%nid,dim=1) ! Find index of nid in thermal_dof_map
    fesim%internals%thermal_dof_map(ll,2) = 1.                                          ! Set the type to 1 to signal temperature bc
    fesim%internals%thermal_dof_map(ll,3) = fesim%thermal%tbcs(ii)%temperature          ! Set the temperature value at the node
 end do

 ! 122. Loop over flux BCs
 do ii=1,size(fesim%thermal%convections)
    node_ids(1) = fesim%thermal%fluxes(ii)%nid1
    node_ids(2) = fesim%thermal%fluxes(ii)%nid2
    node_ids(3) = fesim%thermal%fluxes(ii)%nid3
    node_ids(4) = fesim%thermal%fluxes(ii)%nid4
    node_ids(5) = fesim%thermal%fluxes(ii)%nid5
    node_ids(6) = fesim%thermal%fluxes(ii)%nid6
    node_ids(7) = fesim%thermal%fluxes(ii)%nid7
    node_ids(8) = fesim%thermal%fluxes(ii)%nid8

    do jj=1,8
        ll = findloc(fesim%internals%thermal_dof_map(:,1),node_ids(jj),dim=1)   ! Find index of nid in thermal_dof_map
        if (fesim%internals%thermal_dof_map(ll,2) .NE. 0.) then
            write(*,*) ''
        end if
        fesim%internals%thermal_dof_map(ll,2) = 2.                              ! Set the type to 1 to signal temperature bc
        fesim%internals%thermal_dof_map(ll,3) = fesim%thermal%fluxes(ii)%q      ! Set the temperature value at the node
    end do

 end do

 ! 123. Loop over convection BCs

 ! 124. Loop over MPCs

 ! 13. Loop over thermal_dof_map and create the permutation matrix PermutMat for the reordering of all matrices

! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol3_assemble'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return
end subroutine sol3_assemble
