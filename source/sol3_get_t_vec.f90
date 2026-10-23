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
!> subroutine for creating the global temperature vector q
!
!> @details
!> subroutine for creating the global temperature vector
!
!> @author Fabian Reisener, TU Dresden, WM 2026
!
! =================================================================================================
subroutine sol3_get_t_vec(fesim,T)
use fesimulation_typen
#include "petsc/finclude/petscvec.h"
use petscvec

implicit none

! Input
type(fe_simulation), intent(in) :: fesim
! Output
Vec, intent(out)                :: T
! Internal
PetscErrorCode                  :: ierr
PetscMPIInt                     :: rank
integer                         :: ii, jj, nsize, node_id
double precision                :: tVal

! =================================================================================================
!
! Initialization
!
call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)

! Create T-Vector
nsize = fesim%num_nodes
call VecCreate(PETSC_COMM_WORLD, T, ierr); CHKERRQ(ierr)
call VecSetSizes(T, PETSC_DECIDE, nsize, ierr); CHKERRQ(ierr)
call VecSetFromOptions(T, ierr); CHKERRQ(ierr)

do ii=1,size(fesim%thermal%tbcs)
    tVal = fesim%thermal%tbcs(ii)%temperature
    node_id = fesim%thermal%tbcs(ii)%nid
    call VecSetValue(T,node_id-1,tVal,INSERT_VALUES,ierr); CHKERRQ(ierr)
end do

! Assemble vector
call VecAssemblyBegin(T, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(T, ierr); CHKERRQ(ierr)

return
end subroutine
