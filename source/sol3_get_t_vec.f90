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
subroutine sol3_get_t_vec(fesim,Tf,Tc,is_free,is_bc)
use fesimulation_typen
#include "petsc/finclude/petscvec.h"
use petscvec
#include "petsc/finclude/petscis.h"
use petscis

implicit none

! Input + Output
type(fe_simulation), intent(inout)  :: fesim
IS,intent(inout)                    :: is_free,is_bc
! Output
Vec, intent(out)                :: Tf, Tc
! Internal
PetscErrorCode                  :: ierr
PetscMPIInt                     :: rank
PetscInt                        :: ii, jj, node_id, tc_index, is_bc_size, is_free_size, found
PetscInt, pointer               :: bc_indices(:)
double precision                :: tVal

! =================================================================================================
!
! Initialization
!
call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)

! Create Tf-Vector
call VecCreate(PETSC_COMM_SELF, Tf, ierr); CHKERRQ(ierr)
call ISGetSize(is_free,is_free_size,ierr); CHKERRQ(ierr)
call VecSetSizes(Tf, PETSC_DECIDE, is_free_size, ierr); CHKERRQ(ierr)
call VecSetFromOptions(Tf, ierr); CHKERRQ(ierr)
call VecSet(Tf, 0.0d0, ierr); CHKERRQ(ierr)

! Create Tc-Vector
call VecCreate(PETSC_COMM_SELF, Tc, ierr); CHKERRQ(ierr)
call ISGetSize(is_bc, is_bc_size, ierr); CHKERRQ(ierr)
!call ISView(is_bc, PETSC_VIEWER_STDOUT_SELF, ierr); CHKERRQ(ierr)
call VecSetSizes(Tc, PETSC_DECIDE, is_bc_size, ierr); CHKERRQ(ierr)
call VecSetFromOptions(Tc, ierr); CHKERRQ(ierr)
call VecSet(Tc, 0.0d0, ierr); CHKERRQ(ierr)

! Set the Tc-Vector
call ISGetIndicesF90(is_bc, bc_indices, ierr); CHKERRQ(ierr)
do ii=1,size(fesim%thermal%tbcs)
    tVal = fesim%thermal%tbcs(ii)%temperature
    node_id = fesim%thermal%tbcs(ii)%nid
    found = 0
    ! Find the position of node_id-1 in is_bc
    do jj = 0, is_bc_size-1
        if (bc_indices(jj+1) == node_id-1) then
            call VecSetValue(Tc, jj, tVal, INSERT_VALUES, ierr); CHKERRQ(ierr)
            found = 1
            exit
        end if
    end do

    if (found == 0) then
        write(*,*) 'ERROR: Temperature BC node not found in is_bc: ', node_id
        call ISRestoreIndices(is_bc, bc_indices, ierr); CHKERRQ(ierr)
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, &
                'Temperature BC node not found in is_bc')
    end if

end do
call ISRestoreIndicesF90(is_bc, bc_indices, ierr); CHKERRQ(ierr)

! Assemble Tf-vector
call VecAssemblyBegin(Tf, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(Tf, ierr); CHKERRQ(ierr)

! Assemble Tc-vector
call VecAssemblyBegin(Tc, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(Tc, ierr); CHKERRQ(ierr)

! View the vectors
!call VecView(Tf, PETSC_VIEWER_STDOUT_SELF, ierr); CHKERRQ(ierr)
!call VecView(Tc, PETSC_VIEWER_STDOUT_SELF, ierr); CHKERRQ(ierr)

return
end subroutine
