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
!> \file sol3_get_dof_is.f90
! =================================================================================================
!
!>  \brief      get the index sets which include free and constrained dof indices
!
!>  \details    get the index sets which include free and constrained dof indices
!
!>  \author     Fabian Reisener, WM, 2026
!
! =================================================================================================
subroutine sol3_get_dof_is(fesim,is_free,is_bc)
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscksp.h"
  use petscksp
  use fesimulation_typen
  use internals_typen
!
! =================================================================================================
!
! Data types
!
  implicit none
! Input
  type(fe_simulation),intent(in)  :: fesim
! Output
  IS, intent(out)                 :: is_free, is_bc
! Internal PETSC
  PetscInt, allocatable           :: freeList(:), bcList(:)
  PetscInt                        :: nfree, nbc
  PetscErrorCode                  :: ierr
! Internal other
  integer                         :: ii,jj,ll,nn
  integer                         :: err_code=0
! =================================================================================================
!
! Calculation
! Get free and constrained index sets from thermal_dof_map
  nfree = 0
  nbc   = 0
  do ii = 1, size(fesim%internals%thermal_dof_map, dim=1)
     if (fesim%internals%thermal_dof_map(ii,1) == 0) then
        nfree = nfree + 1
     else
        nbc   = nbc + 1
     end if
  end do

  allocate(freeList(nfree), bcList(nbc))

  nfree = 0
  nbc   = 0
  do ii = 1, size(fesim%internals%thermal_dof_map, dim=1)
     if (fesim%internals%thermal_dof_map(ii,1) == 0) then           ! If type is set to 0, that dof is free
        nfree = nfree + 1
        freeList(nfree) = fesim%internals%thermal_dof_map(ii,2) - 1 ! Shift to 0-based PETSc indexing from Node 1-based indices
     else
        nbc = nbc + 1
        bcList(nbc)     = fesim%internals%thermal_dof_map(ii,2) - 1 ! Shift to 0-based PETSc indexing from Node 1-based indices
     end if
  end do

  call ISCreateGeneral(PETSC_COMM_WORLD, nfree, freeList, PETSC_COPY_VALUES, is_free, ierr); CHKERRQ(ierr)
  call ISCreateGeneral(PETSC_COMM_WORLD, nbc,   bcList,   PETSC_COPY_VALUES, is_bc,   ierr); CHKERRQ(ierr)

! Cleanup
  deallocate(freeList, bcList)
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol3_get_dof_is'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
end if

return
end subroutine
