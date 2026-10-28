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
!> \file sol3_output_vtk.f90
! =================================================================================================
!
!>  \brief      subroutine for writing quantities of thermal solution to vtk
!
!>  \details    subroutine for writing quantities of thermal solution to vtk
!
!>  \author     Fabian Reisener, WM, 2026
!
! =================================================================================================
subroutine sol3_output_vtk(fesim,Tf,Tc,is_free,is_bc)
!
! =================================================================================================
!
! Include
!
use fesimulation_typen
use vtk_variablen
#include "petsc/finclude/petscksp.h"
  use petscksp
#include "petsc/finclude/petscvec.h"
use petscvec
#include "petsc/finclude/petscis.h"
use petscis
!
! =================================================================================================
!
! Data types
!
implicit none
! Input + Output
type(fe_simulation),intent(inout)   :: fesim
! Input
Vec,intent(in)                      :: Tc,Tf
IS,intent(in)                       :: is_free,is_bc
! Output
! Internal
Vec                                 :: T
PetscErrorCode                      :: ierr
integer                             :: ii
PetscViewer                         :: viewer
PetscInt                            :: size_bc,size_free
IS                                  :: is_stride_bc, is_stride_free
VecScatter                          :: scatter_bc, scatter_free
logical                             :: writee = .TRUE.
!
! =================================================================================================
!
! Calculation
! Create the global vector T
call VecCreate(PETSC_COMM_SELF, T, ierr); CHKERRQ(ierr)
call VecSetSizes(T, PETSC_DECIDE, size(fesim%knoten%nodes), ierr); CHKERRQ(ierr)
call VecSetFromOptions(T, ierr); CHKERRQ(ierr)

! Get sizes of ISs
call ISGetSize(is_bc,size_bc,ierr); CHKERRQ(ierr)
call ISGetSize(is_free,size_free,ierr); CHKERRQ(ierr)
! Create stride ISs
call ISCreateStride(PETSC_COMM_SELF,size_bc  ,0,1,is_stride_bc  ,ierr); CHKERRQ(ierr)
call ISCreateStride(PETSC_COMM_SELF,size_free,0,1,is_stride_free,ierr); CHKERRQ(ierr)

! Create scatter_bc
call VecScatterCreate(Tc,is_stride_bc,T,is_bc,scatter_bc,ierr); CHKERRQ(ierr)
call VecScatterBegin(scatter_bc,Tc,T,INSERT_VALUES,SCATTER_FORWARD,ierr); CHKERRQ(ierr)
call   VecScatterEnd(scatter_bc,Tc,T,INSERT_VALUES,SCATTER_FORWARD,ierr); CHKERRQ(ierr)
call VecScatterDestroy(scatter_bc,ierr); CHKERRQ(ierr)
! Create scatter_free
call VecScatterCreate(Tf,is_stride_free,T,is_free,scatter_free,ierr); CHKERRQ(ierr)
call VecScatterBegin(scatter_free,Tf,T,INSERT_VALUES,SCATTER_FORWARD,ierr); CHKERRQ(ierr)
call   VecScatterEnd(scatter_free,Tf,T,INSERT_VALUES,SCATTER_FORWARD,ierr); CHKERRQ(ierr)
call VecScatterDestroy(scatter_free,ierr); CHKERRQ(ierr)

! Write
if (writee) then
    call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'Tc2.txt', viewer, ierr); CHKERRQ(ierr)
    call VecView(Tc, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'Tf2.txt', viewer, ierr); CHKERRQ(ierr)
    call VecView(Tf, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'is_free.txt', viewer, ierr); CHKERRQ(ierr)
    call ISView(is_free, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'is_bc.txt', viewer, ierr); CHKERRQ(ierr)
    call ISView(is_bc, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'T.txt', viewer, ierr); CHKERRQ(ierr)
    call VecView(T, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
end if

! Cleanup
call VecDestroy(T,ierr); CHKERRQ(ierr)

return
end subroutine
