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
Vec                                 :: Ttemp(2)
IS                                  :: ISTemp(2)
PetscErrorCode                      :: ierr
integer                             :: ii
!
! =================================================================================================
!
! Calculation
call VecCreate(PETSC_COMM_WORLD,T,ierr); CHKERRQ(ierr) ! Create T-Vector
call VecSetSizes(T, PETSC_DECIDE, size(fesim%knoten%nodes), ierr); CHKERRQ(ierr)
call VecSetFromOptions(T, ierr); CHKERRQ(ierr)
!call VecSet(T, 0.0d0, ierr); CHKERRQ(ierr)
call VecISCopy(T, is_bc, SCATTER_FORWARD, Tc, ierr); CHKERRQ(ierr)
call VecISCopy(T, is_free, SCATTER_FORWARD, Tf, ierr); CHKERRQ(ierr)

return
end subroutine
