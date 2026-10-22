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
!> \file sol3_get_PermutMat.f90
! =================================================================================================
!
!>  \brief      subroutine for the permutation of thermal matrices
!
!>  \details    subroutine for the permutation of thermal matrices
!
!>  \author     Fabian Reisener, WM, 2026
!
! =================================================================================================

subroutine sol3_get_PermutMat(fesim,C,Cff,Cfc,Ccf,Ccc)
!
! =================================================================================================
!
! Include
!
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
type(fe_simulation),intent(inout)   :: fesim
Mat, intent(in)                     :: C
! Output
Mat, intent(out)                    :: Cff, Cfc, Ccf, Ccc
! Internal PETSC
PetscErrorCode                      :: ierr
!
! =================================================================================================
!
! Initialisation
!
!
! =================================================================================================
!
! Calculation

! Get index sets for free and constrained dofs
  call sol3_get_dof_is(fesim,is_free,is_bc)
! Get submatrices
  call MatCreateSubMatrix(C, is_free, is_free, MAT_INITIAL_MATRIX, Cff, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(C, is_free, is_bc,   MAT_INITIAL_MATRIX, Cfc, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(C, is_bc,   is_free, MAT_INITIAL_MATRIX, Ccf, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(C, is_bc,   is_bc,   MAT_INITIAL_MATRIX, Ccc, ierr); CHKERRQ(ierr)

! Cleanup
  call ISDestroy(is_free, ierr); CHKERRQ(ierr)
  call ISDestroy(is_bc,   ierr); CHKERRQ(ierr)

! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol3_get_permut_mat'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
end if

return
end subroutine
