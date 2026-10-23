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
#include "petsc/finclude/petscksp.h"
  use petscksp
#include "petsc/finclude/petscis.h"
  use petscis
  use fesimulation_typen
!
  implicit none
!
! =================================================================================================
!
! Data types
!
! Input
  type(fe_simulation),intent(inout) :: fesim
! Output
  Mat, intent(out)                  :: Cff, Cfc
  Vec, intent(out)                  :: Tf,Tc,Qf,Qc
! Internal Vectors/Matrices
  Mat                               :: C, Ccf, Ccc
  Vec                               :: T,Q
  Mat                               :: PermutMat
  IS                                :: is_free, is_bc
  PetscErrorCode                    :: ierr
! Internal other
  integer                           :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
!
! =================================================================================================
!
! Interfaces
interface
  subroutine sol3_set_thermal_dof_map(fesim)
    use fesimulation_typen
    implicit none
    type(fe_simulation), intent(inout) :: fesim
  end subroutine sol3_set_thermal_dof_map

  subroutine sol3_get_dof_is(fesim,is_free,is_bc)
    use fesimulation_typen
    use petscis
    implicit none
    type(fe_simulation), intent(in) :: fesim
    IS, intent(out)                 :: is_free, is_bc
  end subroutine sol3_get_dof_is

  subroutine sol3_tsolid20_cond_control(fesim,C)
    use fesimulation_typen
    use petscmat
    implicit none
    type(fe_simulation), intent(in) :: fesim
    Mat,intent(out)                 :: C
  end subroutine sol3_tsolid20_cond_control
end interface
! =================================================================================================
!
! Calculation

! Set the thermal_dof_map
  call sol3_set_thermal_dof_map(fesim)
! Get index sets for free and constrained dofs
  call sol3_get_dof_is(fesim,is_free,is_bc)
! Assemble C-Matrix
  call sol3_tsolid20_cond_control(fesim,C)
! Get C-submatrices
  call MatCreateSubMatrix(C, is_free, is_free, MAT_INITIAL_MATRIX, Cff, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(C, is_free, is_bc,   MAT_INITIAL_MATRIX, Cfc, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(C, is_bc,   is_free, MAT_INITIAL_MATRIX, Ccf, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(C, is_bc,   is_bc,   MAT_INITIAL_MATRIX, Ccc, ierr); CHKERRQ(ierr)

! Cleanup
  call ISDestroy(is_free, ierr); CHKERRQ(ierr)
  call ISDestroy(is_bc,   ierr); CHKERRQ(ierr)
  call MatDestroy(C,ierr); CHKERRQ(ierr)
  call MatDestroy(Ccf,ierr); CHKERRQ(ierr)
  call MatDestroy(Ccc,ierr); CHKERRQ(ierr)
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
