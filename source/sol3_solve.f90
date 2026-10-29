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
!> \file sol3_solve.f90
! =================================================================================================
!
!>  \brief      subroutine for solving the thermal problem
!
!>  \details    subroutine for solving the thermal problem
!
!>  \author     Fabian Reisener, WM, 2026
!
! =================================================================================================
subroutine sol3_solve(fesim,Tf,Tc,Qf,Qc,Cff,Cfc)
!
! =================================================================================================
!
! Include
!
use fesimulation_typen
#include "petsc/finclude/petscksp.h"
  use petscksp
#include "petsc/finclude/petscvec.h"
  use petscvec
#include "petsc/finclude/petscmat.h"
  use petscmat
!
! =================================================================================================
!
! Data types
!
implicit none
! Input + Output
type(fe_simulation),intent(inout)   :: fesim
! Input
Vec,intent(in)                      :: Qf,Qc,Tc
Mat,intent(in)                      :: Cff,Cfc
! Output
Vec,intent(inout)                   :: Tf
! Internal
Vec                                 :: RHS,RHS2 ! Right-hand-side
KSP                                 :: ksp      ! Krylov space
PC                                  :: prec     ! precondition context
PetscErrorCode                      :: ierr
PetscReal                           :: alpha=-1.0, tol=1e-5
PetscViewer                         :: viewer
logical                             :: write=.TRUE.
Vec                                 :: residual
PetscScalar                         :: rnorm,rhsnorm
!
! =================================================================================================
!
! Calculation

! Start the RHS
call VecDuplicate(Qf, RHS, ierr); CHKERRQ(ierr)
call VecCopy(Qf, RHS, ierr);     CHKERRQ(ierr)
RHS = Qf

call MatCreateVecs(Cfc, PETSC_NULL_VEC, RHS2, ierr); CHKERRQ(ierr)
call MatMult(Cfc,Tc,RHS2, ierr); CHKERRQ(ierr)   ! Multiply Cfc with Tc, save to RHS2

! Finish RHS
call VecAXPY(RHS,alpha,RHS2, ierr); CHKERRQ(ierr) ! Subtract RHS2 from RHS(==Qf)

! Create linear solver
call KSPCreate(PETSC_COMM_WORLD,ksp, ierr); CHKERRQ(ierr)
! Set defining Matrix (==A))
call KSPSetOperators(ksp,Cff,Cff, ierr); CHKERRQ(ierr)
! Solver settings (copied from sol1_solve)
call KSPSetType(ksp,KSPPREONLY,ierr); CHKERRQ(ierr)
call KSPGetPC(ksp,prec,ierr); CHKERRQ(ierr)
call PCSetType(prec,PCCHOLESKY,ierr); CHKERRQ(ierr)
call PCFactorSetMatSolverType(prec,MATSOLVERMUMPS,ierr); CHKERRQ(ierr)
call PCFactorSetUpMatSolverType(prec,ierr); CHKERRQ(ierr) ! notwendig, um die Methode MatMumpsSetIcntl ohne Fehler aufrufen zu können

!call KSPSetErrorIfNotConverged(ksp,PETSC_TRUE,ierr); CHKERRQ(ierr) ! Wenn der Workspace für MUMPS nicht ausreicht wirft PETSC 3.7.6 sonst keinen Fehler
!call PCSetErrorIfFailure(prec,PETSC_TRUE,ierr); CHKERRQ(ierr) ! Wenn der Workspace für MUMPS nicht ausreicht wirft PETSC 3.7.6 sonst keinen Fehler

call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
! Solve
call MatCreateVecs(Cff, PETSC_NULL_VEC, Tf, ierr); CHKERRQ(ierr)
call KSPSolve(ksp,RHS,Tf, ierr); CHKERRQ(ierr)

! calculate residual
call MatCreateVecs(Cff, residual, PETSC_NULL_VEC, ierr); CHKERRQ(ierr)
! residual = Cff * Tf - RHS
call MatMult(Cff, Tf, residual, ierr); CHKERRQ(ierr)
call VecAXPY(residual, -1.0d0, RHS, ierr); CHKERRQ(ierr)
! Compute residual norm
call VecNorm(residual, NORM_2, rnorm, ierr); CHKERRQ(ierr)
! Compute RHS norm
call VecNorm(RHS, NORM_2, rhsnorm, ierr); CHKERRQ(ierr)

if (rhsnorm .gt. 0.0d0) then
    write(*,*) 'Relative residual = ', rnorm/rhsnorm
else
    write(*,*) 'Absolute residual = ', rnorm
endif

call VecDestroy(residual, ierr); CHKERRQ(ierr)

if (write) then
    ! Export Tf
    call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'Tf.txt', viewer, ierr)
    call VecView(Tf, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)
endif

! Final Cleanup
call VecDestroy(RHS,ierr); CHKERRQ(ierr)
call VecDestroy(RHS2,ierr); CHKERRQ(ierr)
call KSPDestroy(ksp,ierr); CHKERRQ(ierr)

return
end subroutine
