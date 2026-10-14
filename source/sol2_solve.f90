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
!
!> @details
!
!> @author
!
!> $Id: sol2_solve.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol2_solve(eps, ksp, KaaS, KgaaS, numEigVal)
!
! Use
!
!
! =================================================================================================
!
! Include
!
#include "slepc/finclude/slepceps.h"
  use slepceps
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!

  Mat            :: KaaS
  Mat            :: KgaaS
  EPS            :: eps
  PetscErrorCode :: ierr
  PetscInt       :: nev
  ST             :: spt
  KSP            :: ksp
#ifndef mod_slepc
  PC             :: prec
#endif

  integer        :: numEigVal
  integer        :: err_code=0
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation
!
! Buckling solution
!

  !** Create eigensolver context
  call EPSCreate(PETSC_COMM_WORLD,eps,ierr); CHKERRQ(ierr)

  ! ** Set operators. In this case, it is a standard eigenvalue problem
  call EPSSetOperators(eps,KgaaS,KaaS,ierr); CHKERRQ(ierr)
  call EPSSetProblemType(eps,EPS_GHEP,ierr); CHKERRQ(ierr)

  !Durch das negative Reziproke wird dieser zum kleinsten positiven Eigenwert
  call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL,ierr); CHKERRQ(ierr)
  nev = numEigVal
  call EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE,ierr); CHKERRQ(ierr)

  ! ** Set solver parameters at runtime
  call EPSSetFromOptions(eps,ierr); CHKERRQ(ierr)

  call EPSGetST(eps,spt,ierr); CHKERRQ(ierr)

#ifdef mod_slepc
! Der nächste Abschnitt ist mit der modifzizierten Version von SLEPC zuverwenden
!_______________________________________________________________________________
  call STSetKSP(spt,ksp,ierr); CHKERRQ(ierr)
#else
! Der nächste Abschnitt ist mit der Standardversion von SLEPC zuverwenden
!________________________________________________________________________
  call KSPDestroy(ksp,ierr); CHKERRQ(ierr)
  call STGetKSP(spt,ksp,ierr); CHKERRQ(ierr)
  call KSPSetType(ksp,KSPPREONLY,ierr); CHKERRQ(ierr)
  call KSPGetPC(ksp,prec,ierr); CHKERRQ(ierr)
  call PCSetType(prec,PCCHOLESKY,ierr); CHKERRQ(ierr)
  call PCFactorSetMatSolverType(prec,MATSOLVERMUMPS,ierr); CHKERRQ(ierr)

  call KSPSetErrorIfNotConverged(ksp,PETSC_TRUE,ierr); CHKERRQ(ierr) ! Wenn der Workspace für MUMPS nicht ausreicht wirft PETSC 3.7.6 sonst keinen Fehler
  call PCSetErrorIfFailure(prec,PETSC_TRUE,ierr); CHKERRQ(ierr) ! Wenn der Workspace für MUMPS nicht ausreicht wirft PETSC 3.7.6 sonst keinen Fehler

  call EPSSetUp(eps,ierr); CHKERRQ(ierr)
#endif

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Solve the eigensystem
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call EPSSolve(eps,ierr); CHKERRQ(ierr)

!
! =================================================================================================
!
! Error handling
!
9999 continue

  if (err_code /= 0) then

    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'sol2_solve'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop

  end if

end subroutine sol2_solve
