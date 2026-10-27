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
!> $Id: sol1_solve.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_solve(fesim, ksp, K, u, f, rank, scloop)
!
! Use
!
use konstanten
use fesimulation_typen
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscksp.h"
  use petscksp
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
  type(fe_simulation)                   :: fesim
  PetscMPIInt, intent(in)               :: rank
  integer, intent(in)                   :: scloop
  Mat                                   :: K
  Mat                                   :: mat
  Vec                                   :: u
  Vec                                   :: f

  KSP                                   :: ksp
  PetscErrorCode                        :: ierr
  PC                                    :: prec

  integer                               :: err_code=0
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

  if (fesim%lasten%subcases(scloop)%upmats .eq. .true.) then
    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp,K,K,ierr); CHKERRQ(ierr)

    call KSPSetType(ksp,KSPPREONLY,ierr); CHKERRQ(ierr)
    call KSPGetPC(ksp,prec,ierr); CHKERRQ(ierr)
    call PCSetType(prec,PCCHOLESKY,ierr); CHKERRQ(ierr)
    call PCFactorSetMatSolverType(prec,MATSOLVERMUMPS,ierr); CHKERRQ(ierr)
    call PCFactorSetUpMatSolverType(prec,ierr); CHKERRQ(ierr) ! notwendig, um die Methode MatMumpsSetIcntl ohne Fehler aufrufen zu können

    call KSPSetErrorIfNotConverged(ksp,PETSC_TRUE,ierr); CHKERRQ(ierr) ! Wenn der Workspace für MUMPS nicht ausreicht wirft PETSC 3.7.6 sonst keinen Fehler
    call PCSetErrorIfFailure(prec,PETSC_TRUE,ierr); CHKERRQ(ierr) ! Wenn der Workspace für MUMPS nicht ausreicht wirft PETSC 3.7.6 sonst keinen Fehler
!    call PCFactorGetMatrix(prec,mat,ierr); CHKERRQ(ierr) ! notwendig für die nachfolgende Methode
!    call MatMumpsSetIcntl(mat,14,40,ierr); CHKERRQ(ierr) ! Erhöhung des Workspace für MUMPS, da es sonst Fehler gibt, kann auch über "-mat_mumps_icntl_14 50" direkt beim FiPPS Aufruf gesetzt werden

    call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
  end if

 ! Solution
 ! KSPSolve(ksp, right-hand-side-vector, solution vector)
  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'Loesen des Gleichungssystems'
  call KSPSolve(ksp,f,u,ierr); CHKERRQ(ierr)
  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'Lastfall ', scloop, ' geloest.'

!  call KSPDestroy(ksp,ierr); CHKERRQ(ierr)

!
! =================================================================================================
!
! Error handling
!
9999 continue

  if (err_code /= 0) then

    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'sol1_solve'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop

  end if

end subroutine sol1_solve
