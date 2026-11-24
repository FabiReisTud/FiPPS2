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
subroutine sol3_output(fesim,Tf,Tc,is_free,is_bc)
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
PetscInt                            :: size_bc,size_free,size_all
IS                                  :: is_stride_bc, is_stride_free
VecScatter                          :: scatter_bc, scatter_free
PetscScalar, pointer                :: vec_array(:)
double precision,allocatable        :: T_new(:)
character(100)                      :: fname
double precision                    :: eps
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

!!!!! Write T-field to the thermal state

! Convert T to an array
call VecGetSize(T,size_all,ierr); CHKERRQ(ierr)
allocate(T_new(size_all),vec_array(size_all))
call VecGetArrayF90(T, vec_array, ierr); CHKERRQ(ierr)
T_new = vec_array
call VecRestoreArrayF90(T, vec_array, ierr); CHKERRQ(ierr) ! Restore
! On the first iteration, allocate the temperature fields
if (fesim%thermal%iteration==1) then
    allocate(fesim%thermal%temperature(fesim%num_nodes),fesim%thermal%temperature_old(fesim%num_nodes))
    fesim%thermal%temperature_old = fesim%thermal%tref ! On the first iteration, set the old temperature to the reference temperature
else
    fesim%thermal%temperature_old = fesim%thermal%temperature
end if

!!!!! Check for Convergence
if (fesim%thermal%iterative .AND. fesim%thermal%iteration>=2) then
    eps = dot_product(T_new - fesim%thermal%temperature_old,T_new - fesim%thermal%temperature_old)
    eps = eps / dot_product(fesim%thermal%temperature_old,fesim%thermal%temperature_old)
    write(*,*) 'Thermal temperature change norm = ',eps
    write(*,*) '                      Tolerance = ',fesim%thermal%tolerance
    if (eps < fesim%thermal%tolerance) then
        fesim%thermal%successful = .TRUE.
    end if
elseif (fesim%thermal%iterative .EQ. .FALSE.) then
    fesim%thermal%successful = .TRUE.
end if

!!!!! Determine the new vectors, apply relaxation
if (fesim%thermal%successful .EQ. .FALSE.) then ! If not converged, apply relaxation
    fesim%thermal%temperature = fesim%thermal%relaxation * fesim%thermal%temperature_old + (1 - fesim%thermal%relaxation) * T_new
else ! If converged, dont apply relaxation
    fesim%thermal%temperature = T_new
end if

!!!!! Output to VTK
if (writee) then
    write(fname, '("output_thermal_", I0, ".vtk")') fesim%thermal%iteration
    call vtkoutprep(fname)               ! VTK-Datei öffnen und vorbereiten
    call vtkoutmesh(fesim)                              ! Write mesh
    ! Write T_new to the VTK
    call vtkoutpointscalar(T_new, fesim%num_nodes, 'T_new')
    call vtkoutpointscalar(fesim%thermal%temperature, fesim%num_nodes, 'temperature')
    call vtkoutpointscalar(fesim%thermal%temperature_old, fesim%num_nodes, 'temperature_old')
end if



return
end subroutine
