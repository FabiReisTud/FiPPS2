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
!> control subroutine for computing element conduction matrix for quadratic
!> 20-node solid element and assemble it to the total
!> stiffness matrix
!
!> @details
!> subroutine calculates the element conduction matrix for quadratic 20-node
!> solid element and assembles it to the total stiffness matrix
!
!> @author Fabian Reisener, TU Dresden, WM 2026
!
! =================================================================================================
subroutine tsolid20_cond_control(fesim,CaaS)
! =================================================================================================
! Use
!
use konstanten
use fesimulation_typen
use globale_variablen
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscmat.h"
use petscmat
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation), intent(in)         :: fesim
!
! Output
!

!
! Input + Output
!
Mat                 :: CaaS !< CaaS
!
! Internal
!
integer, dimension(20)              :: node_ids
double precision, dimension(20,3)   :: node_coords
double precision, dimension(20,20)  :: tsolid20_C


integer                             :: jj, kk,lsolid20Num,elemNum
integer                             :: err_code=0

PetscMPIInt                         :: rank
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
!
 call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

 do jj = 1,size(fesim%elemente%lsolid20s,1)
    if (mpi_lsolid20_proc_dist(jj) .NE. rank) cycle

    call tsolid20_cond(fesim,node_coords,lsolid20Num,elemNum,tsolid20_C)
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'tsolid20_cond_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return

end subroutine
