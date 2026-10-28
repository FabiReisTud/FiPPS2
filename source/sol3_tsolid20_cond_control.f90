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
subroutine sol3_tsolid20_cond_control(fesim,C)
! =================================================================================================
! Use
!
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
! Interfaces
interface
  subroutine sol3_tsolid20_cond(fesim, node_coords, node_ids, lsolid20Num, elemNum, k, tsolid20_C)
    use integration_schemes
    use fesimulation_typen
    implicit none
    type(fe_simulation), intent(in)      :: fesim
    double precision, intent(in)         :: node_coords(20,3)
    integer, intent(in)                  :: node_ids(20)
    integer, intent(in)                  :: lsolid20Num, elemNum
    type(thermal_curve_type), intent(in) :: k
    double precision, intent(out)        :: tsolid20_C(20,20)
  end subroutine sol3_tsolid20_cond
end interface
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
Mat,intent(out)                 :: C !< C
!
! Input + Output
!
!
! Internal
!
integer, dimension(20)              :: node_ids
double precision, dimension(20,3)   :: node_coords
double precision, dimension(20,20)  :: tsolid20_C

integer                             :: ii,jj, kk,ll,e,lsolid20Num,elemNum
PetscInt                            :: nnodes_per_elem,nnodes_global
PetscInt, allocatable               :: d_nnz(:), o_nnz(:)
integer                             :: err_code=0
integer                             :: pid, tmidid
type(thermal_curve_type)            :: k

PetscMPIInt                         :: rank
PetscErrorCode                      :: ierr
!
! =================================================================================================
!
! Initialisation
nnodes_per_elem=20
nnodes_global = size(fesim%knoten%nodes)
allocate(d_nnz(nnodes_global), o_nnz(nnodes_global))
d_nnz = 60
o_nnz = 60
!
!
! =================================================================================================
!
! Calculation
!
 call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

! Create global conductivity matrix C
call MatCreateAIJ(PETSC_COMM_WORLD, &
                  PETSC_DECIDE, PETSC_DECIDE, &
                  nnodes_global, nnodes_global, &
                  10, d_nnz, &
                  10, o_nnz, &
                  C, ierr); CHKERRQ(ierr)
call MatSetOption(C, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
!call MatSetOption(C, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr); CHKERRQ(ierr)
!call MatSetOption(C, MAT_ROW_ORIENTED, PETSC_FALSE, ierr); CHKERRQ(ierr)

! Iter through all the elements
 do jj = 1,size(fesim%elemente%lsolid20s,1)
    if (mpi_lsolid20_proc_dist(jj) .NE. rank) cycle
    lsolid20Num = fesim%elemente%lsolid20s(jj)%eid
    pid = fesim%elemente%lsolid20s(jj)%pid ! PropertyID des Elements
    tmidid = fesim%thermal%propthermals(pid)%tmidid ! ID des Thermischen Materials initialisieren
    k = fesim%thermal%materials_td(tmidid)%k! Wärmeleitkoeffizientenkurve holen

    if (tmidid .EQ. 0) then
       write(*,*) 'for lsolid20 number:                  ', lsolid20Num
       write(*,*) 'element number:                       ', elemNum
       write(*,*) 'the thermal material could not be found'
       err_code = 2
       goto 9999
    end if

    ! Get node_ids and node_coords
    do ii=1,20
        node_ids(ii) = fesim%elemente%lsolid20s(jj)%nids(ii)
        node_coords(ii,:) = fesim%knoten%nodes(node_ids(ii))%coords
    end do

    ! Get local tsolid20_C-Matrix
    call sol3_tsolid20_cond(fesim,node_coords,node_ids,lsolid20Num,elemNum,k,tsolid20_C)
    call MatSetValues(C,20,node_ids-1,20,node_ids-1,tsolid20_C,ADD_VALUES,ierr); CHKERRQ(ierr)
 end do

! Finalize global assembly
 call MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
!
! Cleanup
deallocate(d_nnz, o_nnz)
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
