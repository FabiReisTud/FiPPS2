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
PetscInt, dimension(20)             :: node_ids
double precision, dimension(20,3)   :: node_coords
double precision, dimension(20,20)  :: tsolid20_C

PetscInt                            :: ii,jj, kk,ll,e,lsolid20Num,elemNum,nnodes_per_elem
PetscInt                            :: err_code=0
PetscInt                            :: pid, tmidid
type(thermal_curve_type)            :: k

PetscMPIInt                         :: rank
PetscErrorCode                      :: ierr
!
! =================================================================================================
!
! Initialisation
nnodes_per_elem=20
!
!
! =================================================================================================
!
! Calculation
!
 call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

! Create global conductivity matrix C
 call MatCreateAIJ( PETSC_COMM_WORLD,               &
                    PETSC_DECIDE, PETSC_DECIDE,     &
                    nglobal, nglobal,               &
                    nnodes_per_elem, PETSC_NULL_INTEGER, &
                    nnodes_per_elem, PETSC_NULL_INTEGER, &
                    C, ierr)

! Iter through all the elements
 do jj = 1,size(fesim%elemente%lsolid20s,1)
    if (mpi_lsolid20_proc_dist(jj) .NE. rank) cycle

    pid = fesim%elemente%lsolid20s(jj)%pid ! PropertyID des Elements
    tmidid = 0 ! ID des Thermischen Materials initialisieren
    ! Wärmeleitkoeffizientenkurve holen
    do kk=1,size(fesim%thermal%propthermals)
        if (fesim%thermal%propthermals(kk)%propid == pid) then
            tmidid = fesim%thermal%propthermals(kk)%tmidid ! MatThermalTd gefunden! --> TMIDID speichern!
            do ll=1,size(fesim%thermal%materials_td)
                if (fesim%thermal%materials_td(ll)%mid == tmidid) then
                    k = fesim%thermal%materials_td(ll)%k ! k gefunden! --> Speichern!
                    exit
                end if
            end do
            if (tmidid .NE. 0) then
                exit
            end if
        end if
    end do

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
        ll = 0
        do kk=1,size(fesim%knoten%nodes)
            if (fesim%knoten%nodes(kk)%nid == node_ids(ii)) then
                ll = kk
                exit
            end if
        end do

        if (ll == 0) then
            write(*,*) 'Could not find NID', node_ids(ii), 'in list of global nodes fesim%knoten%nodes.'
            err_code = 2
            goto 9999
        end if

        node_coords(ii,:) = fesim%knoten%nodes(ll)%coords
    end do

    ! Get local tsolid20_C-Matrix
    call sol3_tsolid20_cond(fesim,node_coords,node_ids,lsolid20Num,elemNum,k,tsolid20_C)
    ! Assemble
    call MatSetValues(  C,20,node_ids-1,20,node_ids-1,lsolid20_C(1:20,1:20,e),ADD_VALUES,ierr)
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
