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
subroutine sol3_get_c_mat(fesim,C)
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
PetscInt                            :: nnodes_per_elem,nnodes_global,mpc_size
PetscInt                            :: d_nz,o_nz
PetscInt, allocatable               :: d_nnz(:), o_nnz(:)
integer                             :: err_code=0
integer                             :: pid, tmidid
type(thermal_curve_type)            :: k
real                                :: t0,t1

!MPC
PetscReal                           :: mpc_pen=0.d0
PetscInt                            :: max_mpc_size
PetscInt,allocatable                :: a_indices(:)
PetscReal,allocatable               :: a_vals(:),arr_penalty(:,:)


logical,allocatable                 :: s_id_seen(:) ! keeps track of whether an s_id has already had its m_ids added to d_nnz

PetscMPIInt                         :: rank
PetscErrorCode                      :: ierr
!
! =================================================================================================
!
! Initialisation
nnodes_per_elem=20
nnodes_global = size(fesim%knoten%nodes)
allocate(d_nnz(nnodes_global), o_nnz(nnodes_global),s_id_seen(nnodes_global))
d_nz  = min(124,nnodes_global) ! Neighborhood of 5x5x5 nodes - 1
d_nz = 124
d_nnz = d_nz
s_id_seen = .FALSE.
!
!
! =================================================================================================
!
! Calculation
!
! call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

! Add mpcs-independents to estimate of mpc-dependents
do ii=1,size(fesim%randbedingungen%mpcs)
    if (.not. s_id_seen(fesim%randbedingungen%mpcs(ii)%dependend%nid)) then
        s_id_seen(fesim%randbedingungen%mpcs(ii)%dependend%nid) = .true.
        mpc_size = size(fesim%randbedingungen%mpcs(ii)%independend) + 1
        d_nnz(fesim%randbedingungen%mpcs(ii)%dependend%nid) = d_nnz(fesim%randbedingungen%mpcs(ii)%dependend%nid) + mpc_size ! Add an entry for each independent in the row of the dependent
        do jj=1,size(fesim%randbedingungen%mpcs(ii)%independend)
            d_nnz(fesim%randbedingungen%mpcs(ii)%independend(jj)%nid) = d_nnz(fesim%randbedingungen%mpcs(ii)%independend(jj)%nid) + mpc_size  ! Add one entry in the row of indep for the dep
        end do
    end if
end do

o_nz  = d_nz  ! effectively ignored for single-process-run
o_nnz = d_nnz ! effectively ignored for single-process-run


! Create global conductivity matrix C
call cpu_time(t0)
call MatCreateAIJ(PETSC_COMM_SELF, &
                  PETSC_DECIDE, PETSC_DECIDE, &
                  nnodes_global, nnodes_global, &
                  d_nz, d_nnz, &
                  o_nz, o_nnz, &
                  C, ierr); CHKERRQ(ierr)
call MatSetOption(C, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr); CHKERRQ(ierr)

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

    ! Get penalty
    mpc_pen = max(mpc_pen,maxval(abs(tsolid20_C)))
    call MatSetValues(C,20,node_ids-1,20,node_ids-1,tsolid20_C,ADD_VALUES,ierr); CHKERRQ(ierr)
 end do

write(*,*) 'Adding the Penalty-MPC-Contributions to the C-Matrix'
!!!!!!!!!!! Add Penalty-MPC-contributions to the C-Vector
mpc_pen = mpc_pen * 10000 ! increase max element stiffness by xxx
s_id_seen = .false.
! Determine max mpc size
max_mpc_size = 1
do ii = 1, size(fesim%randbedingungen%mpcs)
    max_mpc_size = max(max_mpc_size, size(fesim%randbedingungen%mpcs(ii)%independend) + 1)
end do
! Allocate
allocate(a_vals(max_mpc_size), a_indices(max_mpc_size), arr_penalty(max_mpc_size, max_mpc_size))

! Iter through all mpcs
do ii=1,size(fesim%randbedingungen%mpcs)
    if (.not. s_id_seen(fesim%randbedingungen%mpcs(ii)%dependend%nid)) then
        s_id_seen(fesim%randbedingungen%mpcs(ii)%dependend%nid) = .true.
        mpc_size = size(fesim%randbedingungen%mpcs(ii)%independend) + 1
        a_vals(1:mpc_size) = 0.d0
        a_indices(1:mpc_size) = 0
        arr_penalty(1:mpc_size, 1:mpc_size) = 0.d0

        a_vals(1) = fesim%randbedingungen%mpcs(ii)%dependend%fac
        a_indices(1) = fesim%randbedingungen%mpcs(ii)%dependend%nid

        do jj=1,size(fesim%randbedingungen%mpcs(ii)%independend)
            a_indices(jj+1) = fesim%randbedingungen%mpcs(ii)%independend(jj)%nid
            a_vals(jj+1)    = fesim%randbedingungen%mpcs(ii)%independend(jj)%fac
        end do
        a_indices = a_indices - 1 ! Set to 0-based indices for petsc

        arr_penalty(1:mpc_size,1:mpc_size) = mpc_pen *  spread(a_vals(1:mpc_size), dim=2, ncopies=mpc_size) * &
                                                        spread(a_vals(1:mpc_size), dim=1, ncopies=mpc_size)

        call MatSetValues(C, mpc_size, a_indices(1:mpc_size), mpc_size, a_indices(1:mpc_size), arr_penalty(1:mpc_size, 1:mpc_size), ADD_VALUES, ierr); CHKERRQ(ierr)
    end if
end do
deallocate(a_vals,a_indices,arr_penalty)

write(*,*) 'Finished the Penalty-MPC-Contributions to the C-Matrix'


! Finalize global assembly
 call MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
!
! Cleanup
deallocate(d_nnz, o_nnz)

call cpu_time(t1)
print *, 'Elapsed CPU time during C-Assembly = ', t1 - t0, ' seconds'
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
