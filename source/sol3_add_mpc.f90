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
subroutine sol3_add_mpc(fesim,C)
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
#include "petsc/finclude/petscvec.h"
use petscvec
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Interfaces
interface
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
!
! Input + Output
!
Mat,intent(inout)                   :: C !< C
!
! Internal
!
integer                             :: s_nid, mpc_size, max_mpc_size, m_size, m_nid
logical, allocatable                :: s_nid_seen(:) ! Checks, whether a thermal mpc has been added for that s_nid already
integer,allocatable                 :: a_indices(:)
double precision                    :: m_weight, C_max
double precision                    :: pen ! = penalty factor
double precision,allocatable        :: a_vals(:)
PetscScalar, pointer                :: diag(:)
double precision,dimension(:,:),allocatable :: arr_penalty
Vec                                 :: diag_vec
PetscViewer                         :: viewer
real                                :: t0,t1,dt
PetscInt                            :: nrows,ncols

integer                             :: ii,jj,kk,ll


PetscMPIInt                         :: rank
PetscErrorCode                      :: ierr
type(MatInfo)                             :: mat_info(MAT_INFO_SIZE)
!
! =================================================================================================
!
! Initialisation
!
allocate(diag(fesim%num_nodes),s_nid_seen(fesim%num_nodes))
s_nid_seen = .false.
!
! =================================================================================================
!
! Calculation
!
! Get information about the matrix
call MatGetInfo(C, MAT_LOCAL, mat_info, ierr);CHKERRQ(ierr)
!write(*,*) 'Allocated: ', mat_info
write(*,*) 'Allocated: ', mat_info(MAT_INFO_NZ_ALLOCATED)
write(*,*) 'Used:      ', mat_info(MAT_INFO_NZ_USED)
write(*,*) 'Unneeded:  ', mat_info(MAT_INFO_NZ_UNNEEDED)
write(*,*) 'Mallocs:   ', mat_info(MAT_INFO_MALLOCS)
call MatGetSize(C, nrows, ncols, ierr);CHKERRQ(ierr)
write(*,*) 'Nrows:     ', nrows
write(*,*) 'Ncols:     ', ncols

call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'cbeformpc.txt', viewer, ierr);CHKERRQ(ierr)
call MatView(C, viewer, ierr);CHKERRQ(ierr)
call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

! Get penalty factor
pen = 0. ! penalty multiplier
call MatCreateVecs(C, PETSC_NULL_VEC, diag_vec, ierr); CHKERRQ(ierr)
call MatGetDiagonal(C,diag_vec,ierr); CHKERRQ(ierr)
call VecGetArrayF90(diag_vec, diag, ierr); CHKERRQ(ierr)
do ii=1,fesim%num_nodes
    pen = max(pen,abs(diag(ii)))
end do
call VecRestoreArrayF90(diag_vec, diag, ierr); CHKERRQ(ierr)
call VecDestroy(diag_vec, ierr); CHKERRQ(ierr)
pen = pen * 10000 ! Increase by a large number to ensure that penalty is large enough

max_mpc_size = 1
do ii = 1, size(fesim%randbedingungen%mpcs)
    max_mpc_size = max(max_mpc_size, size(fesim%randbedingungen%mpcs(ii)%independend) + 1)
end do
allocate(a_vals(max_mpc_size), a_indices(max_mpc_size), arr_penalty(max_mpc_size, max_mpc_size))


! Find the ids all mpcs with unique dependend nodes
dt = 0.0
do ii=1,size(fesim%randbedingungen%mpcs)
    if (.not. s_nid_seen(fesim%randbedingungen%mpcs(ii)%dependend%nid)) then ! Only add penalties, if an mpc with that s_nid hasnt been processed already
        s_nid_seen(fesim%randbedingungen%mpcs(ii)%dependend%nid) = .true.
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

!        do jj=1,mpc_size
!            do kk=1,mpc_size
!                arr_penalty(jj,kk) = arr_penalty(jj,kk) + pen * a_vals(jj) * a_vals(kk)
!            end do
!        end do
        arr_penalty(1:mpc_size,1:mpc_size) = pen * spread(a_vals(1:mpc_size), dim=2, ncopies=mpc_size) * &
                                                   spread(a_vals(1:mpc_size), dim=1, ncopies=mpc_size)

        call MatSetValues(C, mpc_size, a_indices(1:mpc_size), mpc_size, a_indices(1:mpc_size), arr_penalty(1:mpc_size, 1:mpc_size), ADD_VALUES, ierr); CHKERRQ(ierr)
    end if
end do
deallocate(a_vals,a_indices,arr_penalty)

call MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
call MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
call cpu_time(t1)
print *, 'Elapsed CPU time = ', dt, ' seconds'

! Get information about the matrix
call MatGetInfo(C, MAT_LOCAL, mat_info, ierr);CHKERRQ(ierr)
!write(*,*) 'Allocated: ', mat_info
write(*,*) 'Allocated: ', mat_info(MAT_INFO_NZ_ALLOCATED)
write(*,*) 'Used:      ', mat_info(MAT_INFO_NZ_USED)
write(*,*) 'Unneeded:  ', mat_info(MAT_INFO_NZ_UNNEEDED)
write(*,*) 'Mallocs:   ', mat_info(MAT_INFO_MALLOCS)
call MatGetSize(C, nrows, ncols, ierr);CHKERRQ(ierr)
write(*,*) 'Nrows:     ', nrows
write(*,*) 'Ncols:     ', ncols

call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'caftermpc.txt', viewer, ierr);CHKERRQ(ierr)
call MatView(C, viewer, ierr);CHKERRQ(ierr)
call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

end subroutine
