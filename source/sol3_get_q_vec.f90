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
!> subroutine for creating the global thermal flux vector q
!
!> @details
!> subroutine for creating the global thermal flux vector q
!
!> @author Fabian Reisener, TU Dresden, WM 2026
!
! =================================================================================================
subroutine sol3_get_q_vec(fesim,q,qf,qc,is_free,is_bc)
!=================================================================================================
use fesimulation_typen
#include "petsc/finclude/petscvec.h"
use petscvec
#include "petsc/finclude/petscmat.h"
use petscmat
#include "petsc/finclude/petscis.h"
use petscis

use globale_variablen
use konstanten
use integration_schemes

implicit none

! Input
type(fe_simulation), intent(in) :: fesim
IS                              :: is_free,is_bc
! Output
Vec, intent(out)                :: q,qf,qc
! Internal
PetscErrorCode                  :: ierr
PetscMPIInt                     :: rank
integer                         :: ii, jj,ip,np, node_id
integer                         :: nsize
double precision                :: q_flux,q_local(8)
double precision, allocatable   :: q_array(:)
double precision                :: T_node
integer, allocatable            :: q_indices(:)
integer, dimension(8)           :: node_ids
PetscInt                        :: is_bc_size, is_free_size
! Internal - Numerical Integration
double precision                :: w_ip,xi_ip(4),eta_ip(4),w_xi(4),w_eta(4)
double precision                :: Ni(8),dNidxi(8),dNideta(8),detJ
double precision                :: node_coords(8,3)

interface
    subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta,d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)

        double precision, intent(in) :: xi, eta

        double precision, dimension(8), optional, intent(out) :: Ni
        double precision, dimension(8), optional, intent(out) :: dNidxi, dNideta
        double precision, dimension(8), optional, intent(out) :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
    end subroutine quad8_ansatzfunction_xieta
end interface
! =================================================================================================
!
! Initialization
!
call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)

! Allocate q vector with size = number of thermal DOFs (number of nodes for temperature)
nsize = fesim%num_nodes
call VecCreate(PETSC_COMM_WORLD, q, ierr); CHKERRQ(ierr)
call VecSetSizes(q, PETSC_DECIDE, nsize, ierr); CHKERRQ(ierr)
call VecSetFromOptions(q, ierr); CHKERRQ(ierr)
call VecSet(q, 0.0d0, ierr); CHKERRQ(ierr)

! Allocate temporary arrays for setting values
allocate(q_array(nsize))
allocate(q_indices(nsize))

! Initialize arrays
q_array = 0.0d0
do ii = 1, nsize
    q_indices(ii) = ii - 1  ! PETSc zero-based indexing
end do

! =================================================================================================
!
! Assemble flux BC contributions
!
do ii = 1, size(fesim%thermal%fluxes)
    q_flux = fesim%thermal%fluxes(ii)%q
    node_ids(1) = fesim%thermal%fluxes(ii)%nid1
    node_ids(2) = fesim%thermal%fluxes(ii)%nid2
    node_ids(3) = fesim%thermal%fluxes(ii)%nid3
    node_ids(4) = fesim%thermal%fluxes(ii)%nid4
    node_ids(5) = fesim%thermal%fluxes(ii)%nid5
    node_ids(6) = fesim%thermal%fluxes(ii)%nid6
    node_ids(7) = fesim%thermal%fluxes(ii)%nid7
    node_ids(8) = fesim%thermal%fluxes(ii)%nid8

    !!! Compute local flux contributions
    q_local = 0.d0 ! Init the local q-contribution-vector
    ! Get coordinates of the eight surface nodes
    do jj = 1, 8
        node_coords(jj,:) = fesim%knoten%nodes(node_ids(jj))%coords
    end do

    call integration_points_2d( 2, "Gauss",   xi_ip, w_xi, 2, "Gauss",  eta_ip, w_eta)
    do ip=1,4   ! Iter through Integration Points
        w_ip = w_xi(ip)*w_eta(ip)
        call quad8_ansatzfunction_xieta(xi_ip(ip), eta_ip(ip), Ni, dNidxi,dNideta)
        detJ = j_s(xi_ip(ip),eta_ip(ip),node_coords) ! area scaling factor
        q_local = q_local + Ni * q_flux * detJ * w_ip
    end do

    do jj = 1, 8
        q_array(node_ids(jj)) = q_array(node_ids(jj)) + q_local(jj)
    end do
end do

! =================================================================================================
!
! Assemble convection BC contributions
!
do ii = 1, size(fesim%thermal%convections)
    node_ids(1) = fesim%thermal%convections(ii)%nid1
    node_ids(2) = fesim%thermal%convections(ii)%nid2
    node_ids(3) = fesim%thermal%convections(ii)%nid3
    node_ids(4) = fesim%thermal%convections(ii)%nid4
    node_ids(5) = fesim%thermal%convections(ii)%nid5
    node_ids(6) = fesim%thermal%convections(ii)%nid6
    node_ids(7) = fesim%thermal%convections(ii)%nid7
    node_ids(8) = fesim%thermal%convections(ii)%nid8

    ! Compute local flux contributions
    q_local = 0.d0

    do jj = 1, 8
        q_array(node_ids(jj)) = q_array(node_ids(jj)) + q_local(jj)
    end do
end do

! =================================================================================================
!
! Set values into PETSc vector q
!
call VecSetValues(q, nsize, q_indices, q_array, INSERT_VALUES, ierr); CHKERRQ(ierr)

! Assemble vector
call VecAssemblyBegin(q, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(q, ierr); CHKERRQ(ierr)

! Get Qf
call ISGetSize(is_free,is_free_size,ierr); CHKERRQ(ierr)
call VecCreate(PETSC_COMM_SELF,qf,ierr); CHKERRQ(ierr)
call VecSetSizes(qf, PETSC_DECIDE, is_free_size, ierr); CHKERRQ(ierr)
call VecSetFromOptions(qf, ierr); CHKERRQ(ierr)
call VecGetSubvector(q,is_free,qf,ierr); CHKERRQ(ierr)
call VecAssemblyBegin(qf, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(qf, ierr); CHKERRQ(ierr)
!call VecView(qf, PETSC_VIEWER_STDOUT_SELF, ierr); CHKERRQ(ierr)

! Get Qc
call ISGetSize(is_bc,is_bc_size,ierr); CHKERRQ(ierr)
call VecCreate(PETSC_COMM_SELF,qc,ierr); CHKERRQ(ierr)
call VecSetSizes(qc, PETSC_DECIDE, is_bc_size, ierr); CHKERRQ(ierr)
call VecSetFromOptions(qc, ierr); CHKERRQ(ierr)
call VecGetSubvector(q,is_bc,qc,ierr); CHKERRQ(ierr)
call VecAssemblyBegin(qc, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(qc, ierr); CHKERRQ(ierr)
!call VecView(qc, PETSC_VIEWER_STDOUT_SELF, ierr); CHKERRQ(ierr)

! Deallocate temporaries
deallocate(q_array)
deallocate(q_indices)

return

contains

function j_s(xi, eta, coords) result(detJ)
    implicit none

    double precision, intent(in) :: xi, eta
    double precision, intent(in) :: coords(8,3)
    double precision :: detJ

    double precision :: Ni(8)
    double precision :: dNidxi(8), dNideta(8)
    double precision :: dx_dxi(3), dx_deta(3)
    double precision :: cross(3)

    call quad8_ansatzfunction_xieta( &
        xi, eta, Ni, dNidxi, dNideta)

    ! Tangent vectors of the surface
    dx_dxi  = matmul(dNidxi,  coords)
    dx_deta = matmul(dNideta, coords)

    ! Cross product
    cross(1) = dx_dxi(2)*dx_deta(3) &
             - dx_dxi(3)*dx_deta(2)

    cross(2) = dx_dxi(3)*dx_deta(1) &
             - dx_dxi(1)*dx_deta(3)

    cross(3) = dx_dxi(1)*dx_deta(2) &
             - dx_dxi(2)*dx_deta(1)

    ! Surface Jacobian determinant (area scaling factor)
    detJ = sqrt(sum(cross**2))

end function j_s


end subroutine sol3_get_q_vec
