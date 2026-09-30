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
!> @brief
!
!
!> @details
!
!> @author
!
!
!> $Id: element_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module element_typen

  implicit none
  
type elemente_type
  type (beam2_type), allocatable      :: beam2s (:)
  type (quad8_type), allocatable      :: quad8s (:)
  type (Lsolid20_type), allocatable   :: lsolid20s(:)
end type elemente_type

!
!--------------------------------------------------------------------------------------------------
!
! BEAM2
!
type beam2_type
  integer                       :: eid      ! Element-ID            (Integer>0)
  integer                       :: pid      ! Property-ID           (Integer>0)
  integer, dimension(1:2)       :: nids     ! Node-IDs              (Integer>0)
  double precision, dimension(3):: xi       ! Components of orientation vector at node 1 in local coordinates (Real)
  integer                       :: n0       ! Node for orientation vector, vector from node 1 to g0 (Integer>0, n0.NE.nid1, n0.NE.nid2)
                                            ! either x1,x2,x3 or n0 can be specified; x1,x2,x3 is default method
                                            ! n0 method not implemented yet (22.06.2010)
! for internal usage
  integer                       :: int_pid  ! Interne Property-Nummer (entspricht dem Index im Property-Array)
end type beam2_type

!
!--------------------------------------------------------------------------------------------------
!
! QUAD8
!
type quad8_type
  integer                       :: eid      ! Element-ID            (Integer>0)
  integer                       :: pid      ! Property-ID           (Integer>0)
  integer, dimension(1:8)       :: nids     ! Node-IDs              (Integer>0)
                                            ! 1-4 Eckknoten
                                            ! 5-8 Knoten Seitenmitte
  double precision	            :: theta    ! Material angle        (Real)
  double precision              :: offset   ! Material offset from FE reference plane  (Real)
  
! for internal usage
  integer                       :: propType ! Property-Type: 1-pshell; 2-pcomp
  integer                       :: int_pid
  double precision              :: area
end type quad8_type

!
!--------------------------------------------------------------------------------
!
! Layered Solid with 20 Nodes
!
type Lsolid20_type
  integer                       :: eid      ! Element-ID   (Integer>0)
  integer                       :: pid      ! Property-ID  (Integer>0)
  integer, dimension(1:20)      :: nids     ! Node-IDs     (Integer>0)
                                            ! 1-8 Vertex nodes
                                            ! 9-20 Midside nodes
! for internal usage
  integer                       :: int_pid
end type Lsolid20_type

logical                                 :: mpi_type_inited = .false.
integer                                 :: mpi_beam2_type
integer                                 :: mpi_quad8_type
integer                                 :: mpi_lsolid20_type

contains

  subroutine element_mpi_init()

#include "petsc/finclude/petscsys.h"

#if !defined (PETSC_HAVE_MPIUNI)

    use petscsys
  
    implicit none
    
    INTEGER                             :: count
    INTEGER, dimension(:), allocatable  :: array_of_blockl
    INTEGER, dimension(:), allocatable  :: array_of_types
    INTEGER(kind=MPI_ADDRESS_KIND), dimension(:), allocatable   :: array_of_disp  
    INTEGER(kind=MPI_ADDRESS_KIND), dimension(:), allocatable   :: address
    
    INTEGER                             :: ii
    INTEGER                             :: ierror
    
    TYPE(beam2_type)                     :: beam2
    TYPE(quad8_type)                    :: quad8
    TYPE(lsolid20_type)                 :: lsolid20
    
    !***************************************************************************
    ! beam2 type
    !***************************************************************************

    count = 6
  
    ALLOCATE(array_of_blockl(1:count)); ALLOCATE(array_of_types(1:count))				
    ALLOCATE(array_of_disp(1:count)); ALLOCATE(address(1:count))
    
    ! eid
    array_of_types (1) = MPI_INTEGER
    array_of_blockl(1) = 1
    CALL MPI_GET_ADDRESS(beam2%eid, address(1), ierror) 	
    
    ! pid
    array_of_types (2) = MPI_INTEGER
    array_of_blockl(2) = 1
    CALL MPI_GET_ADDRESS(beam2%pid, address(2), ierror)
    
    ! nids
    array_of_types (3) = MPI_INTEGER
    array_of_blockl(3) = 2
    CALL MPI_GET_ADDRESS(beam2%nids(1), address(3), ierror)
    
    ! theta
    array_of_types (4) = MPI_DOUBLE_PRECISION
    array_of_blockl(4) = 3
    CALL MPI_GET_ADDRESS(beam2%xi, address(4), ierror)
    
    ! offset
    array_of_types (5) = MPI_INTEGER
    array_of_blockl(5) = 1
    CALL MPI_GET_ADDRESS(beam2%n0, address(5), ierror)
    
    ! int_pid
    array_of_types (6) = MPI_INTEGER
    array_of_blockl(6) = 1
    CALL MPI_GET_ADDRESS(beam2%int_pid, address(6), ierror)
    
    array_of_disp(1) = 0

    DO ii = 2, count
       
      array_of_disp(ii) = address(ii) - address(1)

    END DO

    CALL MPI_TYPE_CREATE_STRUCT(count, array_of_blockl, array_of_disp, array_of_types, &
    &                           mpi_beam2_type, ierror)
    
    CALL MPI_TYPE_COMMIT(mpi_beam2_type, ierror)							
    
    DEALLOCATE(array_of_blockl, array_of_types, array_of_disp, address)
    
    !***************************************************************************
    ! quad8 type
    !***************************************************************************

    count = 8												
  
    ALLOCATE(array_of_blockl(1:count)); ALLOCATE(array_of_types(1:count))				
    ALLOCATE(array_of_disp(1:count)); ALLOCATE(address(1:count))
    
    ! eid
    array_of_types (1) = MPI_INTEGER
    array_of_blockl(1) = 1
    CALL MPI_GET_ADDRESS(quad8%eid, address(1), ierror) 	
    
    ! pid
    array_of_types (2) = MPI_INTEGER
    array_of_blockl(2) = 1
    CALL MPI_GET_ADDRESS(quad8%pid, address(2), ierror)
    
    ! nids
    array_of_types (3) = MPI_INTEGER
    array_of_blockl(3) = 8
    CALL MPI_GET_ADDRESS(quad8%nids(1), address(3), ierror)
    
    ! theta
    array_of_types (4) = MPI_DOUBLE_PRECISION
    array_of_blockl(4) = 1
    CALL MPI_GET_ADDRESS(quad8%theta, address(4), ierror)
    
    ! offset
    array_of_types (5) = MPI_DOUBLE_PRECISION
    array_of_blockl(5) = 1
    CALL MPI_GET_ADDRESS(quad8%offset, address(5), ierror)
    
    ! propType
    array_of_types (6) = MPI_INTEGER
    array_of_blockl(6) = 1
    CALL MPI_GET_ADDRESS(quad8%propType, address(6), ierror)
    
    ! int_pid
    array_of_types (7) = MPI_INTEGER
    array_of_blockl(7) = 1
    CALL MPI_GET_ADDRESS(quad8%int_pid, address(7), ierror)
    
    ! area
    array_of_types (8) = MPI_DOUBLE_PRECISION
    array_of_blockl(8) = 1
    CALL MPI_GET_ADDRESS(quad8%area, address(8), ierror)
    
    array_of_disp(1) = 0										

    DO ii = 2, count											
       
      array_of_disp(ii) = address(ii) - address(1)							
        												
    END DO

    CALL MPI_TYPE_CREATE_STRUCT(count, array_of_blockl, array_of_disp, array_of_types,  	   &	
    &				mpi_quad8_type, ierror)
    
    CALL MPI_TYPE_COMMIT(mpi_quad8_type, ierror)							
    
    DEALLOCATE(array_of_blockl, array_of_types, array_of_disp, address)

    !***************************************************************************
    ! lsolid20 type
    !***************************************************************************

    count = 4                                               
  
    ALLOCATE(array_of_blockl(1:count)); ALLOCATE(array_of_types(1:count))               
    ALLOCATE(array_of_disp(1:count)); ALLOCATE(address(1:count))
    
    ! eid
    array_of_types (1) = MPI_INTEGER
    array_of_blockl(1) = 1
    CALL MPI_GET_ADDRESS(lsolid20%eid, address(1), ierror)  
    
    ! pid
    array_of_types (2) = MPI_INTEGER
    array_of_blockl(2) = 1
    CALL MPI_GET_ADDRESS(lsolid20%pid, address(2), ierror)
    
    ! nids
    array_of_types (3) = MPI_INTEGER
    array_of_blockl(3) = 20
    CALL MPI_GET_ADDRESS(lsolid20%nids(1), address(3), ierror)
    
    ! int_pid
    array_of_types (4) = MPI_INTEGER
    array_of_blockl(4) = 1
    CALL MPI_GET_ADDRESS(lsolid20%int_pid, address(4), ierror)
    
    array_of_disp(1) = 0                                        

    DO ii = 2, count                                            
       
      array_of_disp(ii) = address(ii) - address(1)                          
                                                        
    END DO

    CALL MPI_TYPE_CREATE_STRUCT(count, array_of_blockl, array_of_disp, array_of_types,         &    
    &               mpi_lsolid20_type, ierror)
    
    CALL MPI_TYPE_COMMIT(mpi_lsolid20_type, ierror)                         
    
    DEALLOCATE(array_of_blockl, array_of_types, array_of_disp, address)
#endif
  end subroutine element_mpi_init

  subroutine bcast_element_typen(elems,is_beam2,is_quad8,is_lsolid20)
  
    #include "petsc/finclude/petscsys.h"
    use petscsys
    use globale_variablen
  
    implicit none
    
    type(elemente_type) :: elems
    logical             :: is_beam2
    logical             :: is_quad8
    logical             :: is_lsolid20
    
    PetscMPIInt         :: rank
    PetscErrorCode      :: ierr
    integer             :: beam2_num, quad8_num, lsolid20_num
    
    #if !defined (PETSC_HAVE_MPIUNI)

    if (mpi_type_inited .eq. .false.) then
      call element_mpi_init()
      mpi_type_inited = .true.
    end if
                    
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
    
    ! beam2
    if (is_beam2 .eq. .true.) then
      if (rank .eq. 0) beam2_num = size(elems%beam2s,1) 
      call MPI_Bcast (beam2_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(elems%beam2s(beam2_num))
      call MPI_Bcast (elems%beam2s, beam2_num, MPI_BEAM2_TYPE, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end if

    ! quad8s
    if (is_quad8 .eq. .true.) then
      if (rank .eq. 0) quad8_num = size(elems%quad8s,1) 
      call MPI_Bcast (quad8_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(elems%quad8s(quad8_num))
      call MPI_Bcast (elems%quad8s, quad8_num, MPI_QUAD8_TYPE, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end if
    
    ! lsolid20
    if (is_lsolid20 .eq. .true.) then
      if (rank .eq. 0) lsolid20_num = size(elems%lsolid20s,1) 
      call MPI_Bcast (lsolid20_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(elems%lsolid20s(lsolid20_num))
      call MPI_Bcast (elems%lsolid20s, lsolid20_num, MPI_LSOLID20_TYPE, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end if
    #endif
  end subroutine bcast_element_typen
  
  subroutine free_mem_element(elemente)
  
    implicit none
    
    type(elemente_type) :: elemente
    
    if (allocated(elemente%quad8s)    .eq. .true.) deallocate(elemente%quad8s)
    if (allocated(elemente%beam2s)    .eq. .true.) deallocate(elemente%beam2s)
    if (allocated(elemente%lsolid20s) .eq. .true.) deallocate(elemente%lsolid20s)
  
  end subroutine free_mem_element
    
end module element_typen
