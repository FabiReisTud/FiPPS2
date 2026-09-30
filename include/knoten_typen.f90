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
!> $Id: knoten_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module knoten_typen

  implicit none
  
  type nodes_type
    type(node_type),allocatable     :: nodes(:)
  end type nodes_type

! =================================================================================================
!
! Nodes
!
type node_type
  integer                           :: nid      ! Node-ID                   (Integer>0)
  integer                           :: cid      ! Coordinate System-ID      (Integer>0)
  double precision, dimension(3)    :: coords   ! Coordinates x,y,z         (Real)
end type node_type

! =================================================================================================
!
! MPI-Type Pointer
!
logical                                 :: mpi_type_inited = .false.
integer                                 :: mpi_node_type

contains

! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================

  subroutine knoten_mpi_init()
    
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
    
    TYPE(node_type)                     :: node
    
    !***************************************************************************
    ! node type
    !***************************************************************************

    count = 3												
  
    ALLOCATE(array_of_blockl(1:count)); ALLOCATE(array_of_types(1:count))				
    ALLOCATE(array_of_disp(1:count)); ALLOCATE(address(1:count))
    
    ! nid
    array_of_types (1) = MPI_INTEGER									
    array_of_blockl(1) = 1								
    CALL MPI_GET_ADDRESS(node%nid, address(1), ierror)
    
    ! cid
    array_of_types (2) = MPI_INTEGER									
    array_of_blockl(2) = 1								
    CALL MPI_GET_ADDRESS(node%cid, address(2), ierror)  	
    
    ! coords
    array_of_types (3) = MPI_DOUBLE_PRECISION
    array_of_blockl(3) = 3
    CALL MPI_GET_ADDRESS(node%coords(1), address(3), ierror)
    
    array_of_disp(1) = 0										

    DO ii = 2, count											
       
      array_of_disp(ii) = address(ii) - address(1)							
        												
    END DO

    CALL MPI_TYPE_CREATE_STRUCT(count, array_of_blockl, array_of_disp, array_of_types,  	   &	
    &				mpi_node_type, ierror)
    
    CALL MPI_TYPE_COMMIT(mpi_node_type, ierror) 						
    
    DEALLOCATE(array_of_blockl, array_of_types, array_of_disp, address)
#endif
  end subroutine knoten_mpi_init
  
! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================  
  
  subroutine bcast_knoten(nodes,is_node)
  
    #include "petsc/finclude/petscsys.h"
    use petscsys
    use globale_variablen
  
    implicit none
    
    type(nodes_type)    :: nodes
    logical             :: is_node
    
    PetscMPIInt         :: rank
    PetscErrorCode      :: ierr
    integer             :: node_num
    
    #if !defined (PETSC_HAVE_MPIUNI)

    if (mpi_type_inited .eq. .false.) then
      call knoten_mpi_init()
      mpi_type_inited = .true.
    end if

    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
    
    ! nodes
    if (is_node .eq. .true.) then
      if (rank .eq. 0) node_num = size(nodes%nodes,1) 
      call MPI_Bcast (     node_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) then
        if (allocated(nodes%nodes) .eq. .true.) deallocate(nodes%nodes)
        allocate(nodes%nodes(node_num))
      end if
      call MPI_Bcast (     nodes%nodes, node_num, MPI_NODE_TYPE, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end if
  #endif
  end subroutine bcast_knoten
  
  subroutine free_mem_knoten(knoten)
  
    implicit none
    
    type(nodes_type) :: knoten
    
    if (allocated(knoten%nodes) .eq. .true.) deallocate(knoten%nodes)
  
  end subroutine free_mem_knoten
  
end module knoten_typen
