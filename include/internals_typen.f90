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
!
!> $Id: internals_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
module internals_typen

implicit none

type internals_type

  ! Hilfsgrößen, um Aufwand zu sparen
  integer                               :: dim_dof              ! Anzahl der tatsächlichen Freiheiten (ohne gesperrte)
  integer, dimension(:), allocatable    :: num_dof_vec
  logical                               :: failed               ! Zusatzflag, wenn bereits etwas versagt ist, können damit nachfolgende Berechnungen abgeschaltet werden
  type(meshcoupling_type), dimension(:), allocatable :: aeroElem2structNode
  type(meshcoupling_type), dimension(:), allocatable :: structElem2aeroNode
  integer, dimension(:), allocatable    :: structQuad8IDs
  integer, dimension(:), allocatable    :: structBeam2IDs
  double precision, dimension(:,:), allocatable :: aeroDispl
  integer                               :: indexBeforeAeroP2Load
  integer                               :: indexBeforeAeroP8Load
  integer                               :: indexBeforeAeroLoad
  integer                               :: highestLID
  integer                               :: highestLIDtemp
  double precision,dimension(:,:),allocatable :: thermal_dof_map    ! Array für alle thermischen Randbedingung, Format: (NID,TYPE,VAL1,VAL2)

end type internals_type

type meshcoupling_type

  integer                               :: nodeID               ! Knoten, auf den die Informationen übertragen werden soll (muss im Element liegen)
  integer                               :: elemID               ! Element, von dem die Informationen übernommen werden soll
  double precision                      :: xi                   ! xi-Koordinate des Knotens im Element
  double precision                      :: eta                  ! eta-Koordinate des Knotens im Element

end type meshcoupling_type

contains

  subroutine bcast_internals(internals,num_dof)

#include "petsc/finclude/petscsys.h"
    use petscsys

    implicit none

    PetscMPIInt     :: rank
    PetscErrorCode  :: ierr
    type(internals_type) :: internals
    INTEGER         :: num_dof

#if !defined (PETSC_HAVE_MPIUNI)

    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

    call MPI_Bcast (internals%dim_dof    ,       1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    if (rank .ne. 0) then
        allocate ( internals%num_dof_vec(num_dof) )
    end if
    call MPI_Bcast (internals%num_dof_vec, num_dof, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (internals%failed     ,       1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

#endif

  end subroutine bcast_internals

  subroutine free_mem_internals(internals)

    implicit none

    type(internals_type) :: internals

    if (allocated(internals%num_dof_vec) .eq. .true.) deallocate(internals%num_dof_vec)
    if (allocated(internals%aeroElem2structNode) .eq. .true.) deallocate(internals%aeroElem2structNode)
    if (allocated(internals%structElem2aeroNode) .eq. .true.) deallocate(internals%structElem2aeroNode)
    if (allocated(internals%structQuad8IDs) .eq. .true.) deallocate(internals%structQuad8IDs)
    if (allocated(internals%structBeam2IDs) .eq. .true.) deallocate(internals%structBeam2IDs)
    if (allocated(internals%aeroDispl) .eq. .true.) deallocate(internals%aeroDispl)

  end subroutine free_mem_internals

end module internals_typen
