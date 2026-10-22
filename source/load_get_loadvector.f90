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
!> @author Martin Rädel, TU Dresden, Diplomarbeit, 08.12.2009
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter, 09.01.2010
!
!> $Id: load_get_loadvector.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine load_get_loadvector(fesim, FaS, Fout, scloop)

  use konstanten
  use fesimulation_typen
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscvec.h"
  use petscvec
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
  type(fe_simulation)                           :: fesim
  double precision, dimension(:), allocatable   :: Fa           ! Kraftvektor ohne gesperrte Freiheiten
  double precision, dimension(fesim%num_dof)    :: Fout         ! gesperrter Freiheitenvektor
  integer                                       :: ii

  PetscInt                                      :: pii
  PetscErrorCode                                :: ierr
  PetscScalar                                   :: v
  Vec                                           :: FaS          ! Kraftvektor ohne gesperrte Freiheiten als Sparse-Vektor
  integer                                       :: scloop

  integer                                       :: err_code=0
!
! =================================================================================================
!
! Initialisation
!

  allocate(Fa(fesim%internals%dim_dof))
  Fa = 0.d0
!
! =================================================================================================
!
! Calculation

  Fout = 0.d0

  if (fesim%is_load == .true.) then

    call load_loadvector (fesim, Fa, Fout, scloop)

  end if

  if ((fesim%ausgabe%outputVTK .eq. .true.) .or. (fesim%ausgabe%outputUser .eq. .true.)) then
    do ii = 1,fesim%num_nodes
        if (fesim%knoten%nodes(ii)%cid .ne. 0) then
            Fout((ii-1)*ndof+1:(ii-1)*ndof+3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat,Fout((ii-1)*ndof+1:(ii-1)*ndof+3))
            Fout((ii-1)*ndof+4:(ii-1)*ndof+6) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat,Fout((ii-1)*ndof+4:(ii-1)*ndof+6))
        end if
    end do
    call vtkoutstatic(fesim, Fout, 'Kraftrandbed', .false.)
    call vtkoutstatic(fesim, Fout, 'Momentrandbe', .true.)
  endif
  !
  ! =================================================================================================

  ! convert Fa to sparse
  do ii = 1, fesim%internals%dim_dof
    if (Fa(ii) .NE. 0.d0) then
      v = Fa(ii)
      pii = ii-1
      call VecSetValue(FaS, pii, v, INSERT_VALUES, ierr); CHKERRQ(ierr)
    end if
  end do

  deallocate(Fa)
!
! =================================================================================================
!
! Error handling
!
9999 continue

  if (err_code /= 0) then

    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'load_get_loadvector'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop

  end if

end subroutine load_get_loadvector
