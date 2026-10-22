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
!> \file assemble_aa.f90
! =================================================================================================
!
!>  \brief      subroutine for assembling element stiffness matrix to aa-Part of the total
!>              system matrix
!
!>  \details    subroutine assembles element stiffness matrix of shell element to aa-Part
!>              of the total system stiffness matrix
!
!>  \author     Martin Rädel, TU Dresden, Diplomarbeit  08.12.2009
!>  \author     Andreas Hauffe,TU Dresden, WM           11.10.2012
!
!> $Id: assemble_aa.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
!
! =================================================================================================

subroutine assemble_aa (fesim, num_nodes_el, node_ids, Kel, KaaS)
!
use konstanten, ONLY : ndof
use fesimulation_typen

#include "petsc/finclude/petscmat.h"
use petscmat


!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
type(fe_simulation)                                                             :: fesim
integer, intent(in)                                                             :: num_nodes_el !< number of element nodes
integer, dimension(num_nodes_el), intent(in)                                    :: node_ids     !< IDs of the element nodes
double precision,dimension(num_nodes_el*ndof,num_nodes_el*ndof),intent(inout)   :: Kel          !< element stiffness matrix
!
! Output
!

!
! Input + Output
!
Mat, intent(inout)                                                              :: KaaS         !< AA-Part of of the global stiffness matrix (sparse - PETSC)
!
! inner
!
integer                                                                         :: ii,jj,iii,jjj,mm,nn
integer, dimension(num_nodes_el*ndof)                                           :: index
integer                                                                         :: err_code=0

PetscInt                                                                        :: pii,pjj
PetscErrorCode                                                                  :: ierr
PetscScalar                                                                     :: v

logical                                                                         :: transNodeCoords
double precision, dimension(num_nodes_el*ndof,num_nodes_el*ndof)                :: transMat
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
! upper triangle only
!
  transNodeCoords = .false.
  do ii = 1,num_nodes_el
    if (fesim%knoten%nodes(node_ids(ii))%cid .ne. 0) then ! check if any node has a custom Coordinate System
      transNodeCoords = .true.
      exit
    end if
  end do

  if (transNodeCoords .eq. .true.) then
    call nodalCoordSysTransMat(fesim,num_nodes_el,node_ids,transMat)    ! if custom coordinate system,
    Kel = matmul(matmul(transpose(transMat),Kel),transMat)              !  transform the element stiffness matrix to that coordinate system
  end if

! Determine whether the node part of spc or mpc by looking it up in the num_dof_vec
! index>0 --> normal equation number
! index=0 --> spc
! index<0 --> mpc
do ii = 1,num_nodes_el
  do jj = 1,ndof
    index((ii-1)*ndof+jj) = fesim%internals%num_dof_vec((node_ids(ii)-1)*ndof+jj)
  end do
end do

! Assemble Matrix
do ii = 1,size(Kel,1)
  iii = index(ii)
  ! if first dof is not an SPC
  if ( iii .NE. 0) then
    do jj = 1, size(Kel,1)
      jjj = index(jj)
      ! if second dof is not an SPC
      if (jjj .NE. 0) then
        ! if matrix entry is not equal zero
        if (Kel(ii,jj) .NE. 0.d0) then
          ! if at least one dof is an dependend of an MPC
          if (iii .LT. 0 .OR. jjj .LT. 0) then
            ! if only the first dof is a dependend dof of an MPC
            if (iii .LT. 0 .AND. jjj .GT. 0) then
              pjj = jjj-1
              do mm = 1,size(fesim%randbedingungen%mpcs(-iii)%independend,1)
                if (fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof .gt. 0) then !Freiheit ist durch SPC gesperrt und bereits aus dem Gleichungssystem ausgebaut
                  pii = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-iii)%independend(mm)%nid-1)*ndof+fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof)-1
                  v = Kel(ii,jj)*fesim%randbedingungen%mpcs(-iii)%independend(mm)%fac
                  if (pjj >= pii) then
                    call MatSetValue(KaaS, pii, pjj, v, ADD_VALUES, ierr); CHKERRQ(ierr)
                  end if
                end if
              end do
            ! if only the second dof is a dependend dof of an MPC
            else if (iii .GT. 0 .AND. jjj .LT. 0) then
              pii = iii-1
              do nn = 1,size(fesim%randbedingungen%mpcs(-jjj)%independend,1)
                if (fesim%randbedingungen%mpcs(-jjj)%independend(nn)%dof .gt. 0) then !Freiheit ist durch SPC gesperrt und bereits aus dem Gleichungssystem ausgebaut
                  pjj = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-jjj)%independend(nn)%nid-1)*ndof+fesim%randbedingungen%mpcs(-jjj)%independend(nn)%dof)-1
                  v = Kel(ii,jj)*fesim%randbedingungen%mpcs(-jjj)%independend(nn)%fac
                  if (pjj >= pii) then
                    call MatSetValue(KaaS, pii, pjj, v, ADD_VALUES, ierr); CHKERRQ(ierr)
                  end if
                end if
              end do
            ! if both dofs are dependend dofs of an MPC
            else
              do mm = 1,size(fesim%randbedingungen%mpcs(-iii)%independend,1)
                if (fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof .gt. 0) then !Freiheit ist durch SPC gesperrt und bereits aus dem Gleichungssystem ausgebaut
                  pii = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-iii)%independend(mm)%nid-1)*ndof+fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof)-1
                  do nn = 1,size(fesim%randbedingungen%mpcs(-jjj)%independend,1)
                    if (fesim%randbedingungen%mpcs(-jjj)%independend(nn)%dof .gt. 0) then !Freiheit ist durch SPC gesperrt und bereits aus dem Gleichungssystem ausgebaut
                      pjj = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-jjj)%independend(nn)%nid-1)*ndof+fesim%randbedingungen%mpcs(-jjj)%independend(nn)%dof)-1
                      v = Kel(ii,jj)*fesim%randbedingungen%mpcs(-iii)%independend(mm)%fac*fesim%randbedingungen%mpcs(-jjj)%independend(nn)%fac
                      if (pjj >= pii) then
                        call MatSetValue(KaaS, pii, pjj, v, ADD_VALUES, ierr); CHKERRQ(ierr)
                      end if
                    end if
                  end do
                end if
              end do
            end if
          else
            v = Kel(ii,jj)
            pii = iii-1
            pjj = jjj-1
            if (pjj >= pii) then
              call MatSetValue(KaaS, pii, pjj, v, ADD_VALUES, ierr); CHKERRQ(ierr)
            end if
          end if
        end if
      end if
    end do
  end if
end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'assemble_aa'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return

end subroutine assemble_aa
