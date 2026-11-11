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
!> subroutine for computing load-vectors for different load cases
!
!> @details
!> subroutine calculates the resulting load vectors for various kinds of single
!> loads and stores them in an one dimensional array for each loadcase
!> subroutine can calculate loadvectors for multiple loadcases and stores them
!> in consecutive arrays
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit, 07.12.2009
!
!> $Id: load_loadvector.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine load_loadvector (fesim, Fa, Fout, scloop)

use fesimulation_typen
use konstanten
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation)                                             :: fesim
integer, intent(in)                                             :: scloop
!
! Output
!
double precision, dimension(fesim%internals%dim_dof), intent(out) :: fa
double precision, dimension(fesim%num_dof)                        :: Fout
!
! Input + Output
!

!
! inner
!
double precision, dimension(:,:), allocatable                   :: lvec
integer, dimension(:), allocatable                              :: lvec_lcids
integer                                                         :: pos

integer                                                         :: ii,jj,hh,neededLCNum
integer                                                         :: err_code=0

! MPC - DISTRIBUTION
integer                                                         :: indepNid,indepDof,indepInd
double precision                                                :: depForce,indepFac,indepAdd

! Temporal
double precision, dimension(:), allocatable                     :: pElem

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
neededLCNum = 1
allocate(lvec_lcids(neededLCNum))

lvec_lcids(1) = fesim%lasten%subcases(scloop)%loadid
!
! =================================================================================================
!
! Calculate loadvectors for loadcases
!

! allocate and 0's

allocate(lvec(1:fesim%num_dof,1:neededLCNum))

lvec = 0.D0

allocate(pElem(fesim%num_elements))
pElem = 0.d0

!
do hh=1,size(fesim%lasten%loads,1)                                  ! Loop over loads

    do pos = 1, size(lvec_lcids,1)                                  ! Loop over needed LoadCases

        if (fesim%lasten%loads(hh)%lcid .eq. lvec_lcids(pos)) then

            ! Forces

            if (fesim%is_force == .true.) then                      ! if forces exist

                call load_force(fesim,hh,lvec,pos)

            end if

            ! Moments

            if (fesim%is_moment == .true.) then                     ! if moments exist

                call load_moment(fesim,hh,lvec,pos)

            end if

            ! p2loads auf beam2s

            if (fesim%is_p2load == .true.) then                     ! if p2load exists

                call load_p2load(fesim,hh,lvec,pos,pElem)

            end if

            ! p8loads auf quad8s

            if (fesim%is_p8load == .true.) then                     ! if p8load exist

                call load_p8load(fesim,hh,lvec,pos,pElem)

            end if

            ! p20loads auf lsolid20s

            if (fesim%is_p20load == .true.) then                    ! if p20load exists

                call load_p20load(fesim,hh,lvec,pos,pElem)

            end if

        end if

    end do

end do

if (fesim%ausgabe%outputVTK .eq. .true.) then
  call vtkoutcelldata(pElem, fesim%num_elements, "Pressure                      ")
endif

deallocate(pElem)

if ((fesim%is_temperature == .true.) .or. (fesim%is_beam2temp == .true.) .or. (fesim%is_quad8temp == .true.) .or. (fesim%is_lsolid20temp == .true.)) then ! if temperature loads exist
   do pos=1,size(lvec_lcids,1)
      call temperature_loads(fesim,lvec_lcids(pos), pos, lvec)
   end do
end if

if (fesim%lasten%subcases(scloop)%upstress == .true.) then          ! if initial stress shall be applied
    ! Die Kräfte werden auf den ersten Lastfall dazuaddiert
    call initstress_loads(fesim,lvec,1)
end if

!
! build loadvector according to subcases
!
! find out index of load-id in lvec_lcids
do jj=1,size(lvec_lcids,1)
    do hh = 1,fesim%num_dof
        if (lvec(hh,jj) .eq. 0.d0) cycle
        if (fesim%internals%num_dof_vec(hh) .GT. 0) then
            Fa(fesim%internals%num_dof_vec(hh))=Fa(fesim%internals%num_dof_vec(hh))+lvec(hh,jj)
            if ((fesim%ausgabe%outputVTK .eq. .true.) .or. (fesim%ausgabe%outputUser .eq. .true.)) then
                Fout(hh) = Fout(hh) + lvec(hh,jj)
            endif
        else if (fesim%internals%num_dof_vec(hh) .EQ. 0) then
            if ((fesim%ausgabe%outputVTK .eq. .true.) .or. (fesim%ausgabe%outputUser .eq. .true.)) then
                Fout(hh) = Fout(hh) + lvec(hh,jj)
            endif
        else if (fesim%internals%num_dof_vec(hh) .LT. 0) then
            ! OLD IMPLEMENTATION
!            write(*,*) 'Please define no nodal forces for nodes'
!            write(*,*) 'that are dependend nodes of an MPC, exit program'
!            write(*,*) 'NodeID ', hh/ndof
!            Hier den lvec verteilen
!            err_code=1
!            goto 9999

            ! DISTRIBUTE
            depForce = lvec(hh,jj) ! Get force at dependend node
!            write(*,*) '======================================='
!            write(*,*) '      hh=',hh
!            write(*,*) 'fesim%internals%num_dof_vec(',hh,')=',fesim%internals%num_dof_vec(hh)
!            write(*,*) 'depForce=',depForce
            if ((fesim%ausgabe%outputVTK .eq. .true.) .or. (fesim%ausgabe%outputUser .eq. .true.)) then
                Fout(hh) = Fout(hh) + depForce
            endif

            do ii=1,size(fesim%randbedingungen%mpcs(-fesim%internals%num_dof_vec(hh))%independend)
                ! Get values of independend node
                indepNid = fesim%randbedingungen%mpcs(-fesim%internals%num_dof_vec(hh))%independend(ii)%nid
!                write(*,*) 'indepNid=',indepNid
                indepFac = fesim%randbedingungen%mpcs(-fesim%internals%num_dof_vec(hh))%independend(ii)%fac
!                write(*,*) 'indepFac=',indepFac
                indepDof = fesim%randbedingungen%mpcs(-fesim%internals%num_dof_vec(hh))%independend(ii)%dof
!                write(*,*) 'indepDof=',indepDof
                ! Calculate
                indepInd = (indepNid-1)*ndof + indepDof ! Index of Independend-Node-Dof-Entry
!                write(*,*) 'indepInd=',indepInd
                indepAdd = depForce * indepFac          ! Force to be added to independed dof
!                write(*,*) 'indepAdd=',indepAdd
                ! Add
!                write(*,*) 'Fa(indepInd)=',Fa(indepInd),' (before adding)'
                Fa(indepInd) = Fa(indepInd) + indepAdd
!                write(*,*) 'Fa(indepInd)=',Fa(indepInd),' (after adding)'
                if ((fesim%ausgabe%outputVTK .eq. .true.) .or. (fesim%ausgabe%outputUser .eq. .true.)) then
                    Fout(indepInd) = Fout(indepInd) + indepAdd
                end if
            end do
        end if
    end do

!    write(*,*)''
!    write(*,*)'================================='
!    do hh=1,size(Fa)
!        if (Fa(hh)*Fa(hh) .GT. 0.01) then ! Nur über 0.1 N ausschreiben
!            write(*,*)'hh=',hh,' Fa=',Fa(hh),' num_dof_vec=',fesim%internals%num_dof_vec(hh)
!        end if
!    end do
!    write(*,*)'================================='
!    do hh=1,size(Fout)
!        if (Fout(hh)*Fout(hh) .GT. 0.01) then ! Nur über 0.1 N ausschreiben
!            write(*,*)'hh=',hh,' Fout=',Fout(hh),' num_dof_vec=',fesim%internals%num_dof_vec(hh)
!        end if
!    end do
!    write(*,*)'================================='

end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,'(A)',advance='YES')    'load_loadvector'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return

end subroutine load_loadvector
