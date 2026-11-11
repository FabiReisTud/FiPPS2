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
!> control subroutine for inserting multi point/freedom constraints
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 10.10.2012
!
!> $Id: boundary_mpc_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine boundary_mpc_control (fesim, fixedU, scloop)

use konstanten
use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
type(fe_simulation)                             :: fesim
integer, intent(in)                             :: scloop
!
! Output
!
logical, dimension(fesim%num_dof),intent(out)   :: fixedU
!
! inner
!
integer                                         :: gg,ii,jj
integer                                         :: err_code=0

! Loop over MPCs
if (fesim%is_mpcadd == .true.) then
    do gg=1,size(fesim%randbedingungen%mpcadds,1)
        if (fesim%lasten%subcases(scloop)%mpcaddid == fesim%randbedingungen%mpcadds(gg)%scid .or. &
          & fesim%randbedingungen%mpcadds(gg)%scid == -10) then                                         ! interne MPC durch Kontakte

            ! Multi point constraints

            if (fesim%is_mpc == .true.) then
                !do ii=1,size(mpcs,1)
                ii = fesim%randbedingungen%mpcadds(gg)%sid

                if (fesim%randbedingungen%mpcadds(gg)%sid .eq. fesim%randbedingungen%mpcs(ii)%sid) then
                    if (fixedU((fesim%randbedingungen%mpcs(ii)%dependend%nid-1)*ndof+fesim%randbedingungen%mpcs(ii)%dependend%dof) == .TRUE.) then

                        write(*,*) 'fesim%randbedingungen%mpcadds(gg)%sid=',fesim%randbedingungen%mpcadds(gg)%sid
                        write(*,*) '   fesim%randbedingungen%mpcs(ii)%sid=',fesim%randbedingungen%mpcs(ii)%sid
                        write(*,*) 'You tried to set an MPC on a DoF that'
                        write(*,*) 'has also been locked been via an SPC.'
                        err_code=1
                        goto 9999

                    else

                        fixedU((fesim%randbedingungen%mpcs(ii)%dependend%nid-1)*ndof+fesim%randbedingungen%mpcs(ii)%dependend%dof) = .TRUE.

                    end if
                    do jj=1,size(fesim%randbedingungen%mpcs(ii)%independend,1)
                        fesim%randbedingungen%mpcs(ii)%independend(jj)%fac = fesim%randbedingungen%mpcs(ii)%independend(jj)%fac / (-fesim%randbedingungen%mpcs(ii)%dependend%fac)

                        ! Freiheit des Independend Knoten auf einen negativen Wert setzen, falls diese Freiheit bereits durch ein SPC gesperrt ist und später nicht weiter betrachtet werden darf.
                        if (fixedU((fesim%randbedingungen%mpcs(ii)%independend(jj)%nid-1)*ndof+fesim%randbedingungen%mpcs(ii)%independend(jj)%dof) .EQ. .TRUE.) then
                            fesim%randbedingungen%mpcs(ii)%independend(jj)%dof = -fesim%num_dof-1
                        end if
                    end do
                    fesim%randbedingungen%mpcs(ii)%dependend%fac = 1.d0
                end if
                !end do
            end if
        end if
    end do
end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'boundary_mpc_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return

end subroutine boundary_mpc_control
