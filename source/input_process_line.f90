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
!> subroutine for file-handling
!
!> @details
!> Subroutine, for "reading" words from line of words and blanks
!> stores words in array
!> Subroutine can handle multiple words per line seperated by blanks
!> lines can start with blanks
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit 07.12.2009
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 07.11.2012
!
!> $Id: input_process_line.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine input_process_line (fesim,line,kk,array,n)

use fesimulation_typen
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
character(*), intent(in)    :: line     !< line to be processed
integer, intent(in)         :: n        !< values can not be changed in subroutine
!
! Output
!
 character(*), dimension(n) :: array    !< array of found words
!
! Input+Output
!
integer, intent(inout)      :: kk       !< word index
type(fe_simulation)         :: fesim    !< index in the words array
!
! Inner
!
integer                     :: begin_index, end_index,jj,ll
character(50)               :: word     ! current word
character(50)               :: varname
integer                     :: err_code=0
logical                     :: flag
!
! =================================================================================================
!

begin_index = 0                                 ! start index of word
end_index   = 0                                 ! end index of word
jj = 1                                          ! character position in the line
flag = .false.
err_code = 0

! loop over line characters

do while (jj <= len(line))
!do jj=1,len(line)
   
   ! avoid prefaced blanks in the beginning of a line
   
   if ((begin_index == 0) .and. (line(jj:jj) /= ' ')) then
      begin_index = jj
   end if
   
   ! start with first actual character thats not a blank
   
   if (begin_index > 0) then
      
     ! avoid blanks between words in one line and postpositioned blanks
     
     if (line(jj:jj) == ' ') then
       end_index = jj - 1
     end if
      
     if (jj == len(line)) then
       end_index = jj
     end if
      
     ! word processing
      
     if (end_index > 0) then
      
       kk = kk + 1
       word = line(begin_index:end_index)
       array(kk) = word
    
       do ll = 1,len(word)
         if (word(ll:ll) == '=') then
           flag = .true.
           varname = word(1:(ll-1))
           exit
         end if
      end do
    
      if (flag == .false.) then
    
      ! check for Keywords
    
        if (word == 'node') then
            fesim%is_node  = .true.
        else if (word == 'quad8') then
            fesim%is_quad8 = .true.
        else if (word == 'lsolid20') then
            fesim%is_lsolid20 = .true.
        else if (word == 'load') then
            fesim%is_load  = .true.
        else if (word == 'force') then
            fesim%is_force = .true.
        else if (word == 'moment') then
            fesim%is_moment= .true.
        else if (word == 'p8load') then
            fesim%is_p8load = .true.
        else if (word == 'p20load') then
            fesim%is_p20load = .true.
        else if (word == 'aeroload2d') then
            fesim%is_aeroload2d = .true.
        else if (word == 'aeroload3d') then
            fesim%is_aeroload3d = .true.
        else if (word == 'temperature') then
            fesim%is_temperature = .true.
        else if (word == 'beam2temp') then
            fesim%is_beam2temp = .true.
        else if (word == 'quad8temp') then
            fesim%is_quad8temp = .true.
        else if (word == 'lsolid20temp') then
            fesim%is_lsolid20temp = .true.
        else if (word == 'p2load') then
            fesim%is_p2load= .true.
        else if (word == 'mat1') then
            fesim%is_mat1  = .true.
        else if (word == 'mat8') then
            fesim%is_mat8  = .true.
        else if (word == 'mat20') then
            fesim%is_mat20 = .true.
        else if (word == 'pshell') then
            fesim%is_pshell= .true.
        else if (word == 'pcomp') then
            fesim%is_pcomp = .true.
        else if (word == 'plsolid') then
            fesim%is_plsolid= .true.
        else if (word == 'spcadd') then
            fesim%is_spcadd= .true.
        else if (word == 'spc1') then
            fesim%is_spc1  = .true.
        else if (word == 'coord') then
            fesim%is_coord = .true.
        else if (word == 'subcase') then
            fesim%is_subcase = .true.
        else if (word == 'lam8') then
            fesim%is_lam8 = .true.
        else if (word == 'lam20') then
            fesim%is_lam20 = .true.
        else if (word == 'coupling') then
            fesim%is_coupling = .true.
        else if (word == 'beam2') then
            fesim%is_beam2 = .true.
        else if (word == 'pbeam') then
            fesim%is_pbeam = .true.
        else if (word == 'mpc') then
            fesim%is_mpc = .true.
        else if (word == 'mpcadd') then
            fesim%is_mpcadd = .true.
        else if (word == 'contact_node_beam2') then
            fesim%is_contact_node_beam2 = .true.
        else if (word == 'contact_node_quad8') then
            fesim%is_contact_node_quad8 = .true.
        else if (word == 'contact_node_lsolid20') then
            fesim%is_contact_node_lsolid20 = .true.
        else if (word == 'failure') then
            fesim%is_failure = .true.

        !!!!!!!!!!! Thermal: Start
        else if (word == 'thermal') then
            fesim%is_thermal = .true.
        else if (word == 'matthermal') then
            fesim%is_mat_thermal = .true.
        else if (word == 'matthermaltd') then
            fesim%is_mat_thermal_td = .true.
        else if (word == 'tbc') then
            fesim%is_tbc = .true.
        else if (word == 'fluxbc') then
            fesim%is_flux_bc = .true.
        else if (word == 'convectionbc') then
            fesim%is_convection_bc = .true.
        !!!!!!!!!!! Thermal: End
            
        else if (word == 'outputvtk') then
            fesim%ausgabe%outputVTK = .true.
        else if (word == 'outputuser') then
            fesim%ausgabe%outputUser = .true.
        else if (word == 'outputshort') then
            fesim%ausgabe%outputShort = .true.
        else if (word == 'outputkoopt') then
            fesim%ausgabe%outputKoopt = .true.
        else if (word == 'outputadvila') then
            fesim%ausgabe%outputAdviLa = .true.
        else if (word == 'outputoptitube') then
            fesim%ausgabe%outputOptitube = .true.
        else if (word == 'outputglawi') then
            fesim%ausgabe%outputGlawi = .true.
        else if (word == 'outputhymowi') then
            fesim%ausgabe%outputHyMoWi = .true.
        else if (word == 'outputelemcoord') then
            fesim%ausgabe%outputElemCoord = .true.
        else if (word == 'outputboundcond') then
            fesim%ausgabe%outputBoundCond = .true.
        else if (word == 'outputapamepressures') then
            fesim%ausgabe%outputApamePressures = .true.
        else if (word == 'skipfailed') then
            fesim%skipFailed = .true.
        else if (word == 'calculatetse') then
            fesim%calculateTSE = .true.
        else if (word == 'calculateelementaltse') then
            fesim%calculateElementalTSE = .true.
        else if (word == 'calculateReactForce') then
            fesim%calculateReactForce = .true.
        else if (word == 'globalReactForce') then
            fesim%globalReactForce = .true.
        end if
            
      else if (flag == .true.) then
            
        if (varname == 'sol') then
            fesim%sol = CharToInt(word)
        else if (varname == 'numEigVal') then
            fesim%numEigVal = CharToInt(word)
        else if (varname == 'blocksize') then
            fesim%blocksize = CharToInt(word)
        else if (varname == 'shellResPos') then
            fesim%shellResPos = CharToInt(word)
        else if (varname == 'xmin') then
            fesim%ausgabe%xmin = CharToDou(word)
        else if (varname == 'xmax') then
            fesim%ausgabe%xmax = CharToDou(word)
        else if (varname == 'ymin') then
            fesim%ausgabe%ymin = CharToDou(word)
        else if (varname == 'ymax') then
            fesim%ausgabe%ymax = CharToDou(word)
        else if (varname == 'zmin') then
            fesim%ausgabe%zmin = CharToDou(word)
        else if (varname == 'zmax') then
            fesim%ausgabe%zmax = CharToDou(word)
        else if (varname == 'aeroDisplEps') then
            fesim%aeroDisplEps = CharToDou(word)
        end if
            
      else
        write(*,*) 'Error processing flag'
        err_code=1
        goto 9999
      end if
         
      ! re-initialize indices
      begin_index = 0
      end_index   = 0
         
    end if
  end if
   
  ! update index
   
  jj = jj+1
   
end do

!
! =================================================================================================
!
! Error handling
!

9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'input_process_line'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

contains

 function CharToInt(word)
 
   implicit none
   
   character(50)               :: word          ! current word
   integer                     :: CharToInt
   
   do ll=1,len(word)
     if (word(ll:ll) == '=') then
       READ(word((ll+1):len(word)),'(I)') CharToInt
     end if
   end do
 
 end function CharToInt

 function CharToDou(word)
 
   implicit none
   
   character(50)               :: word          ! current word
   double precision            :: CharToDou
   
   do ll=1,len(word)
     if (word(ll:ll) == '=') then
       READ(word((ll+1):len(word)),'(F)') CharToDou
     end if
   end do
 
 end function CharToDou
 
 
end subroutine input_process_line
