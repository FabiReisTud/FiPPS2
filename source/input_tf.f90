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
!> subroutine for handling input files
!
!> @details
!> Subroutine for creating fields from input data stored in input files
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit 07.12.2009
!
!> $Id: input_tf.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine input_tf (fesim)

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
! 
! Output
! 
!
! Input+Output
!
type(fe_simulation)     :: fesim
! 
! Inner
! 
integer, parameter      :: n=255                            ! maximum number of words
character(50)           :: fname, line, control_array(n)    ! filename with maximum length of (XX) characters
integer                 :: io_error, read_error             ! Error-handling parameters
integer                 :: rowcount
integer                 :: nr_lines, nr_words               ! Steuer-Datei Parameter
integer                 :: ii,jj
integer                 :: elcount                          ! Counts overall element number
integer                 :: err_code=0
character(3)            :: angleType
integer                 :: num_master_dofs
character(17)           :: form
character(27)           :: form2
integer                 :: count_tresca, count_mises, count_maxpstress, count_puck, count_hill, count_norris
integer                 :: count_fibre, count_maxstrain, count_cuntze
integer                 :: count_maxstrain3d, count_tsaiwu3d
character(400)          :: line_long
character(10)           :: failureType
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
! open input-entry control-file
fname='control.fipps'
open(unit=19,file=fname,status='old',action='read', iostat=io_error)

! error checking
if (io_error /= 0) then
        write(*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
end if

! read 1st actual value
read (19,*,iostat=read_error) rowcount        ! Zeilenanzahl aus erster Zeile
! error handling
if (read_error /= 0) then
   write(*,*) 'Error reading file ', fname
   err_code=1
   goto 9999
end if

! read control file
nr_lines = 0                                                        ! number of lines in control file
nr_words = 0                                                        ! number of words found

do while (.true.)
        nr_lines=nr_lines+1
        if (nr_lines .GT. rowcount) GOTO 100
        read(unit=19,fmt='(A)',iostat=read_error,end=100) line        ! read file and jump to marker end=XXX when last low reached
        if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                err_code=1
                goto 9999
        end if
        line=trim(line)                                                ! Remove trailing blank characters of string
        ! search for words in the line and store them in array
        call input_process_line(fesim,line,nr_words,control_array,n)
end do

100 continue
close(unit=19,iostat=io_error)
if (io_error /= 0) then
        write(*,*) 'Error closing file', fname
        err_code=1
        goto 9999
end if

! Keyword processing

if ((fesim%is_aeroload2d .eq. .true.) .and. (fesim%is_aeroload3d .eq. .true.)) then
        write(*,*) "Keywords 'aeroload2d' and 'aeroload3d' may not"
        write(*,*) "be used in combination."
        err_code=1
        goto 9999
end if

if ((fesim%is_temperature .eq. .true.) .and. (fesim%is_beam2temp .eq. .true.)) then
        write(*,*) "Keywords 'temperature' and 'beam2temp' may not"
        write(*,*) "be used in combination."
        err_code=1
        goto 9999
end if

if ((fesim%is_temperature .eq. .true.) .and. (fesim%is_quad8temp .eq. .true.)) then
        write(*,*) "Keywords 'temperature' and 'quad8temp' may not"
        write(*,*) "be used in combination."
        err_code=1
        goto 9999
end if

if ((fesim%is_temperature .eq. .true.) .and. (fesim%is_lsolid20temp .eq. .true.)) then
        write(*,*) "Keywords 'temperature' and 'lsolid20temp' may not"
        write(*,*) "be used in combination."
        err_code=1
        goto 9999
end if

if (fesim%is_node == .false.) then
        write(*,*) 'No nodes defined on ', fname
        err_code=1
        goto 9999
else
        fname='nodes.fipps'
        open (unit=20,file=fname,status='old',action='read', iostat=io_error)
                ! unit=x, x=Dateinummer (Ganzzahl,>10,da oft Nummern <10 fix zugeordnet, z.B. der Standardein-, ausgabe)
                ! file=x, x=Dateiname
                ! status, old-Datei existiert bereits. default: 'unknown'
                ! action, read-lesen, write-schreiben, readwrite-lesen+schreiben
        if (io_error == 0) then
                ! read 1st actual value
                read (20,*,iostat=read_error) rowcount        ! Zeilenanzahl aus erster Zeile
                ! error handling
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                ! allocate memory for type and number of rows
                allocate(fesim%knoten%nodes(1:rowcount))
                ! set fixed values (Panelmodeller)
                do ii=1,rowcount
                   fesim%knoten%nodes(ii)%nid = ii
                end do
                ! read values
                do ii=1,rowcount
                        read (20,'(I10,3E24.16)',iostat=read_error) fesim%knoten%nodes(ii)%cid, (fesim%knoten%nodes(ii)%coords(jj), jj=1,3)
!                        read (20,'(I10,I10,3E23.16)',iostat=read_error) fesim%knoten%nodes(ii)%nid, fesim%knoten%nodes(ii)%cid,(fesim%knoten%nodes(ii)%fesim%koordinatensysteme%coords(jj), jj=1,3)
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                ! Close file
                close (20, iostat=io_error)
                fesim%num_nodes = size(fesim%knoten%nodes,1)
                ! Error handling
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_load == .false.) then
        write(*,*) 'No loads defined on ', fname
        err_code=1
        goto 9999
else
        fname='loads.fipps'
        open (unit=21,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (21,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%lasten%loads(1:rowcount))
                ! set fixed values
                do ii=1,rowcount
                   fesim%lasten%loads(ii)%sfac  = 1.D0
                end do
                ! read values
                do ii=1,rowcount
!                        read (21,'(I10,2E23.16,I10)',iostat=read_error) fesim%lasten%loads(ii)%lcid, fesim%lasten%loads(ii)%sfac, fesim%lasten%loads(ii)%sfaci, fesim%lasten%loads(ii)%lidi
                        read (21,'(I10,E23.16,I10)',iostat=read_error) fesim%lasten%loads(ii)%lcid, fesim%lasten%loads(ii)%sfaci, fesim%lasten%loads(ii)%lidi
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (21, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_temperature /= .false.) then
    fname='temperatures.fipps'
    open (unit=22,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (22,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%lasten%nodetemps(1:rowcount))
        do ii=1,rowcount
            read (22,'(2I10,E23.16)',iostat=read_error) fesim%lasten%nodetemps(ii)%lid, fesim%lasten%nodetemps(ii)%nid, fesim%lasten%nodetemps(ii)%temp
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (22, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_beam2temp /= .false.) then
    fname='beam2temps.fipps'
    open (unit=22,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (22,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%lasten%beam2temps(1:rowcount))
        do ii=1,rowcount
            read (22,'(2I10,E23.16)',iostat=read_error) fesim%lasten%beam2temps(ii)%lid, fesim%lasten%beam2temps(ii)%eid, fesim%lasten%beam2temps(ii)%temp
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (22, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_quad8temp /= .false.) then
    fname='quad8temps.fipps'
    open (unit=22,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (22,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%lasten%quad8temps(1:rowcount))
        do ii=1,rowcount
            read (22,'(2I10,E23.16)',iostat=read_error) fesim%lasten%quad8temps(ii)%lid, fesim%lasten%quad8temps(ii)%eid, fesim%lasten%quad8temps(ii)%temp
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (22, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_lsolid20temp /= .false.) then
    fname='lsolid20temps.fipps'
    open (unit=22,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (22,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%lasten%lsolid20temps(1:rowcount))
        do ii=1,rowcount
            read (22,'(2I10,E23.16)',iostat=read_error) fesim%lasten%lsolid20temps(ii)%lid, fesim%lasten%lsolid20temps(ii)%eid, fesim%lasten%lsolid20temps(ii)%temp
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (22, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_p2load /= .false.) then
    fname='p2loads.fipps'
        open (unit=28,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (28,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%lasten%p2loads(1:rowcount))
                do ii=1,rowcount
                        read (28,'(3I10,2E23.16,L1,I10)',iostat=read_error) fesim%lasten%p2loads(ii)%lid,fesim%lasten%p2loads(ii)%eid1,fesim%lasten%p2loads(ii)%dir, &
                                      & (fesim%lasten%p2loads(ii)%pi(jj), jj=1,2),fesim%lasten%p2loads(ii)%thru,fesim%lasten%p2loads(ii)%eid2
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (28, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_spcadd == .false.) then
        write(*,*) 'No spcs defined on ', fname
        err_code=1
        goto 9999
else
        fname='spcadd.fipps'
        open (unit=22,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (22,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%randbedingungen%spcadds(1:rowcount))
                do ii=1,rowcount
                        read (22,'(2I10)',iostat=read_error) fesim%randbedingungen%spcadds(ii)%scid, fesim%randbedingungen%spcadds(ii)%sid
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (22, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if


if (fesim%is_mpcadd == .true.) then
        fname='mpcadd.fipps'
        open (unit=22,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (22,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%randbedingungen%mpcadds(1:rowcount))
                do ii=1,rowcount
                        read (22,'(2I10)',iostat=read_error) fesim%randbedingungen%mpcadds(ii)%scid, fesim%randbedingungen%mpcadds(ii)%sid
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (22, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

! Elements

elcount=0

if (fesim%is_beam2 /= .false.) then
        fname='beam2.fipps'
        open (unit=24,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (24,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%elemente%beam2s(1:rowcount))
                ! set fixed values
                do ii=1,rowcount
                   elcount = elcount + 1
                   fesim%elemente%beam2s(ii)%eid    = elcount
                end do
                ! read values
                do ii=1,rowcount
                        read (24,'(3I10,3E23.16,I10)',iostat=read_error) fesim%elemente%beam2s(ii)%pid,(fesim%elemente%beam2s(ii)%nids(jj),jj=1,2), (fesim%elemente%beam2s(ii)%xi(jj), jj=1,3), fesim%elemente%beam2s(ii)%n0
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=2
                                goto 9999
                        end if
                end do
                close (24, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_quad8 /= .false.) then
        fname='quad8.fipps'
        open (unit=24,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (24,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%elemente%quad8s(1:rowcount))
                ! set fixed values
                do ii=1,rowcount
                   elcount = elcount + 1
                   fesim%elemente%quad8s(ii)%eid    = elcount
                   fesim%elemente%quad8s(ii)%offset = 0.D0
                end do
                ! read values
                do ii=1,rowcount
                        read (24,'(9I10,E23.16)',iostat=read_error) fesim%elemente%quad8s(ii)%pid,(fesim%elemente%quad8s(ii)%nids(jj),jj=1,8), fesim%elemente%quad8s(ii)%theta
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=2
                                goto 9999
                        end if
                end do
                close (24, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_lsolid20 /= .false.) then
  fname='lsolid20.fipps'
  open (unit=24,file=fname,status='old',action='read', iostat=io_error)
  if (io_error == 0) then
    read (24,*,iostat=read_error) rowcount
    if (read_error /= 0) then
      write(*,*) 'Error reading file ', fname
      write(*,*) 'There was a problem reading the number of rows of the file'
      err_code=1
      goto 9999
    end if
    allocate(fesim%elemente%lsolid20s(1:rowcount))
    ! set fixed values
    do ii=1,rowcount
      elcount = elcount + 1
      fesim%elemente%lsolid20s(ii)%eid    = elcount
    end do
    ! read values
    do ii=1,rowcount
      read (24,'(21I10)',iostat=read_error) fesim%elemente%lsolid20s(ii)%pid&
          &,(fesim%elemente%lsolid20s(ii)%nids(jj),jj=1,20)
      if (read_error /= 0) then
        write(*,*) 'Error reading file ', fname
        write(*,*) 'There was an error processing line', ii
        err_code=2
        goto 9999
      end if
    end do
    close (24, iostat=io_error)
    if (io_error /= 0) then
      write(*,*) 'Error closing file ', fname
      err_code=1
      goto 9999
    end if
  else
    write (*,*) 'Error opening file ', fname
    err_code=1
    goto 9999
  end if
end if

if (fesim%is_force /= .false.) then
        fname='forces.fipps'
        open (unit=25,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (25,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%lasten%forces(1:rowcount))
                do ii=1,rowcount
                        read (25,'(2I10,4E23.16)',iostat=read_error) fesim%lasten%forces(ii)%lid, fesim%lasten%forces(ii)%nid, fesim%lasten%forces(ii)%fac, (fesim%lasten%forces(ii)%ni(jj), jj=1,3)
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                write(*,*) 'iostat ', read_error
                                err_code=1
                                goto 9999
                        end if
                end do
                close (25, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_moment /= .false.) then
        fname='moments.fipps'
        open (unit=26,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (26,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%lasten%moments(1:rowcount))
                do ii=1,rowcount
                        read (26,'(2I10,4E23.16)',iostat=read_error) fesim%lasten%moments(ii)%lid,fesim%lasten%moments(ii)%nid,fesim%lasten%moments(ii)%fac,(fesim%lasten%moments(ii)%ni(jj), jj=1,3)
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (26, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_p8load /= .false.) then
        fname='p8loads.fipps'
        open (unit=28,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (28,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%lasten%p8loads(1:rowcount))
                do ii=1,rowcount
                        read (28,'(3I10,4E23.16,L1,I10)',iostat=read_error) fesim%lasten%p8loads(ii)%lid,fesim%lasten%p8loads(ii)%eid1,fesim%lasten%p8loads(ii)%cid,(fesim%lasten%p8loads(ii)%pi(jj), jj=1,4), &
                                                                          & fesim%lasten%p8loads(ii)%thru,fesim%lasten%p8loads(ii)%eid2
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (28, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_p20load /= .false.) then
    fname='p20loads.fipps'
    open (unit=28,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (28,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%lasten%p20loads(1:rowcount))
        do ii=1,rowcount
            read (28,'(2I10,E23.16,I10,L1,I10)',iostat=read_error) fesim%lasten%p20loads(ii)%lid,fesim%lasten%p20loads(ii)%eid1,fesim%lasten%p20loads(ii)%p,fesim%lasten%p20loads(ii)%surf, &
                                      & fesim%lasten%p20loads(ii)%thru,fesim%lasten%p20loads(ii)%eid2
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (28, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_aeroload2d /= .false.) then
    fname='aeroload2ds.fipps'
    open (unit=28,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (28,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%lasten%aeroload2ds(1:rowcount))
        do ii=1,rowcount
            read (28,'(2I10,E23.16)',iostat=read_error) fesim%lasten%aeroload2ds(ii)%lid,fesim%lasten%aeroload2ds(ii)%mthd,fesim%lasten%aeroload2ds(ii)%dfac
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
            if (fesim%lasten%aeroload2ds(ii)%dfac .LE. 0.d0 .OR. fesim%lasten%aeroload2ds(ii)%dfac .GT. 1.d0) then
                fesim%lasten%aeroload2ds(ii)%dfac = 1.d0
            end if
        end do
        close (28, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_aeroload3d /= .false.) then
    fname='aeroload3ds.fipps'
    open (unit=28,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (28,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%lasten%aeroload3ds(1:rowcount))
        do ii=1,rowcount
            read (28,'(2I10,E23.16)',iostat=read_error) fesim%lasten%aeroload3ds(ii)%lid,fesim%lasten%aeroload3ds(ii)%mthd,fesim%lasten%aeroload3ds(ii)%dfac
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
            if (fesim%lasten%aeroload3ds(ii)%dfac .LE. 0.d0 .OR. fesim%lasten%aeroload3ds(ii)%dfac .GT. 1.d0) then
                fesim%lasten%aeroload3ds(ii)%dfac = 1.d0
            end if
        end do
        close (28, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_mat1 /= .false.) then
        fname='mat1.fipps'
        open (unit=29,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (29,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%materialien%mat1s(1:rowcount))
                do ii=1,rowcount
                        read (29,'(I10,7E23.16,4I10)',iostat=read_error) fesim%materialien%mat1s(ii)%mid, fesim%materialien%mat1s(ii)%ym, fesim%materialien%mat1s(ii)%sm, fesim%materialien%mat1s(ii)%nu, fesim%materialien%mat1s(ii)%rho, fesim%materialien%mat1s(ii)%ath,&
                                                                       & fesim%materialien%mat1s(ii)%tref, fesim%materialien%mat1s(ii)%ge, fesim%materialien%mat1s(ii)%fid(1), fesim%materialien%mat1s(ii)%fid(2), fesim%materialien%mat1s(ii)%fid(3), fesim%materialien%mat1s(ii)%fid(4)
                        do jj = 1,4
                          if (fesim%materialien%mat1s(ii)%fid(jj) .LT. 1) then
                            fesim%materialien%mat1s(ii)%fid(jj) = 0
                          end if
                        end do
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (29, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_mat8 /= .false.) then
        fname='mat8.fipps'
        open (unit=30,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (30,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%materialien%mat8s(1:rowcount))
                do ii=1,rowcount
                        read (30,'(I10,11E23.16,4I10)',iostat=read_error)  fesim%materialien%mat8s(ii)%mid,fesim%materialien%mat8s(ii)%ym11,fesim%materialien%mat8s(ii)%ym22,&
                        & fesim%materialien%mat8s(ii)%nu12,fesim%materialien%mat8s(ii)%sm12,fesim%materialien%mat8s(ii)%sm13,&
                        & fesim%materialien%mat8s(ii)%sm23,fesim%materialien%mat8s(ii)%rho,fesim%materialien%mat8s(ii)%ath11,&
                        & fesim%materialien%mat8s(ii)%ath22,fesim%materialien%mat8s(ii)%tref,fesim%materialien%mat8s(ii)%ge,&
                        & fesim%materialien%mat8s(ii)%fid(1), fesim%materialien%mat8s(ii)%fid(2), fesim%materialien%mat8s(ii)%fid(3),&
                        & fesim%materialien%mat8s(ii)%fid(4)
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                        do jj = 1,4
                          if (fesim%materialien%mat8s(ii)%fid(jj) .LT. 1) then
                            fesim%materialien%mat8s(ii)%fid(jj) = 0
                          end if
                        end do
                end do
                close (30, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_mat20 /= .false.) then
  fname='mat20.fipps'
  open (unit=30,file=fname,status='old',action='read', iostat=io_error)
  if (io_error == 0) then
    read (30,*,iostat=read_error) rowcount
    if (read_error /= 0) then
      write(*,*) 'Error reading file ', fname
      write(*,*) 'There was a problem reading the number of rows of the file'
      err_code=1
      goto 9999
    end if
    allocate(fesim%materialien%mat20s(1:rowcount))
    do ii=1,rowcount
      read (30,'(I10,13E23.16,4I10)',iostat=read_error) fesim%materialien%mat20s(ii)%mid,fesim%materialien%mat20s(ii)&
          &%ym11, fesim%materialien%mat20s(ii)%ym22,  fesim%materialien%mat20s(ii)%ym33,  fesim%materialien%mat20s(ii)%nu12,  fesim%materialien%mat20s(ii)&
          &%nu13, fesim%materialien%mat20s(ii)%nu23,  fesim%materialien%mat20s(ii)%sm12,  fesim%materialien%mat20s(ii)%sm13,  fesim%materialien%mat20s(ii)&
          &%sm23, fesim%materialien%mat20s(ii)%ath11, fesim%materialien%mat20s(ii)%ath22, fesim%materialien%mat20s(ii)%ath33, fesim%materialien%mat20s(ii)%rho, &
          &fesim%materialien%mat20s(ii)%fid(1), fesim%materialien%mat20s(ii)%fid(2), fesim%materialien%mat20s(ii)%fid(3), fesim%materialien%mat20s(ii)%fid(4)
      if (read_error /= 0) then
        write(*,*) 'Error reading file ', fname
        write(*,*) 'There was an error processing line', ii
        err_code=1
        goto 9999
      end if
    end do
    close (30, iostat=io_error)
    if (io_error /= 0) then
      write(*,*) 'Error closing file ', fname
      err_code=1
      goto 9999
    end if
  else
    write (*,*) 'Error opening file ', fname
    err_code=1
    goto 9999
  end if
end if

if (fesim%is_failure /= .false.) then
        fname='failure.fipps'
        open (unit=38,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (38,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                count_tresca = 0
                count_mises = 0
                count_maxpstress = 0
                count_puck = 0
                count_hill = 0
                count_norris = 0
                count_fibre = 0
                count_maxstrain = 0
                count_cuntze = 0
                count_maxstrain3d = 0
                count_tsaiwu3d = 0
                do ii=1,rowcount
                    read (38,'(10X,A10)',iostat=read_error) failureType
                    if (TRIM(ADJUSTL(failureType)) == "tresca") then
                        count_tresca = count_tresca + 1
                    else if (TRIM(ADJUSTL(failureType)) == "mises") then
                        count_mises = count_mises + 1
                    else if (TRIM(ADJUSTL(failureType)) == "maxpstress") then
                        count_maxpstress = count_maxpstress + 1
                    else if (TRIM(ADJUSTL(failureType)) == "puck") then
                        count_puck = count_puck + 1
                    else if (TRIM(ADJUSTL(failureType)) == "hill") then
                        count_hill = count_hill + 1
                    else if (TRIM(ADJUSTL(failureType)) == "norris") then
                        count_norris = count_norris + 1
                    else if (TRIM(ADJUSTL(failureType)) == "fibre") then
                        count_fibre = count_fibre + 1
                    else if (TRIM(ADJUSTL(failureType)) == "maxstrain") then
                        count_maxstrain = count_maxstrain + 1
                    else if (TRIM(ADJUSTL(failureType)) == "cuntze") then
                        count_cuntze = count_cuntze + 1
                    else if (TRIM(ADJUSTL(failureType)) == "maxstrain3") then
                        count_maxstrain3d = count_maxstrain3d + 1
                    else if (TRIM(ADJUSTL(failureType)) == "tsaiwu3d") then
                        count_tsaiwu3d = count_tsaiwu3d + 1
                    else
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'unknown failure criteria name'
                        err_code=2
                        goto 9999
                    end if
                end do
                
                if (count_tresca .gt. 0) allocate(fesim%versagenskriterien%failTrescas(1:count_tresca))
                if (count_mises .gt. 0) allocate(fesim%versagenskriterien%failMises(1:count_mises))
                if (count_maxpstress .gt. 0) allocate(fesim%versagenskriterien%failMaxprincstresses(1:count_maxpstress))
                if (count_puck .gt. 0) allocate(fesim%versagenskriterien%failPucks(1:count_puck))
                if (count_hill .gt. 0) allocate(fesim%versagenskriterien%failHills(1:count_hill))
                if (count_norris .gt. 0) allocate(fesim%versagenskriterien%failNorris(1:count_norris))
                if (count_fibre .gt. 0) allocate(fesim%versagenskriterien%failFibres(1:count_fibre))
                if (count_maxstrain .gt. 0) allocate(fesim%versagenskriterien%failMaxStrains(1:count_maxstrain))
                if (count_cuntze .gt. 0) allocate(fesim%versagenskriterien%failCuntzes(1:count_cuntze))
                if (count_maxstrain3d .gt. 0) allocate(fesim%versagenskriterien%failMaxStrain3Ds(1:count_maxstrain3d))
                if (count_tsaiwu3d .gt. 0) allocate(fesim%versagenskriterien%failTsaiWu3Ds(1:count_tsaiwu3d))
                rewind(38)
                read (38,*) 
                count_tresca = 0
                count_mises = 0
                count_maxpstress = 0
                count_puck = 0
                count_hill = 0
                count_norris = 0
                count_fibre = 0
                count_maxstrain = 0
                count_cuntze = 0
                count_maxstrain3d = 0
                count_tsaiwu3d = 0
                do ii=1,rowcount
                    read (38,'(A400)',iostat=read_error) line_long
                    if (read_error /= 0) then
                      write(*,*) 'Error reading file ', fname
                      write(*,*) 'There was an error processing line', ii
                      err_code=1
                      goto 9999
                    end if
                    read (line_long,'(10X,A10)',iostat=read_error) failureType
                    if (TRIM(ADJUSTL(failureType)) == "tresca") then
                    
                        count_tresca = count_tresca + 1
                        read (line_long,'(I10,10x,1E23.16)',iostat=read_error)  &
                          & fesim%versagenskriterien%failTrescas(count_tresca)%fid,fesim%versagenskriterien%failTrescas(count_tresca)%ys
                          
                    else if (TRIM(ADJUSTL(failureType)) == "mises") then
                    
                        count_mises = count_mises + 1
                        read (line_long,'(I10,10x,1E23.16)',iostat=read_error)  &
                          & fesim%versagenskriterien%failMises(count_mises)%fid,fesim%versagenskriterien%failMises(count_mises)%ys
                          
                    else if (TRIM(ADJUSTL(failureType)) == "maxpstress") then
                    
                        count_maxpstress = count_maxpstress + 1
                        read (line_long,'(I10,10x,2E23.16)',iostat=read_error)  &
                          & fesim%versagenskriterien%failMaxprincstresses(count_maxpstress)%fid,fesim%versagenskriterien%failMaxprincstresses(count_maxpstress)%ys,fesim%versagenskriterien%failMaxprincstresses(count_maxpstress)%ysC
                    
                    else if (TRIM(ADJUSTL(failureType)) == "puck") then
                    
                        count_puck = count_puck + 1
                        read (line_long,'(I10,10x,9E23.16)',iostat=read_error)  &
                          & fesim%versagenskriterien%failPucks(count_puck)%fid,fesim%versagenskriterien%failPucks(count_puck)%RParTen,fesim%versagenskriterien%failPucks(count_puck)%RParCom, &
                          & fesim%versagenskriterien%failPucks(count_puck)%RNorTen,fesim%versagenskriterien%failPucks(count_puck)%RNorCom,fesim%versagenskriterien%failPucks(count_puck)%RShear, &
                          & fesim%versagenskriterien%failPucks(count_puck)%Pspd,fesim%versagenskriterien%failPucks(count_puck)%Pspz,fesim%versagenskriterien%failPucks(count_puck)%a0, &
                          & fesim%versagenskriterien%failPucks(count_puck)%lambdamin
                          
                    else if (TRIM(ADJUSTL(failureType)) == "hill") then
                    
                        count_hill = count_hill + 1
                        read (line_long,'(I10,10x,6E23.16)',iostat=read_error)  &
                          & fesim%versagenskriterien%failHills(count_hill)%fid,fesim%versagenskriterien%failHills(count_hill)%RParTen,fesim%versagenskriterien%failHills(count_hill)%RParCom, &
                          & fesim%versagenskriterien%failHills(count_hill)%RNorTen,fesim%versagenskriterien%failHills(count_hill)%RNorCom,fesim%versagenskriterien%failHills(count_hill)%RShear, &
                          & fesim%versagenskriterien%failHills(count_hill)%F12star
                          
                    else if (TRIM(ADJUSTL(failureType)) == "norris") then
                    
                        count_norris = count_norris + 1
                        read (line_long,'(I10,10x,3E23.16)',iostat=read_error)  &
                          & fesim%versagenskriterien%failNorris(count_norris)%fid,fesim%versagenskriterien%failNorris(count_norris)%RPar, &
                          & fesim%versagenskriterien%failNorris(count_norris)%RNor, &
                          & fesim%versagenskriterien%failNorris(count_norris)%RShear
                    
                    else if (TRIM(ADJUSTL(failureType)) == "fibre") then
                    
                        count_fibre = count_fibre + 1
                        read (line_long,'(I10,10x,2E23.16)',iostat=read_error)  &
                          & fesim%versagenskriterien%failFibres(count_fibre)%fid,fesim%versagenskriterien%failFibres(count_fibre)%RParTen,fesim%versagenskriterien%failFibres(count_fibre)%RParCom
                    
                    else if (TRIM(ADJUSTL(failureType)) == "maxstrain") then

                        count_maxstrain = count_maxstrain + 1
                        read (line_long,'(I10,10x,5E23.16,L1)',iostat=read_error)  &
                          & fesim%versagenskriterien%failMaxStrains(count_maxstrain)%fid,fesim%versagenskriterien%failMaxStrains(count_maxstrain)%epsParTen, &
                          & fesim%versagenskriterien%failMaxStrains(count_maxstrain)%epsParCom, fesim%versagenskriterien%failMaxStrains(count_maxstrain)%epsNorTen, &
                          & fesim%versagenskriterien%failMaxStrains(count_maxstrain)%epsNorCom,fesim%versagenskriterien%failMaxStrains(count_maxstrain)%epsShear, &
                          & fesim%versagenskriterien%failMaxStrains(count_maxstrain)%useGlobal
                    
                    else if (TRIM(ADJUSTL(failureType)) == "cuntze") then
                    
                        count_cuntze = count_cuntze + 1
                        read (line_long,'(I10,10x,9E23.16)',iostat=read_error)  &
                        & fesim%versagenskriterien%failCuntzes(count_cuntze)%fid,fesim%versagenskriterien%failCuntzes(count_cuntze)%RParTen,fesim%versagenskriterien%failCuntzes(count_cuntze)%RParCom, &
                        & fesim%versagenskriterien%failCuntzes(count_cuntze)%RNorTen,fesim%versagenskriterien%failCuntzes(count_cuntze)%RNorCom,fesim%versagenskriterien%failCuntzes(count_cuntze)%RShear, &
                        & fesim%versagenskriterien%failCuntzes(count_cuntze)%muNorPar,fesim%versagenskriterien%failCuntzes(count_cuntze)%m
                    
                    else if (TRIM(ADJUSTL(failureType)) == "maxstrain3") then

                        count_maxstrain3d = count_maxstrain3d + 1
                        read(line_long,'(I10,10x,9E23.16,L1)',iostat=read_error)   &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%fid, &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps11Ten,fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps11Com, &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps22Ten,fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps22Com, &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps33Ten,fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps33Com, &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps12Shear, &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps13Shear, &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%eps23Shear, &
                        & fesim%versagenskriterien%failMaxStrain3Ds(count_maxstrain3d)%useGlobal
                    
                    else if (TRIM(ADJUSTL(failureType)) == "tsaiwu3d") then
                    
                        count_tsaiwu3d = count_tsaiwu3d + 1
                        read(line_long,'(I10,10x,12E23.16)',iostat=read_error)   &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%fid, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R11Ten,fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R11Com, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R22Ten,fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R22Com, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R33Ten,fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R33Com, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R12Shear, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R13Shear, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%R23Shear, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%coupl12, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%coupl13, &
                        & fesim%versagenskriterien%failTsaiWu3Ds(count_tsaiwu3d)%coupl23

                    else
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'unknown failre criteria name in line', ii, '  ', TRIM(ADJUSTL(failureType))
                        err_code=2
                        goto 9999
                    end if
                end do
                close (38, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_pshell /= .false.) then
        fname='pshell.fipps'
        open (unit=31,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (31,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%eigenschaften%pshells(1:rowcount))
                do ii=1,rowcount
                        read (31,'(2I10,E23.16,I10,E23.16,I10,4E23.16,I10)',iostat=read_error) fesim%eigenschaften%pshells(ii)%pid, fesim%eigenschaften%pshells(ii)%mid1, fesim%eigenschaften%pshells(ii)%mt, &
                                        & fesim%eigenschaften%pshells(ii)%mid2, fesim%eigenschaften%pshells(ii)%bmr, fesim%eigenschaften%pshells(ii)%mid3, fesim%eigenschaften%pshells(ii)%tst, fesim%eigenschaften%pshells(ii)%nsm, fesim%eigenschaften%pshells(ii)%z1, fesim%eigenschaften%pshells(ii)%z2, &
                                        & fesim%eigenschaften%pshells(ii)%mid4
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                write(*,*) 'Errorcode: ', read_error
                                err_code=1
                                goto 9999
                        end if
                end do
                close (31, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_pcomp /= .false.) then
        fname='pcomp.fipps'
        open (unit=32,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (32,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%eigenschaften%pcomps(1:rowcount))
                ! set fixed values
                do ii=1,rowcount
                   fesim%eigenschaften%pcomps(ii)%lay = 0                                ! Value is calculated from lam8-cards
                end do
                ! read values
                do ii=1,rowcount
                        read (32,'(2I10,3E23.16,A2)',iostat=read_error) fesim%eigenschaften%pcomps(ii)%pid,fesim%eigenschaften%pcomps(ii)%lamid,fesim%eigenschaften%pcomps(ii)%offset,&
                                                                                & fesim%eigenschaften%pcomps(ii)%nsm,fesim%eigenschaften%pcomps(ii)%sb,fesim%eigenschaften%pcomps(ii)%ft
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (32, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_pbeam /= .false.) then
        fname='pbeam.fipps'
        open (unit=32,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (32,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%eigenschaften%pbeams(1:rowcount))
                do ii=1,rowcount
                        read (32,'(2I10,8E23.16,A3,E23.16)',iostat=read_error) fesim%eigenschaften%pbeams(ii)%pid,fesim%eigenschaften%pbeams(ii)%mid,fesim%eigenschaften%pbeams(ii)%AA,fesim%eigenschaften%pbeams(ii)%I11,fesim%eigenschaften%pbeams(ii)%I22,&
                                                                             & fesim%eigenschaften%pbeams(ii)%I12,fesim%eigenschaften%pbeams(ii)%It,fesim%eigenschaften%pbeams(ii)%t1,fesim%eigenschaften%pbeams(ii)%t2,fesim%eigenschaften%pbeams(ii)%angle,angleType,fesim%eigenschaften%pbeams(ii)%nsm
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                        if (angleType == 'deg') then
                          fesim%eigenschaften%pbeams(ii)%angle=fesim%eigenschaften%pbeams(ii)%angle/180.D0*acos(-1.d0)
                        else if (angleType .NE. 'rad') then
                             write(*,*) 'wrong input on value atype on pbeam',ii
                             err_code=1
                             goto 9999
                        end if
                end do
                close (32, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_plsolid /= .false.) then
  fname='plsolid.fipps'
  open (unit=32,file=fname,status='old',action='read', iostat=io_error)
  if (io_error == 0) then
    read (32,*,iostat=read_error) rowcount
    if (read_error /= 0) then
      write(*,*) 'Error reading file ', fname
      write(*,*) 'There was a problem reading the number of rows of the file'
      err_code=1
      goto 9999
    end if
    allocate(fesim%eigenschaften%plsolids(1:rowcount))
    ! set fixed values
    do ii=1,rowcount
      fesim%eigenschaften%plsolids(ii)%lay = 0
      fesim%eigenschaften%plsolids(ii)%cid = 0
      fesim%eigenschaften%plsolids(ii)%globOut = .FALSE.
      fesim%eigenschaften%plsolids(ii)%resLay  = 0
    end do
    ! read values
    do ii=1,rowcount
      read (32,'(4I10,L1)',iostat=read_error)&
          & fesim%eigenschaften%plsolids(ii)%pid, fesim%eigenschaften%plsolids(ii)%lamid, fesim%eigenschaften%plsolids(ii)%cid, fesim%eigenschaften%plsolids(ii)%resLay, fesim%eigenschaften%plsolids(ii)%globOut
      if (read_error /= 0) then
        write(*,*) 'Error reading file ', fname
        write(*,*) 'There was an error processing line', ii
        err_code=1
        goto 9999
      end if
    end do
    close (32, iostat=io_error)
    if (io_error /= 0) then
      write(*,*) 'Error closing file ', fname
      err_code=1
      goto 9999
    end if
  else
    write (*,*) 'Error opening file ', fname
    err_code=1
    goto 9999
  end if
end if

if (fesim%is_spc1 /= .false.) then
        fname='spc1.fipps'
        open (unit=33,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (33,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%randbedingungen%spc1s(1:rowcount))
                do ii=1,rowcount
                        read (33,'(3I10,L1,I10)',iostat=read_error) fesim%randbedingungen%spc1s(ii)%sid, fesim%randbedingungen%spc1s(ii)%dof, fesim%randbedingungen%spc1s(ii)%n1, fesim%randbedingungen%spc1s(ii)%thru, fesim%randbedingungen%spc1s(ii)%nn
!                        read (33,'(3I10,L1,2I10)',iostat=read_error) fesim%randbedingungen%spc1s(ii)%sid,fesim%randbedingungen%spc1s(ii)%dof,fesim%randbedingungen%spc1s(ii)%n1,fesim%randbedingungen%spc1s(ii)%thru,fesim%randbedingungen%spc1s(ii)%nn,fesim%randbedingungen%spc1s(ii)%inc
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (33, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_mpc /= .false.) then
        fname='mpc.fipps'
        open (unit=33,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (33,*,iostat=read_error) rowcount ! Soll zunächst die Anzahl an MPCs sein und NICHT die Anzahl der Zeilen
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%randbedingungen%mpcs(1:rowcount))
                do ii=1,rowcount
                        fesim%randbedingungen%mpcs(ii)%sid = ii
                        read (33,'(3I10,E23.16)',iostat=read_error) num_master_dofs, fesim%randbedingungen%mpcs(ii)%dependend%nid, fesim%randbedingungen%mpcs(ii)%dependend%dof, fesim%randbedingungen%mpcs(ii)%dependend%fac
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                        fesim%randbedingungen%mpcs(ii)%mpc_type = 0 ! Set as external MPC
                        allocate(fesim%randbedingungen%mpcs(ii)%independend(1:num_master_dofs))
                        write(form,'(A1,I2.2,A14)') '(', num_master_dofs, '(2I10,E23.16))'
                        read (33,form,iostat=read_error) (fesim%randbedingungen%mpcs(ii)%independend(jj)%nid, fesim%randbedingungen%mpcs(ii)%independend(jj)%dof, fesim%randbedingungen%mpcs(ii)%independend(jj)%fac,jj=1,num_master_dofs)
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (33, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_coord /= .false.) then
        fname='coord.fipps'
        open (unit=34,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (34,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%koordinatensysteme%coords(1:rowcount))
                do ii=1,rowcount
                        read (34,'(I10,9E23.16)',iostat=read_error) fesim%koordinatensysteme%coords(ii)%cid, (fesim%koordinatensysteme%coords(ii)%transMat(jj,1), jj=1,3), &
                                                                                                             (fesim%koordinatensysteme%coords(ii)%transMat(jj,2), jj=1,3), &
                                                                                                             (fesim%koordinatensysteme%coords(ii)%transMat(jj,3), jj=1,3)
                        fesim%koordinatensysteme%coords(ii)%transMat(1:3,1) = fesim%koordinatensysteme%coords(ii)%transMat(1:3,1)/(sqrt(dot_product(fesim%koordinatensysteme%coords(ii)%transMat(1:3,1),fesim%koordinatensysteme%coords(ii)%transMat(1:3,1))))
                        fesim%koordinatensysteme%coords(ii)%transMat(1:3,2) = fesim%koordinatensysteme%coords(ii)%transMat(1:3,2)/(sqrt(dot_product(fesim%koordinatensysteme%coords(ii)%transMat(1:3,2),fesim%koordinatensysteme%coords(ii)%transMat(1:3,2))))
                        fesim%koordinatensysteme%coords(ii)%transMat(1:3,3) = fesim%koordinatensysteme%coords(ii)%transMat(1:3,3)/(sqrt(dot_product(fesim%koordinatensysteme%coords(ii)%transMat(1:3,3),fesim%koordinatensysteme%coords(ii)%transMat(1:3,3))))
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (34, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_subcase /= .false.) then
        fname='subcase.fipps'
        open (unit=35,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (35,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%lasten%subcases(1:rowcount))
                do ii=1,rowcount
                
                        read (35,'(4I10,5L1)',iostat=read_error) fesim%lasten%subcases(ii)%scid, fesim%lasten%subcases(ii)%spcaddid, fesim%lasten%subcases(ii)%loadid, fesim%lasten%subcases(ii)%mpcaddid, &
                                                                   & fesim%lasten%subcases(ii)%skipBuckling, fesim%lasten%subcases(ii)%upgeom, fesim%lasten%subcases(ii)%upstress, &
                                                                   & fesim%lasten%subcases(ii)%output, fesim%lasten%subcases(ii)%readApameInput
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                        if (fesim%lasten%subcases(ii)%upgeom .eq. .true. .or. fesim%lasten%subcases(ii)%upstress .eq. .true.) then
                          fesim%is_multistep = .true.
                        end if
                end do
                close (35, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_lam8 /= .false.) then
        fname='lam8.fipps'
        open (unit=36,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (36,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%laminate%lam8s(1:rowcount))
                do ii=1,rowcount
                        read (36,'(3I10,2E23.16,A3)',iostat=read_error) fesim%laminate%lam8s(ii)%lamid,fesim%laminate%lam8s(ii)%plyid,fesim%laminate%lam8s(ii)%mat8id,fesim%laminate%lam8s(ii)%th,&
                                                                        & fesim%laminate%lam8s(ii)%angle,angleType
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                        if (angleType == 'deg') then
                          fesim%laminate%lam8s(ii)%angle=fesim%laminate%lam8s(ii)%angle/180.D0*acos(-1.d0)
                        else if (angleType .NE. 'rad') then
                             write(*,*) 'wrong input on value atype on lam8 ',ii
                             err_code=1
                             goto 9999
                        end if
                end do
                close (36, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_lam20 /= .false.) then
  fname='lam20.fipps'
  open (unit=36,file=fname,status='old',action='read', iostat=io_error)
  if (io_error == 0) then
    read (36,*,iostat=read_error) rowcount
    if (read_error /= 0) then
      write(*,*) 'Error reading file ', fname
      write(*,*) 'There was a problem reading the number of rows of the file'
      err_code=1
      goto 9999
    end if
    allocate(fesim%laminate%lam20s(1:rowcount))
    do ii=1,rowcount
      if (fesim%is_multistep /= .false.) then
        allocate (fesim%laminate%lam20s(ii)%mat20id(1:size(fesim%lasten%subcases,1)))
        write(form2,'(A21,I2.2,A4)') '(3I10,2E23.16,A3,I10,', size(fesim%lasten%subcases,1)-1, 'I10)'
        read (36,form2,iostat=read_error) fesim%laminate%lam20s(ii)%lamid,fesim%laminate%lam20s(ii)%plyid,fesim%laminate%lam20s(ii)%mat20id(1),fesim%laminate%lam20s(ii)%th,&
                                        & fesim%laminate%lam20s(ii)%angle,angleType,fesim%laminate%lam20s(ii)%nop,(fesim%laminate%lam20s(ii)%mat20id(jj),jj=2,size(fesim%lasten%subcases,1))
      else
        allocate (fesim%laminate%lam20s(ii)%mat20id(1))
        read (36,'(3I10,2E23.16,A3,I10)',iostat=read_error) fesim%laminate%lam20s(ii)%lamid,fesim%laminate%lam20s(ii)%plyid,fesim%laminate%lam20s(ii)%mat20id(1),fesim%laminate%lam20s(ii)%th,&
                                                          & fesim%laminate%lam20s(ii)%angle,angleType,fesim%laminate%lam20s(ii)%nop
      end if
      if (read_error /= 0) then
          write(*,*) 'Error reading file ', fname
          write(*,*) 'There was an error processing line', ii
          err_code=1
          goto 9999
      end if
      if (angleType == 'deg') then
                    fesim%laminate%lam20s(ii)%angle=fesim%laminate%lam20s(ii)%angle/180.D0*acos(-1.d0)
      else if (angleType .NE. 'rad') then
        write(*,*) 'wrong input on value atype on lam20 ',ii
        err_code=1
        goto 9999
      end if
    end do
    close (36, iostat=io_error)
    if (io_error /= 0) then
        write(*,*) 'Error closing file ', fname
        err_code=1
        goto 9999
    end if
  else
    write (*,*) 'Error opening file ', fname
    err_code=1
    goto 9999
  end if
end if



if (fesim%is_coupling /= .false.) then
        fname='couplings.fipps'
        open (unit=37,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (37,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                allocate(fesim%randbedingungen%couplings(1:rowcount))
                do ii=1,rowcount
                        read (37,'(3I10)',iostat=read_error) fesim%randbedingungen%couplings(ii)%cpsid,fesim%randbedingungen%couplings(ii)%dof,fesim%randbedingungen%couplings(ii)%nid
                        if (read_error /= 0) then
                                write(*,*) 'Error reading file ', fname
                                write(*,*) 'There was an error processing line', ii
                                err_code=1
                                goto 9999
                        end if
                end do
                close (37, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_contact_node_beam2 /= .false.) then
        fname='contact_node_beam2.fipps'
        open (unit=37,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (37,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                call read_contact_beam2_node(37, rowcount)
                close (37, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_contact_node_quad8 /= .false.) then
        fname='contact_node_quad8.fipps'
        open (unit=37,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (37,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                call read_contact_quad8_node(37, rowcount)
                close (37, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_contact_node_lsolid20 /= .false.) then
        fname='contact_node_lsolid20.fipps'
        open (unit=37,file=fname,status='old',action='read', iostat=io_error)
        if (io_error == 0) then
                read (37,*,iostat=read_error) rowcount
                if (read_error /= 0) then
                        write(*,*) 'Error reading file ', fname
                        write(*,*) 'There was a problem reading the number of rows of the file'
                        err_code=1
                        goto 9999
                end if
                call read_contact_lsolid20_node(37, rowcount)
                close (37, iostat=io_error)
                if (io_error /= 0) then
                        write(*,*) 'Error closing file ', fname
                        err_code=1
                        goto 9999
                end if
        else
                write (*,*) 'Error opening file ', fname
                err_code=1
                goto 9999
        end if
end if

if (fesim%is_aeroload2d /= .false.) then
    fname='aeroelem2structnode2d.fipps'
    open (unit=37,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        ! Einlesen der Knoteninformationen
        read (37,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%internals%aeroElem2structNode(1:rowcount))
        do ii=1,rowcount
            read (37,'(2I10,E23.16)',iostat=read_error) &
                & fesim%internals%aeroElem2structNode(ii)%nodeID, &
                & fesim%internals%aeroElem2structNode(ii)%elemID, &
                & fesim%internals%aeroElem2structNode(ii)%xi
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        ! Einlesen der betroffenen Elemente
        read (37,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%internals%structBeam2IDs(1:rowcount))
        do ii=1,rowcount
            read (37,'(I10)',iostat=read_error) fesim%internals%structBeam2IDs(ii)
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (37, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if


    fname='structelem2aeronode2d.fipps'
    open (unit=37,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (37,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%internals%structElem2aeroNode(1:rowcount))
        do ii=1,rowcount
            read (37,'(2I10,E23.16)',iostat=read_error) &
                & fesim%internals%structElem2aeroNode(ii)%nodeID, &
                & fesim%internals%structElem2aeroNode(ii)%elemID, &
                & fesim%internals%structElem2aeroNode(ii)%xi
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (37, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_aeroload3d /= .false.) then
    fname='aeroelem2structnode3d.fipps'
    open (unit=37,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        ! Einlesen der Knoteninformationen
        read (37,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%internals%aeroElem2structNode(1:rowcount))
        do ii=1,rowcount
            read (37,'(2I10,2E23.16)',iostat=read_error) &
                & fesim%internals%aeroElem2structNode(ii)%nodeID, &
                & fesim%internals%aeroElem2structNode(ii)%elemID, &
                & fesim%internals%aeroElem2structNode(ii)%xi,     &
                & fesim%internals%aeroElem2structNode(ii)%eta
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        ! Einlesen der betroffenen Elemente
        read (37,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%internals%structQuad8IDs(1:rowcount))
        do ii=1,rowcount
            read (37,'(I10)',iostat=read_error) fesim%internals%structQuad8IDs(ii)
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (37, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if


    fname='structelem2aeronode3d.fipps'
    open (unit=37,file=fname,status='old',action='read', iostat=io_error)
    if (io_error == 0) then
        read (37,*,iostat=read_error) rowcount
        if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was a problem reading the number of rows of the file'
            err_code=1
            goto 9999
        end if
        allocate(fesim%internals%structElem2aeroNode(1:rowcount))
        do ii=1,rowcount
            read (37,'(2I10,2E23.16)',iostat=read_error) &
                & fesim%internals%structElem2aeroNode(ii)%nodeID, &
                & fesim%internals%structElem2aeroNode(ii)%elemID, &
                & fesim%internals%structElem2aeroNode(ii)%xi,     &
                & fesim%internals%structElem2aeroNode(ii)%eta
            if (read_error /= 0) then
                write(*,*) 'Error reading file ', fname
                write(*,*) 'There was an error processing line', ii
                err_code=1
                goto 9999
            end if
        end do
        close (37, iostat=io_error)
        if (io_error /= 0) then
            write(*,*) 'Error closing file ', fname
            err_code=1
            goto 9999
        end if
    else
        write (*,*) 'Error opening file ', fname
        err_code=1
        goto 9999
    end if
end if

if (fesim%is_thermal) then
    fname = 'thermal.fipps'
    open (unit=37,file=fname,status='old',action='read', iostat=io_error)
    if (read_error /= 0) then
        write(*,*) 'Error reading file ', fname
        write(*,*) 'There was a problem reading the number of rows of the file'
        err_code=1
        goto 9999
    end if
    ! Hier noch die ganzen thermal_type-Variablen auslesen
end if
!hier die weiteren Randbedingungen, Materialien einlesen
!diese weiteren karten auch mit der input_process_line-subroutine aus der control auslesen lassen
!deren is_... variablen in fesimulation_typen setzen

!
! =================================================================================================
!
! Error handling
!

9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'input_tf'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

contains

  subroutine read_contact_beam2_node(unit, rowcount)
  
    implicit none
    
    integer, intent(in)                         :: unit, rowcount
    
    type nodeListItem
      integer                                   :: nid                                  ! Node-ID
      double precision                          :: xi                                   ! Elementkoordinate
      type(nodeListItem), pointer               :: next => null()
    end type nodeListItem
    
    type(nodeListItem), pointer                 :: actItem,tmpItem
    
    type nodeListItemArray
        type(nodeListItem), pointer             :: p
    end type nodeListItemArray
    
    type(nodeListItemArray), allocatable        :: beam2ConNodeList(:)
    
    integer                                     :: nodeID, beam2ID, conNum 
    integer                                     :: ii,jj,anz,numContactedElems
    double precision                            :: xi
    
    allocate(beam2ConNodeList(1:size(fesim%elemente%beam2s,1)))
    
    do ii = 1,size(fesim%elemente%beam2s,1)
      beam2ConNodeList(ii)%p => null()
    end do
    actItem => null()
    
    numContactedElems = 0
    
    do ii=1,rowcount
      read (unit,'(2I10,E23.16)',iostat=read_error) beam2ID, nodeID, xi
      if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was an error processing line', ii
            err_code=1
            goto 9999
      end if
      allocate(actItem)
      actItem%nid = nodeID
      actItem%xi  = xi
      if (associated(beam2ConNodeList(beam2ID)%p)) then
        actItem%next => beam2ConNodeList(beam2ID)%p
      else
        numContactedElems = numContactedElems + 1
      end if
      
      beam2ConNodeList(beam2ID)%p => actItem

      actItem => null()
    end do
    
    allocate(fesim%randbedingungen%contact_node_beam2(1:numContactedElems))
    conNum = 0
    
    do ii = 1,size(fesim%elemente%beam2s,1)
      if (associated(beam2ConNodeList(ii)%p) .eq. .FALSE.) cycle
      anz = 1
      actItem => beam2ConNodeList(ii)%p
      do while (associated(actItem%next))
        anz = anz + 1
        actItem => actItem%next
      end do
      conNum = conNum + 1
      fesim%randbedingungen%contact_node_beam2(conNum)%beam2ID = ii
      allocate(fesim%randbedingungen%contact_node_beam2(conNum)%nodeIDs(anz))
      allocate(fesim%randbedingungen%contact_node_beam2(conNum)%xi(anz))
      actItem => beam2ConNodeList(ii)%p
      jj = 0
      do while (associated(actItem))
        jj = jj + 1
        fesim%randbedingungen%contact_node_beam2(conNum)%nodeIDs(jj) = actItem%nid
        fesim%randbedingungen%contact_node_beam2(conNum)%xi(jj)      = actItem%xi
        tmpItem => actItem
        actItem => actItem%next
        deallocate(tmpItem)
      end do
    end do
    
    deallocate(beam2ConNodeList)

  !
  ! =================================================================================================
  !
  ! Error handling
  !
  
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)                      'An error occured in subroutine'
     write(*,*)                      'read_contact_beam2_node'
     write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
     write(*,*)                      'exit program '
     stop
     
  end if
      
  end subroutine read_contact_beam2_node

  subroutine read_contact_quad8_node(unit, rowcount)
  
    implicit none
    
    integer, intent(in)                         :: unit, rowcount
    
    type nodeListItem
      integer                                   :: nid                                  ! Node-ID
      logical                                   :: kpt                                  ! Kontakttyp
      type(nodeListItem), pointer               :: next => null()
    end type nodeListItem
    
    type(nodeListItem), pointer                 :: actItem,tmpItem
    
    type nodeListItemArray
        type(nodeListItem), pointer             :: p
    end type nodeListItemArray
    
    type(nodeListItemArray), allocatable        :: quad8ConNodeList(:)
    
    integer                                     :: nodeID, quad8ID, conNum 
    integer                                     :: ii,jj,anz,numContactedElems
    
    allocate(quad8ConNodeList(1:size(fesim%elemente%quad8s,1)))
    
    do ii = 1,size(fesim%elemente%quad8s,1)
      quad8ConNodeList(ii)%p => null()
    end do
    actItem => null()
    
    numContactedElems = 0
    
    do ii=1,rowcount
      read (unit,'(2I10)',iostat=read_error) quad8ID, nodeID
      if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was an error processing line', ii
            err_code=1
            goto 9999
      end if
      allocate(actItem)
      actItem%nid = nodeID
      if (associated(quad8ConNodeList(quad8ID)%p)) then
        actItem%next => quad8ConNodeList(quad8ID)%p
      else
        numContactedElems = numContactedElems + 1
      end if
      
      quad8ConNodeList(quad8ID)%p => actItem

      actItem => null()
    end do
    
    allocate(fesim%randbedingungen%contact_node_quad8(1:numContactedElems))
    conNum = 0
    
    do ii = 1,size(fesim%elemente%quad8s,1)
      if (associated(quad8ConNodeList(ii)%p) .eq. .FALSE.) cycle
      anz = 1
      actItem => quad8ConNodeList(ii)%p
      do while (associated(actItem%next))
        anz = anz + 1
        actItem => actItem%next
      end do
      conNum = conNum + 1
      fesim%randbedingungen%contact_node_quad8(conNum)%quad8ID = ii
      allocate(fesim%randbedingungen%contact_node_quad8(conNum)%nodeIDs(anz))
      actItem => quad8ConNodeList(ii)%p
      jj = 0
      do while (associated(actItem))
        jj = jj + 1
        fesim%randbedingungen%contact_node_quad8(conNum)%nodeIDs(jj) = actItem%nid
        tmpItem => actItem
        actItem => actItem%next
        deallocate(tmpItem)
      end do
    end do
    
    deallocate(quad8ConNodeList)

  !
  ! =================================================================================================
  !
  ! Error handling
  !
  
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)                      'An error occured in subroutine'
     write(*,*)                      'read_contact_quad8_node'
     write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
     write(*,*)                      'exit program '
     stop
     
  end if
      
  end subroutine read_contact_quad8_node

  subroutine read_contact_lsolid20_node(unit, rowcount)
  
    implicit none
    
    integer, intent(in)                         :: unit, rowcount
    
    type nodeListItem
      integer                                   :: nid                                  ! Node-ID
      logical                                   :: kpt                                  ! Kontakttyp
      type(nodeListItem), pointer               :: next => null()
    end type nodeListItem
    
    type(nodeListItem), pointer                 :: actItem,tmpItem
    
    type nodeListItemArray
        type(nodeListItem), pointer             :: p
    end type nodeListItemArray
    
    type(nodeListItemArray), allocatable        :: lsolid20ConNodeList(:)
    
    integer                                     :: nodeID, lsolid20ID, conNum 
    integer                                     :: ii,jj,anz,numContactedElems
    
    allocate(lsolid20ConNodeList(1:size(fesim%elemente%lsolid20s,1)))
    
    do ii = 1,size(fesim%elemente%lsolid20s,1)
      lsolid20ConNodeList(ii)%p => null()
    end do
    actItem => null()
    
    numContactedElems = 0
    
    do ii=1,rowcount
      read (unit,'(2I10)',iostat=read_error) lsolid20ID, nodeID
      if (read_error /= 0) then
            write(*,*) 'Error reading file ', fname
            write(*,*) 'There was an error processing line', ii
            err_code=1
            goto 9999
      end if
      allocate(actItem)
      actItem%nid = nodeID
      if (associated(lsolid20ConNodeList(lsolid20ID)%p)) then
        actItem%next => lsolid20ConNodeList(lsolid20ID)%p
      else
        numContactedElems = numContactedElems + 1
      end if
      
      lsolid20ConNodeList(lsolid20ID)%p => actItem

      actItem => null()
    end do
    
    allocate(fesim%randbedingungen%contact_node_lsolid20(1:numContactedElems))
    conNum = 0
    
    do ii = 1,size(fesim%elemente%lsolid20s,1)
      if (associated(lsolid20ConNodeList(ii)%p) .eq. .FALSE.) cycle
      anz = 1
      actItem => lsolid20ConNodeList(ii)%p
      do while (associated(actItem%next))
        anz = anz + 1
        actItem => actItem%next
      end do
      conNum = conNum + 1
      fesim%randbedingungen%contact_node_lsolid20(conNum)%lsolid20ID = ii
      allocate(fesim%randbedingungen%contact_node_lsolid20(conNum)%nodeIDs(anz))
      actItem => lsolid20ConNodeList(ii)%p
      jj = 0
      do while (associated(actItem))
        jj = jj + 1
        fesim%randbedingungen%contact_node_lsolid20(conNum)%nodeIDs(jj) = actItem%nid
        tmpItem => actItem
        actItem => actItem%next
        deallocate(tmpItem)
      end do
    end do
    
    deallocate(lsolid20ConNodeList)

  !
  ! =================================================================================================
  !
  ! Error handling
  !
  
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)                      'An error occured in subroutine'
     write(*,*)                      'read_contact_lsolid20_node'
     write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
     write(*,*)                      'exit program '
     stop
     
  end if
      
  end subroutine read_contact_lsolid20_node

end subroutine input_tf
  
  
  
  
  

  

  
  
  
  
