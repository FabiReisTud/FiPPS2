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
!> Berechnung der Fluid-Struktur-Interaktion in 2D oder 3D
!
!> @details
!> Ansteuerung eines zwei- oder dreidimensionalen Stroemungsloesers zur Berechnung der Fluid-
!> Struktur-Interaktion.
!
!> @author Florian Dexl, TU Dresden, WiMi, 10.09.2018
!
!> $Id: aero_coupling.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine aero_coupling (Ausgabe, fesim, scloop, cdi)

use globale_variablen, only : apamePressureOutFile
use fesimulation_typen
use konstanten
!
! =================================================================================================
!
implicit none

interface
    subroutine apame_main(fname,readmodel,displacements,dispfac,pressures,F_total,cdi,output)
        character(len=*), optional                                            :: fname             !< Filename of the apame input file
        logical,                                       intent(in), optional   :: readmodel         !< Flag, to read the input file
        double precision, dimension(:,:),              intent(in), optional   :: displacements     !< Displacements to be added on node coordinates
        double precision,                              intent(in), optional   :: dispfac           !< Factor with which displacements are added
        double precision, dimension(:,:), allocatable, intent(out), optional  :: pressures         !< Manometer pressure for all panels and cases
        double precision, dimension(:), allocatable,   intent(out), optional  :: F_total           !< Sum of pressure times panel area
        double precision, dimension(:), allocatable,   intent(out), optional  :: cdi               !< induzierter Widerstand
        logical,                                       intent(in), optional   :: output            !< Flag, if APAME should write any output
    end subroutine apame_main

    subroutine panel2d_main(fname,readmodel,inviscid,displacements,dispfac,given_dalpha,paneling,pressures,F_total,al,ca,cd,cm,hkmax,xsep,Mcrit,output)
        character(len=*)                                                      :: fname             !< Filename of the PANEL2D input file
        logical,                                       intent(in), optional   :: readmodel         !< Flag, to read the input file
        logical,                                       intent(in), optional   :: inviscid          !< Flag, to force inviscid calculation
        double precision, dimension(:,:),              intent(in), optional   :: displacements     !< Displacements to be added on node coordinates
        double precision,                              intent(in), optional   :: dispfac           !< Factor with which displacements are added
        double precision,                              intent(in), optional   :: given_dalpha      !< Delta Alpha (overrides value from inputfile)
        logical,                                       intent(in), optional   :: paneling          !< Flag, to enable repaneling in Xfoil
        double precision, dimension(:,:), allocatable, intent(out), optional  :: pressures         !< Manometer pressure for all panels and cases
        double precision, dimension(:), allocatable,   intent(out), optional  :: F_total           !< Sum of pressure times panel area
        double precision, dimension(:), allocatable,   intent(out), optional  :: al                !< Anstellwinkel
        double precision, dimension(:), allocatable,   intent(out), optional  :: ca                !< Auftriebsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cd                !< Widerstandsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cm                !< Momentenbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: hkmax             !< Maximaler kinematischer Formfaktor
        double precision, dimension(:), allocatable,   intent(out), optional  :: xsep              !< Position der turbulenten Abloesung
        double precision, dimension(:), allocatable,   intent(out), optional  :: Mcrit             !< Kritische MACH-Zahl
        logical,                                       intent(in), optional   :: output            !< Flag, if PANEL2D should write any output
    end subroutine panel2d_main

    subroutine xfoilwrapper_main(fname,readmodel,inviscid,displacements,dispfac,given_dalpha,paneling,pressures,F_total,al,ca,cd,cm,hkmax,Mcrit,output)
        character(len=*)                                                      :: fname             !< Filename of the XFOIL_WRAPPER input file
        logical,                                       intent(in), optional   :: readmodel         !< Flag, to read the input file
        logical,                                       intent(in), optional   :: inviscid          !< Flag, to force inviscid calculation
        double precision, dimension(:,:),              intent(in), optional   :: displacements     !< Displacements to be added on node coordinates
        double precision,                              intent(in), optional   :: dispfac           !< Factor with which displacements are added
        double precision,                              intent(in), optional   :: given_dalpha      !< Delta Alpha (overrides value from inputfile)
        logical,                                       intent(in), optional   :: paneling          !< Flag, to enable repaneling in Xfoil
        double precision, dimension(:,:), allocatable, intent(out), optional  :: pressures         !< Manometer pressure for all panels and cases
        double precision, dimension(:), allocatable,   intent(out), optional  :: F_total           !< Sum of pressure times panel area
        double precision, dimension(:), allocatable,   intent(out), optional  :: al                !< Anstellwinkel
        double precision, dimension(:), allocatable,   intent(out), optional  :: ca                !< Auftriebsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cd                !< Widerstandsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cm                !< Momentenbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: hkmax             !< Maximaler kinematischer Formfaktor
        double precision, dimension(:), allocatable,   intent(out), optional  :: Mcrit             !< Kritische MACH-Zahl
        logical,                                       intent(in), optional   :: output            !< Flag, if PANEL2D should write any output
    end subroutine xfoilwrapper_main

    subroutine load_add_aeropressure(fesim, scloop, pAeroElem, F_total_aero)
        use fesimulation_typen
        use konstanten
        implicit none
        type(fe_simulation)                                 :: fesim
        integer, intent(in)                                 :: scloop
        double precision, dimension(:,:), intent(in)        :: pAeroElem
        double precision, dimension(:), intent(in)          :: F_total_aero
    end subroutine load_add_aeropressure
end interface
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
logical, intent(in)                                               :: Ausgabe
type(fe_simulation), intent(in)                                   :: fesim
integer, intent(in)                                               :: scloop

!
! Output
!
double precision, allocatable, dimension(:), intent(out)          :: cdi

!
! inner
!
character(len=56)                                                 :: filename

double precision, allocatable, dimension(:,:)                     :: pressures
double precision, allocatable, dimension(:)                       :: F_total

logical                                                           :: found_aeroload

integer                                                           :: aeroloadID
integer                                                           :: err_code=0

integer                                                           :: ii
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

found_aeroload = .false.

aeroloadID = fesim%lasten%subcases(scloop)%aeroloadID

if (aeroloadID .ne. 0) then

  if (fesim%is_aeroload3d .eq. .true.) then

    if (Ausgabe .eq. .true.) write(*,*) 'START - Aerodynamische Rechnung'

    if (fesim%lasten%aeroload3ds(aeroloadID)%mthd .eq. 1) then
      write(filename,'(A9,I4.4)') "apame_sc_", scloop

      if (allocated(fesim%internals%aeroDispl) .eq. .true.) then
        call apame_main(fname=filename,readmodel=.false.,displacements=fesim%internals%aeroDispl, pressures=pressures, F_total=F_total, cdi=cdi, output=fesim%ausgabe%outputVTK)
      else
        write(*,*) 'Originalmodell'
        call apame_main(fname=filename,readmodel=.true.,pressures=pressures,F_total=F_total,cdi=cdi,output=fesim%ausgabe%outputVTK)
      end if

      if (fesim%ausgabe%outputApamePressures .eq. .true.) then
        open (apamePressureOutFile, file = trim(filename) // '_pressure.res', status='UNKNOWN')
        write(apamePressureOutFile,'(E27.17)') F_total(1)
        write(apamePressureOutFile,'(E27.17)') cdi(1)
        write(apamePressureOutFile,'(I10)') size(pressures,1)
        do ii = 1, size(pressures,1)
          write(apamePressureOutFile,'(E27.17)') pressures(ii,1)
        end do
        close(apamePressureOutFile)
      end if
    else if (fesim%lasten%aeroload3ds(aeroloadID)%mthd .eq. 2) then
      write(filename,'(A9,I4.4)') "apame_sc_", scloop
      open (apamePressureOutFile, file = trim(filename) // '_pressure.res', status='old')
      allocate(F_total(1),cdi(1))
      read(apamePressureOutFile,'(E27.17)') F_total(1)
      read(apamePressureOutFile,'(E27.17)') cdi(1)
      read(apamePressureOutFile,'(I10)') ii
      allocate(pressures(ii,1))
      do ii = 1, size(pressures,1)
        read(apamePressureOutFile,'(E27.17)') pressures(ii,1)
      end do
      close(apamePressureOutFile)
    else
      write(*,*) 'Unsupported method for aeroload3d!'
      write(*,*) 'Valid methods are:'
      write(*,*) '  1 - APAME'
      write(*,*) '  2 - APAME saved results'
      STOP
    end if

    if (Ausgabe .eq. .true.) write(*,*) 'START - Hinzufügen von Aero-Lasten'

    call load_add_aeropressure(fesim, scloop, pressures, F_total)

    if (Ausgabe .eq. .true.) write(*,*) 'ENDE  - Hinzufügen von Aero-Lasten'

    if (Ausgabe .eq. .true.) write(*,*) 'ENDE  - Aerodynamische Rechnung'

  end if

  if (fesim%is_aeroload2d .eq. .true.) then

    if (Ausgabe .eq. .true.) write(*,*) 'START - Aerodynamische Rechnung'

    if (fesim%lasten%aeroload2ds(aeroloadID)%mthd .eq. 1) then
      write(filename,'(A11,I4.4)') "panel2d_sc_", scloop

      if (allocated(fesim%internals%aeroDispl) .eq. .true.) then
        call  panel2d_main(fname=filename,readmodel=.false.,inviscid=.true.,displacements=fesim%internals%aeroDispl,pressures=pressures,F_total=F_total,paneling=.false.,cd=cdi,output=fesim%ausgabe%outputVTK)
      else
        write(*,*) 'Originalmodell'
        call panel2d_main(fname=filename,readmodel=.true.,inviscid=.true.,pressures=pressures,F_total=F_total,paneling=.false.,cd=cdi,output=fesim%ausgabe%outputVTK)
      end if
    else if (fesim%lasten%aeroload2ds(aeroloadID)%mthd .eq. 2) then
      write(filename,'(A11,I4.4)') "panel2d_sc_", scloop

      if (allocated(fesim%internals%aeroDispl) .eq. .true.) then
        call xfoilwrapper_main(fname=filename,readmodel=.false.,inviscid=.false.,displacements=fesim%internals%aeroDispl,pressures=pressures,F_total=F_total,paneling=.false.,cd=cdi,output=fesim%ausgabe%outputVTK)
      else
        write(*,*) 'Originalmodell'
        call xfoilwrapper_main(fname=filename,readmodel=.true.,inviscid=.true.,pressures=pressures,F_total=F_total,paneling=.false.,cd=cdi,output=fesim%ausgabe%outputVTK)
      end if
    else if (fesim%lasten%aeroload2ds(aeroloadID)%mthd .eq. 3) then
      write(filename,'(A11,I4.4)') "panel2d_sc_", scloop

      if (allocated(fesim%internals%aeroDispl) .eq. .true.) then
        call xfoilwrapper_main(fname=filename,readmodel=.false.,inviscid=.true.,displacements=fesim%internals%aeroDispl,pressures=pressures,F_total=F_total,paneling=.false.,cd=cdi,output=fesim%ausgabe%outputVTK)
      else
        write(*,*) 'Originalmodell'
        call xfoilwrapper_main(fname=filename,readmodel=.true.,inviscid=.true.,pressures=pressures,F_total=F_total,paneling=.false.,cd=cdi,output=fesim%ausgabe%outputVTK)
      end if
    else
      write(*,*) 'Unsupported method for aeroload2d!'
      write(*,*) 'Valid methods are:'
      write(*,*) '  1 - PANEL2D'
      write(*,*) '  2 - XFOIL viscous'
      write(*,*) '  3 - XFOIL inviscid'
      STOP
    end if

    if (Ausgabe .eq. .true.) write(*,*) 'START - Hinzufügen von Aero-Lasten'

    call load_add_aeropressure(fesim, scloop, pressures, F_total)

    if (Ausgabe .eq. .true.) write(*,*) 'ENDE  - Hinzufügen von Aero-Lasten'

    if (Ausgabe .eq. .true.) write(*,*) 'ENDE  - Aerodynamische Rechnung'

  end if

end if

if (allocated(pressures)) deallocate(pressures)
if (allocated(F_total))   deallocate(F_total)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'aero_coupling'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

return

end subroutine aero_coupling
