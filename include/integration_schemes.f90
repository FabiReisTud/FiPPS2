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
!> $Id: integration_schemes.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module integration_schemes

    implicit none

    contains

! =================================================================================================
!
!> @brief
!> Liefert alle notwendigen Daten für eine Integration in der Fläche
!
!> @details
!> Gibt die Koordinaten der Stützenstellen und die dazugehörigen Gewichte
!> für je eine Gauss- (1,2,3 Stützstellen) oder eine Simpson-Integration
!> (1,3,5,7,9 Stützstellen) für beide Koordinaten zurück
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter, 03.05.2017
!
! =================================================================================================

    subroutine integration_points_2d(nPcoord1, intType1, coord1, wCoord1, &
                                     nPcoord2, intType2, coord2, wCoord2)

      implicit none

      integer, intent(in)                                         :: nPcoord1           !< Anzahl der Integrationspunkte in Richtung der ersten Koordinate
      integer, intent(in)                                         :: nPcoord2           !< Anzahl der Integrationspunkte in Richtung der zweiten Koordinate
      character(len=5), intent(in)                                :: intType1           !< Integrationstype ("Gauss", "Simps") in Richtung der ersten Koordinate
      character(len=5), intent(in)                                :: intType2           !< Integrationstype ("Gauss", "Simps") in Richtung der zweiten Koordinate

      double precision, dimension(nPcoord1*nPcoord2), intent(out) :: coord1             !< Relative Koordinaten der Integrationspunkte in Richtung der ersten Koordinate (-1,1)
      double precision, dimension(nPcoord1*nPcoord2), intent(out) :: coord2             !< Relative Koordinaten der Integrationspunkte in Richtung der zweiten Koordinate (-1,1)
      double precision, dimension(nPcoord1*nPcoord2), intent(out) :: wCoord1            !< Gewichte der Integrationspunkte in Richtung der ersten Koordinate
      double precision, dimension(nPcoord1*nPcoord2), intent(out) :: wCoord2            !< Gewichte der Integrationspunkte in Richtung der zweiten Koordinate

      double precision, dimension(nPcoord1)                       :: coord1_t, wCoord1_t
      double precision, dimension(nPcoord2)                       :: coord2_t, wCoord2_t

      integer                                                     :: index, i1, i2

      if (intType1 .EQ. "Gauss") then
        call gauss_integration(nPcoord1, coord1_t, wCoord1_t)
      else if (intType1 .EQ. "Simps") then
        call simpson_integration(nPcoord1, coord1_t, wCoord1_t)
      else
        write(*,*) "No integration type found."
        STOP
      end if

      if (intType2 .EQ. "Gauss") then
        call gauss_integration(nPcoord2, coord2_t, wCoord2_t)
      else if (intType2 .EQ. "Simps") then
        call simpson_integration(nPcoord2, coord2_t, wCoord2_t)
      else
        write(*,*) "No integration type found."
        STOP
      end if

      index = 1
      do i2 = 1, nPcoord2
        do i1 = 1, nPcoord1
          coord1(index)  = coord1_t(i1)
          wCoord1(index) = wCoord1_t(i1)

          coord2(index)  = coord2_t(i2)
          wCoord2(index) = wCoord2_t(i2)
          index = index + 1
        end do
      end do

    end subroutine integration_points_2d

! =================================================================================================
!
!> @brief
!> Liefert alle notwendigen Daten für eine Integration im Volumen (3 Koordinaten)
!
!> @details
!> Gibt die Koordinaten der Stützenstellen und die dazugehörigen Gewichte
!> für je eine Gauss- (1,2,3 Stützstellen) oder eine Simpson-Integration
!> (1,3,5,7,9 Stützstellen) für alle drei Koordinaten zurück
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter, 03.05.2017
!
! =================================================================================================

    subroutine integration_points_3d(nPcoord1, intType1, coord1, wCoord1, &
                                     nPcoord2, intType2, coord2, wCoord2, &
                                     nPcoord3, intType3, coord3, wCoord3)

      implicit none

      integer, intent(in)                                         :: nPcoord1 !> Anzahl der Integrationspunkte in Richtung der ersten Koordinate
      integer, intent(in)                                         :: nPcoord2 !> Anzahl der Integrationspunkte in Richtung der zweiten Koordinate
      integer, intent(in)                                         :: nPcoord3 !> Anzahl der Integrationspunkte in Richtung der dritten Koordinate
      character(len=5), intent(in)                                :: intType1 !> Integrationstype ("Gauss", "Simps") in Richtung der ersten Koordinate
      character(len=5), intent(in)                                :: intType2 !> Integrationstype ("Gauss", "Simps") in Richtung der zweiten Koordinate
      character(len=5), intent(in)                                :: intType3 !> Integrationstype ("Gauss", "Simps") in Richtung der dritten Koordinate

      double precision, dimension(nPcoord1*nPcoord2*nPcoord3), intent(out) :: coord1 !> Relative Koordinate der Integrationspunkte in Richtung der ersten Koordinate (-1,1)
      double precision, dimension(nPcoord1*nPcoord2*nPcoord3), intent(out) :: coord2 !> Relative Koordinate der Integrationspunkte in Richtung der zweiten Koordinate (-1,1)
      double precision, dimension(nPcoord1*nPcoord2*nPcoord3), intent(out) :: coord3 !> Relative Koordinate der Integrationspunkte in Richtung der dritten Koordinate (-1,1)
      double precision, dimension(nPcoord1*nPcoord2*nPcoord3), intent(out) :: wCoord1 !> Gewichte der Integrationspunkte in Richtung der ersten Koordinate
      double precision, dimension(nPcoord1*nPcoord2*nPcoord3), intent(out) :: wCoord2 !> Gewichte der Integrationspunkte in Richtung der zweiten Koordinate
      double precision, dimension(nPcoord1*nPcoord2*nPcoord3), intent(out) :: wCoord3 !> Gewichte der Integrationspunkte in Richtung der dritten Koordinate

      double precision, dimension(nPcoord1)                       :: coord1_t, wCoord1_t
      double precision, dimension(nPcoord2)                       :: coord2_t, wCoord2_t
      double precision, dimension(nPcoord3)                       :: coord3_t, wCoord3_t

      integer                                                     :: index, i1, i2, i3

      if (intType1 .EQ. "Gauss") then
        call gauss_integration(nPcoord1, coord1_t, wCoord1_t)
      else if (intType1 .EQ. "Simps") then
        call simpson_integration(nPcoord1, coord1_t, wCoord1_t)
      else
        write(*,*) "No integration type found."
        STOP
      end if

      if (intType2 .EQ. "Gauss") then
        call gauss_integration(nPcoord2, coord2_t, wCoord2_t)
      else if (intType2 .EQ. "Simps") then
        call simpson_integration(nPcoord2, coord2_t, wCoord2_t)
      else
        write(*,*) "No integration type found."
        STOP
      end if

      if (intType3 .EQ. "Gauss") then
        call gauss_integration(nPcoord3, coord3_t, wCoord3_t)
      else if (intType3 .EQ. "Simps") then
        call simpson_integration(nPcoord3, coord3_t, wCoord3_t)
      else
        write(*,*) "No integration type found."
        STOP
      end if

      index = 1
      do i3 = 1, nPcoord3
        do i2 = 1, nPcoord2
          do i1 = 1, nPcoord1
            coord1(index)  = coord1_t(i1)
            wCoord1(index) = wCoord1_t(i1)

            coord2(index)  = coord2_t(i2)
            wCoord2(index) = wCoord2_t(i2)

            coord3(index)  = coord3_t(i3)
            wCoord3(index) = wCoord3_t(i3)
            index = index + 1
          end do
        end do
      end do

    end subroutine integration_points_3d

! =================================================================================================
!
!> @brief
!> Liefert alle notwendigen Daten für eine Gausspunktintegration
!
!> @details
!> Gibt die Koordinaten der Stützenstellen und die dazugehörigen Gewichte
!> für eine 1,2 oder 3- Gausspunktintragration zurück
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter, 16.06.2010
!
! =================================================================================================

    subroutine gauss_integration(numPoints, coord, weight)

    !
    ! use
    !

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
    integer, intent(in)                                           :: numPoints
    integer                                                       :: err_code = 0

    double precision                                              :: facDP
    double precision, dimension(numPoints*numPoints), intent(out) :: coord
    double precision, dimension(numPoints*numPoints), intent(out) :: weight
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
    if (numPoints == 3) then

        ! Build vectors with Gauß-point coordinates in xi and eta

        facDP = sqrt(0.6D0)

        coord(1) = -facDP
        coord(2) =  0.D0
        coord(3) =  facDP

        ! Build vectors with Gauss-integration weights

        weight(1) = 5.D0/9.D0
        weight(2) = 8.D0/9.D0
        weight(3) = weight(1)

    else if (numPoints == 2) then

        facDP = 1.D0/sqrt(3.D0)

        coord(1) = -facDP
        coord(2) =  facDP

        weight(1)  = 1.D0
        weight(2)  = weight(1)

    else if (numPoints == 1) then

        facDP      = 0.D0

        coord(1)  = facDP

        weight(1)    = 2.D0

    else

        write(*,*) 'wrong input on parameter numPoints'
        write(*,*) 'values 1, 2 or 3 are allowed to be chosen'
        err_code = 1
        goto 9999

    end if
    !
    ! =================================================================================================
    !
    ! Error handling
    !
    9999 continue

    if (err_code /= 0) then

        write(*,*)                      'An error occured in subroutine'
        write(*,*)                      'gauss_integration'
        write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
        write(*,*)                      'exit program '
        stop

    end if

    end subroutine gauss_integration

! =================================================================================================
!
!> @brief
!> Calculate integration points and weights with respect to the Simpson integration
!
!> @details
!> Returns the coordinates of the integration points and the corredponding weights
!> for the Simpson's integration rule for 9, 7, 5, 3 or 1 points
!
!> @author Florian Dexl, TU Dresden, wissenschaftlicher Mitarbeiter, ??.??.20??
!
! =================================================================================================

    subroutine simpson_integration(numPoints, coord, weight)
    !
    !
    ! use
    !
    !
    ! =================================================================================================
    !
    implicit none
    !
    ! =================================================================================================
    !
    ! Data types
    !
    ! Input
    !
    integer, intent(in)                                    :: numPoints !< number of integration points
    !
    ! Output
    !
    double precision, dimension(numPoints), intent(out)    :: coord     !< array containing coordinates of integration points
    double precision, dimension(numPoints), intent(out)    :: weight    !< weights needed for integration
    !
    ! Internal
    !
    integer                                                :: err_code=0
    !
    ! =================================================================================================
    !
    ! Calculation
    !

    ! Get through thickness integration points by Simpson's rule
    if (numPoints == 1) then

        coord = 0.d0

        weight   = 2.d0

    else if (numPoints == 3) then

        coord(1) = -1.d0
        coord(2) =  0.d0
        coord(3) =  1.d0

        weight(1) = 1.d0/3.d0
        weight(2) = 4.d0/3.d0
        weight(3) = 1.d0/3.d0

    else if (numPoints == 5) then

        coord(1) = -1.d0
        coord(2) = -0.5d0
        coord(3) =  0.d0
        coord(4) =  0.5d0
        coord(5) =  1.d0

        weight(1) = 1.d0/6.d0
        weight(2) = 4.d0/6.d0
        weight(3) = 2.d0/6.d0
        weight(4) = 4.d0/6.d0
        weight(5) = 1.d0/6.d0

    else if (numPoints == 7) then

        coord(1) = -1.d0
        coord(2) = -2.d0/3.d0
        coord(3) = -1.d0/3.d0
        coord(4) =  0.d0
        coord(5) =  1.d0/3.d0
        coord(6) =  2.d0/3.d0
        coord(7) =  1.d0

        weight(1) = 1.d0/9.d0
        weight(2) = 4.d0/9.d0
        weight(3) = 2.d0/9.d0
        weight(4) = 4.d0/9.d0
        weight(5) = 2.d0/9.d0
        weight(6) = 4.d0/9.d0
        weight(7) = 1.d0/9.d0

    else if (numPoints == 9) then

        coord(1) = -1.d0
        coord(2) = -3.d0/4.d0
        coord(3) = -2.d0/4.d0
        coord(4) = -1.d0/4.d0
        coord(5) =  0.d0
        coord(6) =  1.d0/4.d0
        coord(7) =  2.d0/4.d0
        coord(8) =  3.d0/4.d0
        coord(9) =  1.d0

        weight(1) = 1.d0/12.d0
        weight(2) = 4.d0/12.d0
        weight(3) = 2.d0/12.d0
        weight(4) = 4.d0/12.d0
        weight(5) = 2.d0/12.d0
        weight(6) = 4.d0/12.d0
        weight(7) = 2.d0/12.d0
        weight(8) = 4.d0/12.d0
        weight(9) = 1.d0/12.d0

    else

        write(*,*) 'wrong input on parameter numPoints (number of'
        write(*,*) 'integration points)'
        write(*,*) 'values 1, 3, 5, 7 or 9 are allowed to be chosen'
        err_code = 2
        goto 9999

    end if

    !
    ! =================================================================================================
    !
    ! Error handling
    !
    9999 continue
    !
    if (err_code /= 0) then
        !
        write(*,*)                      'An error occured in subroutine'
        write(*,*)                      'simpson_integration'
        write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
        write(*,*)                      'exit program '
        stop
        !
    end if
    !
    return
    !
    end subroutine simpson_integration

! =================================================================================================
!
!> @brief
!> Delivers 14 integration points for fully integrated brick elements
!
!> @details
!> Delivers 14 integration points for fully integrated brick elements
!> Not finished yet
!
!> @author Fabian Reisener, TU Dresden, wissenschaftlicher Mitarbeiter
!
! =================================================================================================
    subroutine integration_points_brick14(coord1,wcoord1,coord2,wcoord2,coord3,wcoord3)
        double precision, dimension(14), intent(out) :: coord1,coord2,coord3    !> Relative Koordinate der Integrationspunkte in Richtung der 1.,2.,3. Koordinate (-1,1)
        double precision, dimension(14), intent(out) :: wcoord1,wcoord2,wcoord3 !> Gewichtung der Integrationspunkte in Richtung der 1.,2.,3. Koordinate


    end subroutine
end module integration_schemes
