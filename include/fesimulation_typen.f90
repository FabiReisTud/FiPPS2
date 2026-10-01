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
!> @brief
!
!> @details
!
!> @author
!
!> $Id: fesimulation_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module fesimulation_typen

use knoten_typen
use eigenschaften_typen
use element_typen
use failurecriteria_typen
use koordinatensystem_typen
use laminat_typen
use last_typen
use material_typen
use randbedingung_typen
use residual_typen
use internals_typen
use ausgabe_typen
use ergebnis_typen
use thermal_typen

implicit none
  
type fe_simulation
  
  ! Notwendige Daten einer FE-Simulation
  type(nodes_type)              :: knoten
  type(elemente_type)           :: elemente
  type(eigenschaften_type)      :: eigenschaften
  type(failurecriteria_type)    :: versagenskriterien
  type(coordsyses_type)         :: koordinatensysteme
  type(laminats_type)           :: laminate
  type(lasten_type)             :: lasten
  type(materials_type)          :: materialien
  type(randbedinungs_type)      :: randbedingungen
  type(residuals_type)          :: residuals
  type(internals_type)          :: internals
  type(ausgabe_type)            :: ausgabe
  type(ergebnis_type)           :: ergebnisse
  type(thermal_state_type)      :: thermal ! Enthält den Zustand der Wärmeleitungslösung

  ! Flags, ob die verschiedenen Punkte ursprünglich eingelesen/angegeben wurden
  logical                       :: is_node
  logical                       :: is_beam2
  logical                       :: is_quad8
  logical                       :: is_lsolid20
  logical                       :: is_load
  logical                       :: is_force
  logical                       :: is_moment
  logical                       :: is_p2load
  logical                       :: is_p8load
  logical                       :: is_p20load
  logical                       :: is_aeroload2d
  logical                       :: is_aeroload3d
  logical                       :: is_mat1
  logical                       :: is_mat8
  logical                       :: is_mat20
  logical                       :: is_pbeam
  logical                       :: is_pshell
  logical                       :: is_pcomp
  logical                       :: is_plsolid
  logical                       :: is_spcadd
  logical                       :: is_spc1
  logical                       :: is_spcd
  logical                       :: is_mpcadd
  logical                       :: is_mpc
  logical                       :: is_temperature
  logical                       :: is_beam2temp
  logical                       :: is_quad8temp
  logical                       :: is_lsolid20temp
  logical                       :: is_subcase
  logical                       :: is_lam8
  logical                       :: is_lam20
  logical                       :: is_coupling
  logical                       :: is_contact_node_beam2
  logical                       :: is_contact_node_quad8
  logical                       :: is_contact_node_lsolid20
  logical                       :: is_failure
  logical                       :: is_coord
  logical                       :: is_multistep
  ! Thermal: Start
  logical                       :: is_thermal
  logical                       :: is_mat_thermal
  logical                       :: is_mat_thermal_td
  logical                       :: is_tbc
  logical                       :: is_flux_bc
  logical                       :: is_convection_bc
  ! Thermal: End
  
  
  ! Hilfsgrößen, um Aufwand zu sparen
  integer                       :: num_dof              ! number of dof in whole structure
  integer                       :: num_nodes            ! number of nodes
  integer                       :: num_elements         ! number of elements
  integer                       :: num_loadcases        ! number of loadcases
  integer                       :: num_subcases         ! number of subcases
  integer                       :: num_cps              ! number of coulings(-sets)
  
  ! Einstellungen der Berechnung und Ausgabe
  integer                       :: sol
  integer                       :: numEigVal = 1        ! number of Eigenvalues
  integer                       :: shellResPos = 0      ! die Ergebnisse der Shell-Ausgabe soll 0-in der Mittelebene, 1-auf der Oberseite, 2-auf der Unterseite erfolgen
  integer                       :: blocksize            ! blocksize for estimation of nonzero elements in the stiffness matrix
  logical                       :: skipFailed           ! true, wenn im Falle eines Strukturversagens oder bei Beulen, alle nachfolgenden Lastfälle nicht mehr berechnet werden sollen

  logical                       :: calculateTSE         ! true, wenn Gesamtdehnungsenergie (Total Strain Energy) berechnet werden soll

  logical                       :: calculateElementalTSE! true, wenn elementweise Gesamtdehnungsenergie (Total Strain Energy) berechnet werden soll
  
  logical                       :: calculateReactForce  ! true, wenn die Reaktionskräfte in einem zweiten Berechnungsschritt bestimmt werden sollen
  logical                       :: globalReactForce     ! true, wenn die Reaktionskräfte im globalen statt im Knotenkoordinatensystem ausgegeben werden sollen
  
  double precision              :: aeroDisplEps         ! Konvergenzkriterium für die Fluid-Strukturkopplung

end type fe_simulation

contains

    !
    !> $Id: fesimulation_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
    !> $Author: s1080304 $
    !> $Revision: 484 $
    !> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
  subroutine init_default(fesim)
    ! =================================================================================================
    !
    !	Header:		
    !
    !	Content:	
    !
    !	Input:		
    !
    !	Output:		
    !
    !	Internal:	
    !
    !	Calls:		
    !
    !	Called by:	
    !
    !	Author:		
    !
    !	Revision:	
    !
    ! =================================================================================================
    !
    ! Use
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

    ! =================================================================================================
    !
    ! Data types
    !
    type(fe_simulation) :: fesim
    integer             :: err_code=0
    !
    ! =================================================================================================
    !
    ! Initialsiation
    !

    fesim%sol           = 1                  ! solution type; 1-linear static, 2-prelinearized buckling
    fesim%numEigVal     = 1                  ! number of Eigenvalues
    fesim%shellResPos   = 0                  ! Position der Ergebnisse für Schalenelemente 0 - MID, 1 - TOP,
    
    fesim%blocksize = HUGE(fesim%blocksize)  ! set the blocksize to the maximum integer number
    
    fesim%is_beam2                 = .FALSE.
    fesim%is_coord                 = .FALSE.
    fesim%is_coupling              = .FALSE.
    fesim%is_force                 = .FALSE.
    fesim%is_lam8                  = .FALSE.
    fesim%is_lam20                 = .FALSE.
    fesim%is_load                  = .FALSE.
    fesim%is_mat1                  = .FALSE.
    fesim%is_mat8                  = .FALSE.
    fesim%is_mat20                 = .FALSE.
    fesim%is_moment                = .FALSE.
    fesim%is_mpcadd                = .FALSE.
    fesim%is_mpc                   = .FALSE.
    fesim%is_node                  = .FALSE.
    fesim%is_p2load                = .FALSE.
    fesim%is_p8load                = .FALSE.
    fesim%is_p20load               = .FALSE.
    fesim%is_aeroload2d            = .FALSE.
    fesim%is_aeroload3d            = .FALSE.
    fesim%is_temperature           = .FALSE.
    fesim%is_beam2temp             = .FALSE.
    fesim%is_quad8temp             = .FALSE.
    fesim%is_lsolid20temp          = .FALSE.
    fesim%is_pbeam                 = .FALSE.
    fesim%is_pcomp                 = .FALSE.
    fesim%is_pshell                = .FALSE.
    fesim%is_plsolid               = .FALSE.
    fesim%is_quad8                 = .FALSE.
    fesim%is_lsolid20              = .FALSE.
    fesim%is_spc1                  = .FALSE.
    fesim%is_spcadd                = .FALSE.
    fesim%is_spcd                  = .FALSE.
    fesim%is_subcase               = .FALSE.
    fesim%is_contact_node_beam2    = .FALSE.
    fesim%is_contact_node_quad8    = .FALSE.
    fesim%is_contact_node_lsolid20 = .FALSE.
    fesim%is_multistep             = .FALSE.
    ! Thermal: START
    fesim%is_thermal               = .FALSE.
    fesim%is_mat_thermal           = .FALSE.
    fesim%is_mat_thermal_td        = .FALSE.
    fesim%is_tbc                   = .FALSE.
    fesim%is_convection_bc         = .FALSE.
    fesim%is_flux_bc               = .FALSE.
    ! Thermal: END

    fesim%skipFailed               = .FALSE.
    fesim%calculateTSE             = .FALSE.
    fesim%calculateElementalTSE    = .FALSE.
    fesim%calculateReactForce      = .FALSE.
    fesim%globalReactForce         = .FALSE.
    
    fesim%aeroDisplEps  = 0.01d0
    
    fesim%internals%failed = .FALSE.
    fesim%internals%indexBeforeAeroP8Load = 0
    fesim%internals%highestLID = 0
    fesim%internals%indexBeforeAeroLoad = 0
    fesim%internals%highestLIDtemp = 0
    
    fesim%ausgabe%outputVTK = .false.
    fesim%ausgabe%outputUser = .false.
    fesim%ausgabe%outputShort = .false.
    fesim%ausgabe%outputKoopt = .false.
    fesim%ausgabe%outputAdviLa = .false.
    fesim%ausgabe%outputOptitube = .false.
    fesim%ausgabe%outputGlawi = .false.
    fesim%ausgabe%outputHyMoWi = .false.
    fesim%ausgabe%outputElemCoord = .false.
    fesim%ausgabe%outputBoundCond = .false.
    fesim%ausgabe%outputApamePressures = .false.
    !
    ! =================================================================================================
    !
    ! Calculation
    !

    !
    ! =================================================================================================
    !
    ! Error handling
    !
    9999 continue

    if (err_code /= 0) then
    
      write(*,*)                      'An error occured in subroutine'
      write(*,*)                      'init_default'
      write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
      write(*,*)                      'exit program '
      stop
    
    end if

  end subroutine init_default

  subroutine bcast_fesim(fesim)
    
#include "petsc/finclude/petscsys.h"
    use petscsys
    
    implicit none

    type(fe_simulation) :: fesim
    PetscMPIInt     :: rank
    PetscErrorCode  :: ierr
    
#if !defined (PETSC_HAVE_MPIUNI)
    
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
    
    call MPI_Bcast (              fesim%num_dof, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%num_nodes, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (         fesim%num_elements, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (        fesim%num_loadcases, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (         fesim%num_subcases, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (              fesim%num_cps, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (                  fesim%sol, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%numEigVal, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (          fesim%shellResPos, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%blocksize, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
     
    
    call MPI_Bcast (              fesim%is_node, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_beam2, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_quad8, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (          fesim%is_lsolid20, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (              fesim%is_load, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_force, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%is_moment, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%is_p2load, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%is_p8load, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (           fesim%is_p20load, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (        fesim%is_aeroload2d, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (        fesim%is_aeroload3d, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (              fesim%is_mat1, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (              fesim%is_mat8, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_mat20, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_pbeam, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%is_pshell, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_pcomp, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (           fesim%is_plsolid, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%is_spcadd, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (              fesim%is_spc1, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (              fesim%is_spcd, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (            fesim%is_mpcadd, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (               fesim%is_mpc, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (       fesim%is_temperature, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (         fesim%is_beam2temp, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (         fesim%is_quad8temp, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (      fesim%is_lsolid20temp, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (           fesim%is_subcase, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (              fesim%is_lam8, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_lam20, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (          fesim%is_coupling, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (fesim%is_contact_node_beam2, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (fesim%is_contact_node_quad8, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (           fesim%is_failure, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (             fesim%is_coord, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (         fesim%is_multistep, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (           fesim%skipFailed, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (         fesim%calculateTSE, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (fesim%calculateElementalTSE, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (  fesim%calculateReactForce, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (     fesim%globalReactForce, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

    ! Thermal: START
    call MPI_Bcast (           fesim%is_thermal, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (       fesim%is_mat_thermal, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (    fesim%is_mat_thermal_td, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (               fesim%is_tbc, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (     fesim%is_convection_bc, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast (           fesim%is_flux_bc, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    ! Thermal: END
    
    
    call bcast_knoten(fesim%knoten,fesim%is_node)
    call bcast_element_typen(fesim%elemente,fesim%is_beam2,fesim%is_quad8,fesim%is_lsolid20)
    call bcast_eigenschaften(fesim%eigenschaften,fesim%is_pbeam,fesim%is_pshell,fesim%is_pcomp,fesim%is_plsolid)
    !versagenskriterien
    call bcast_ksys(fesim%koordinatensysteme,fesim%is_coord)
    !laminate
    call bcast_lasten(fesim%lasten)
    call bcast_material(fesim%materialien,fesim%is_mat1,fesim%is_mat8,fesim%is_mat20)
    call bcast_rb(fesim%randbedingungen,fesim%is_mpc)
    call bcast_ergebnisse(fesim%ergebnisse)
    !residuals
    !thermal
    call bcast_thermal(fesim%thermal,fesim%is_thermal,fesim%is_mat_thermal,fesim%is_mat_thermal_td,fesim%is_tbc,fesim%is_convection_bc,fesim%is_flux_bc)

#endif
    
  end subroutine bcast_fesim
  
  subroutine free_mem_fesim(fesim)
  
    implicit none 
    
    type(fe_simulation) :: fesim
    
    call free_mem_eigenschaften(fesim%eigenschaften)
    call free_mem_element(fesim%elemente)
    call free_mem_fail(fesim%versagenskriterien)
    call free_mem_knoten(fesim%knoten)
    call free_mem_ksys(fesim%koordinatensysteme)
    call free_mem_laminate(fesim%laminate)
    call free_mem_lasten(fesim%lasten)
    call free_mem_material(fesim%materialien)
    call free_mem_rb(fesim%randbedingungen)
    call free_mem_residual(fesim%residuals)
    call free_mem_ergebnisse(fesim%ergebnisse)

  end subroutine free_mem_fesim

end module fesimulation_typen
