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
subroutine get_inverse3(Kmat, numDof, invKmat, check_zeros)
  
#include "petsc/finclude/petscksp.h"
  use petscksp

  implicit none
  
  Mat, intent(in)                               :: Kmat
  integer, intent(in)                           :: numDof
  logical, intent(in)                           :: check_zeros
  
  Mat, intent(out)                              :: invKmat
  
  Mat                                           :: tempKmat, KmatFak, Ident, invKmatDense, Kmat_ones
  IS                                            :: perm
  MatFactorInfo                                 :: matfacinfo(MAT_FACTORINFO_SIZE)
  PetscScalar                                   :: petScalar
  integer, dimension(:), allocatable            :: d_nzz 
  integer                                       :: mm, nn, ncols
  integer, dimension(:), allocatable            :: cols
  double precision, dimension(:), allocatable   :: values
  PetscErrorCode                                :: ierr
  
  Vec                                           :: invMatAsVec
  PetscScalar, pointer                          :: xx_v(:)
  
  KSP                                           :: ksp
  PC                                            :: prec
  
  Integer                                       :: ii, jj
  PetscInt                                      :: pii
  PetscScalar                                   :: v
  
!   integer, parameter                            :: fileunit = 83
!   integer                                       :: entrycounter
!   double precision                              :: val
  
  if (check_zeros .eq. .true.) then
  
    ! Erzeugen K einer Kopie von K
    WRITE(*,*) 'Erzeugen einer Kopie von K'
    
    ! Bestimmen der Anzahl der Nichtnullelement in den Zeilen der Matrix
    allocate(d_nzz(numDof), cols(numDof), values(numDof))
    d_nzz = 0
    do mm = 1, numDof
      call MatGetRow(Kmat, mm-1, ncols, cols, values, ierr); CHKERRQ(ierr)
      do nn = 1, ncols
        if (values(nn) .NE. 0.d0) d_nzz(mm) = d_nzz(mm) + 1
      end do
      call MatRestoreRow(Kmat, mm-1, ncols, cols, values, ierr); CHKERRQ(ierr)
    end do

    ! Anlegen der Matrix
    WRITE(*,*) 'Anlegen der Matrix'
    
    call MatCreate(PETSC_COMM_WORLD,tempKmat,ierr); CHKERRQ(ierr)
    call MatSetSizes(tempKmat,PETSC_DECIDE,PETSC_DECIDE,numDof,numDof,ierr); CHKERRQ(ierr)
    call MatSetType(tempKmat,MATMPIAIJ,ierr); CHKERRQ(ierr)
    call MatSetFromOptions(tempKmat,ierr); CHKERRQ(ierr)
    call MatSetUp(tempKmat,ierr); CHKERRQ(ierr)
    call MatMPIAIJSetPreallocation(tempKmat, 0, d_nzz, 0, PETSC_NULL_INTEGER_ARRAY, ierr); CHKERRQ(ierr)
   
    d_nzz = 1
     
    call MatCreate(PETSC_COMM_WORLD,Kmat_ones,ierr); CHKERRQ(ierr)
    call MatSetSizes(Kmat_ones,PETSC_DECIDE,PETSC_DECIDE,numDof,numDof,ierr); CHKERRQ(ierr)
    call MatSetType(Kmat_ones,MATMPIAIJ,ierr); CHKERRQ(ierr)
    call MatSetFromOptions(Kmat_ones,ierr); CHKERRQ(ierr)
    call MatSetUp(Kmat_ones,ierr); CHKERRQ(ierr)
    call MatMPIAIJSetPreallocation(Kmat_ones, 0, d_nzz, 0, PETSC_NULL_INTEGER_ARRAY, ierr); CHKERRQ(ierr)
   
   
    ! Setzen der Werte
    WRITE(*,*) 'Setzen der Werte'
   
    do mm = 1, numDof
      call MatGetRow(Kmat, mm-1, ncols, cols, values, ierr); CHKERRQ(ierr)
      if (ncols .EQ. 0) then
        call MatSetValue(Kmat_ones, mm-1, mm-1, 1.d0, ADD_VALUES, ierr); CHKERRQ(ierr)
      end if
      do nn = 1, ncols
        if (values(nn) .NE. 0.d0) then
          call MatSetValue(tempKmat, mm-1, cols(nn), values(nn), ADD_VALUES, ierr); CHKERRQ(ierr)
        end if
      end do
      call MatRestoreRow(Kmat, mm-1, ncols, cols, values, ierr); CHKERRQ(ierr)
    end do
    
    deallocate(d_nzz, cols, values)
   
    ! Assemblieren
    WRITE(*,*) 'Assemblieren'
    
    call MatAssemblyBegin(tempKmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call MatAssemblyEnd  (tempKmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    
    call MatAssemblyBegin(Kmat_ones,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call MatAssemblyEnd  (Kmat_ones,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    
    call MatAXPY(tempKmat, 1.d0, Kmat_ones, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr)
    
  else 
  
    tempKmat = Kmat
  
  end if

  ! Berechnen von inv(Kii)
  WRITE(*,*) 'Berechnen von inv(K)'
  
  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
  call KSPSetOperators(ksp,tempKmat,tempKmat,ierr); CHKERRQ(ierr)
  
  call KSPSetType(ksp,KSPPREONLY,ierr); CHKERRQ(ierr)
  call KSPGetPC(ksp,prec,ierr); CHKERRQ(ierr)
  call PCSetType(prec,PCCHOLESKY,ierr); CHKERRQ(ierr)
  call PCFactorSetMatSolverType(prec,MATSOLVERMUMPS,ierr); CHKERRQ(ierr)

  call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
  
  allocate(d_nzz(numDof))
  
  d_nzz = 0
  
!   open(unit=fileunit, file='temp_inverse_entries.dat', status='unknown', action='write', form="unformatted")
!   entrycounter = 0
  
  do ii = 1, numDof
  
    call VecCreate(PETSC_COMM_WORLD,invMatAsVec,ierr); CHKERRQ(ierr)
    call VecSetSizes(invMatAsVec,PETSC_DECIDE,numDof,ierr); CHKERRQ(ierr)
    call VecSetFromOptions(invMatAsVec,ierr); CHKERRQ(ierr)
    
    v = 1.d0
    pii = ii-1
    call VecSetValue(invMatAsVec, pii, v, INSERT_VALUES, ierr); CHKERRQ(ierr)
    
    call VecAssemblyBegin(invMatAsVec,ierr); CHKERRQ(ierr);
    call VecAssemblyEnd(invMatAsVec,ierr); CHKERRQ(ierr);
      
    call KSPSolve(ksp,invMatAsVec,invMatAsVec,ierr); CHKERRQ(ierr)
   
    call VecGetArrayF90(invMatAsVec,xx_v,ierr); CHKERRQ(ierr)

    do jj = 1, numDof

      if (xx_v(jj) .NE. 0.d0) then
      
        d_nzz(jj) = d_nzz(jj) + 1
        
!         write(fileunit) ii, jj, xx_v(jj)
!         
!         entrycounter = entrycounter + 1
        
      end if

    end do

    call VecRestoreArrayF90(invMatAsVec,xx_v,ierr); CHKERRQ(ierr)
    call VecDestroy(invMatAsVec,ierr); CHKERRQ(ierr)

    IF (MODULO(ii, 100) .EQ. 0) WRITE(*,*) 'DOF ', ii, ' von ', numDof, ' fertig'

  end do
  
!   close(fileunit)

  call MatCreate(PETSC_COMM_WORLD,invKmat,ierr); CHKERRQ(ierr)
  call MatSetSizes(invKmat,PETSC_DECIDE,PETSC_DECIDE,numDof,numDof,ierr); CHKERRQ(ierr)
  call MatSetType(invKmat,MATMPIAIJ,ierr); CHKERRQ(ierr)
  call MatSetFromOptions(invKmat,ierr); CHKERRQ(ierr)
  call MatSetUp(invKmat,ierr); CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(invKmat, 0, d_nzz, 0, PETSC_NULL_INTEGER_ARRAY, ierr); CHKERRQ(ierr)
  
  WRITE(*,*) 'Inverse bilden'
  WRITE(*,*) 'numDof ', numDof
  
!   open(unit=fileunit, file='temp_inverse_entries.dat', status='unknown', action='read', form="unformatted")
!   
!   do mm = 1, entrycounter
!   
!     read(fileunit) ii, jj, val
!     
!     call MatSetValue(invKmat, jj-1, ii-1, val, ADD_VALUES, ierr); CHKERRQ(ierr)
!   
!   end do
!   
!   close(fileunit)
  
  do ii = 1, numDof
  
    call VecCreate(PETSC_COMM_WORLD,invMatAsVec,ierr); CHKERRQ(ierr)
    call VecSetSizes(invMatAsVec,PETSC_DECIDE,numDof,ierr); CHKERRQ(ierr)
    call VecSetFromOptions(invMatAsVec,ierr); CHKERRQ(ierr)
    
    v = 1.d0
    pii = ii-1
    call VecSetValue(invMatAsVec, pii, v, INSERT_VALUES, ierr); CHKERRQ(ierr)
    
    call VecAssemblyBegin(invMatAsVec,ierr); CHKERRQ(ierr);
    call VecAssemblyEnd(invMatAsVec,ierr); CHKERRQ(ierr);
      
    call KSPSolve(ksp,invMatAsVec,invMatAsVec,ierr); CHKERRQ(ierr)
   
    call VecGetArrayF90(invMatAsVec,xx_v,ierr); CHKERRQ(ierr)

    do jj = 1, numDof

      if (xx_v(jj) .NE. 0.d0) call MatSetValue(invKmat, jj-1, ii-1, xx_v(jj), ADD_VALUES, ierr); CHKERRQ(ierr)

    end do

    call VecRestoreArrayF90(invMatAsVec,xx_v,ierr); CHKERRQ(ierr)
    call VecDestroy(invMatAsVec,ierr); CHKERRQ(ierr)

    IF (MODULO(ii, 100) .EQ. 0) WRITE(*,*) 'DOF ', ii, ' von ', numDof, ' fertig'

  end do
  
  call KSPDestroy(ksp,ierr); CHKERRQ(ierr)

  ! Assemblieren
  WRITE(*,*) 'Assemblieren'
  
  call MatAssemblyBegin(invKmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd  (invKmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  
  if (check_zeros .eq. .true.) then
  
    call MatAXPY(invKmat, -1.d0, Kmat_ones, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr)
    call MatDestroy(Kmat_ones, ierr); CHKERRQ(ierr)

  end if
  
end subroutine get_inverse3
