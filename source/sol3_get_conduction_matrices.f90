!> @brief
!> control subroutine for obtaining the global conduction matrix from
!> the different elements and their conduction matrices

subroutine sol3_get_conduction_matrices(fesim,CaaS)
    !
  use fesimulation_typen
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscmat.h"
  use petscmat
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
  type(fe_simulation), intent(in)               :: fesim
!
! Inner
!
  integer                                       :: err_code=0
!
! Output
!
  Mat                                           :: CaaS         ! Konduktionsmatrix ohne gesperrte Freiheiten als Sparse-Matrix

  PetscMPIInt                                   :: rank
  PetscErrorCode                                :: ierr
!
! =================================================================================================
!
! Initialisation
!
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
!
! =================================================================================================

! =================================================================================================
!
! Calculation
!
  if (fesim%is_lsolid20 == .true. .AND. fesim%is_thermal == .true.) then
    call tsolid20_cond_control(fesim, CaaS)
  end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol3_get_cond_matrices'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop

end if

end subroutine
