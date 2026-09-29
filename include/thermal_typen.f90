!
! FiPPS2 thermal extension
! Steady thermal conduction and thermo-mechanical coupling
!
module thermal_typen

  implicit none

  type thermal_curve_type
    integer :: n = 0
    double precision, allocatable :: T(:)
    double precision, allocatable :: value(:)
  end type thermal_curve_type

  type thermal_material_type ! Temperaturunabhängiges thermisches Material
    integer :: mid = 0
    double precision :: tref = 0.0d0
    double precision :: k
  end type thermal_material_type

  type thermal_material_td_type ! Temperaturabhängig thermisches Material
    integer :: mid = 0
    double precision :: tref = 0.0d0
    type(thermal_curve_type) :: k
  end type thermal_material_td_type

  type thermal_tbc_type
    integer :: lid = 0  ! LoadID
    integer :: nid = 0  ! NodeID
    double precision :: temperature = 0.d0 ! Temperature
  end type thermal_tbc_type

  type thermal_flux_bc_type
    integer :: lid = 0              ! LoadID
    integer :: eid = 0              ! ID of element
    integer :: face = 0             ! ID of element face
    double precision :: q = 0.d0    ! Heat flux across face
  end type thermal_flux_bc_type

  type thermal_convection_bc_type
    integer :: lid = 0              ! LoadID
    integer :: eid = 0              ! ID of element
    integer :: face = 0             ! ID of element face
    double precision :: h = 0.d0    ! Film coefficient, q = Q/A = h * (T - T_amb)
    double precision :: T_amb = 0.d0    ! Ambient temperature at Face
  end type thermal_convection_bc_type

  type thermal_state_type
    logical :: enabled = .false. ! Is thermal calculation supposed to be performed?
    logical :: coupled_to_structure = .false.   ! Translate calculated temperatures to stress analysis
    logical :: thermal_only = .false.   ! Do only thermal?
    logical :: iterative = .false.             ! Does the model include temperature dependent materials --> if so, iterative solve is necessary
    logical :: successful = .false. ! war die thermische Berechnung erfolgreich
    integer :: max_iterations = 50      ! Max Iterations for iterative (nonLinear)
    double precision :: tolerance = 1.d-8   ! Convergence tolerance
    double precision :: relaxation = 0.8d0  ! Relaxation factor for NL line search
    double precision :: tref = 293.15d0     ! Reference temperature --> Material also has a tref??
    double precision, allocatable :: temperature(:)
    double precision, allocatable :: temperature_old(:)
    type(thermal_material_type), allocatable :: materials(:)
    type(thermal_material_td_type), allocatable :: materials_td(:)
    type(thermal_tbc_type), allocatable :: tbcs(:)
    type(thermal_convection_bc_type), allocatable :: convections(:)
    type(thermal_flux_bc_type), allocatable :: fluxes(:)
  end type thermal_state_type

contains

  subroutine thermal_curve_add(curve, temperature, value)
    type(thermal_curve_type), intent(inout) :: curve
    double precision, intent(in) :: temperature, value
    double precision, allocatable :: ttmp(:), vtmp(:)
    integer :: i, j

    if (.not.allocated(curve%T)) then
      allocate(curve%T(1), curve%value(1))
      curve%T(1) = temperature
      curve%value(1) = value
      curve%n = 1
      return
    end if

    allocate(ttmp(curve%n+1), vtmp(curve%n+1))
    ttmp(1:curve%n) = curve%T
    vtmp(1:curve%n) = curve%value
    ttmp(curve%n+1) = temperature
    vtmp(curve%n+1) = value
    deallocate(curve%T, curve%value)
    allocate(curve%T(curve%n+1), curve%value(curve%n+1))
    curve%T = ttmp
    curve%value = vtmp
    deallocate(ttmp, vtmp)
    curve%n = curve%n + 1

    ! Sorting
    do i = 1, curve%n-1
      do j = i+1, curve%n
        if (curve%T(j) < curve%T(i)) then
          call swap_dp(curve%T(i), curve%T(j))
          call swap_dp(curve%value(i), curve%value(j))
        end if
      end do
    end do
  end subroutine thermal_curve_add

  subroutine swap_dp(a,b)
    double precision, intent(inout) :: a,b
    double precision :: tmp
    tmp = a; a = b; b = tmp
  end subroutine swap_dp

  double precision function thermal_curve_eval(curve, temperature) ! returns the value of a curve for a given temperature
    type(thermal_curve_type), intent(in) :: curve
    double precision, intent(in) :: temperature
    integer :: i

    if (curve%n <= 0) then
      thermal_curve_eval = 0.d0
      return
    end if
    if (curve%n == 1) then
      thermal_curve_eval = curve%value(1)
      return
    end if

    if (temperature <= curve%T(1)) then
      thermal_curve_eval = curve%value(1)
      return
    end if
    if (temperature >= curve%T(curve%n)) then
      thermal_curve_eval = curve%value(curve%n)
      return
    end if

    do i = 1, curve%n-1
      if (temperature <= curve%T(i+1)) then
        thermal_curve_eval = curve%value(i) + &
             (curve%value(i+1)-curve%value(i)) * &
             (temperature-curve%T(i))/(curve%T(i+1)-curve%T(i))
        return
      end if
    end do
    thermal_curve_eval = curve%value(curve%n)
  end function thermal_curve_eval

  double precision function thermal_curve_integral(curve, t0, t1)
    type(thermal_curve_type), intent(in) :: curve
    double precision, intent(in) :: t0,t1
    double precision :: a,b,fa,fb,lo,hi,dt
    integer :: i
    logical :: reverse

    if (curve%n <= 0) then
      thermal_curve_integral = 0.d0
      return
    end if

    reverse = .false.
    lo = t0; hi = t1
    if (hi < lo) then
      reverse = .true.
      lo = t1; hi = t0
    end if

    if (abs(hi-lo) <= tiny(1.d0)) then
      thermal_curve_integral = 0.d0
      return
    end if

    thermal_curve_integral = 0.d0
    a = lo
    fa = thermal_curve_eval(curve,a)

    do i = 1, curve%n-1
      if (curve%T(i+1) <= lo) cycle
      if (curve%T(i) >= hi) exit

      b = min(hi, max(lo,curve%T(i+1)))
      if (b <= a) cycle
      fb = thermal_curve_eval(curve,b)
      dt = b-a
      thermal_curve_integral = thermal_curve_integral + 0.5d0*(fa+fb)*dt
      a = b
      fa = fb
      if (a >= hi) exit
    end do

    if (a < hi) then
      fb = thermal_curve_eval(curve,hi)
      thermal_curve_integral = thermal_curve_integral + 0.5d0*(fa+fb)*(hi-a)
    end if

    if (reverse) thermal_curve_integral = -thermal_curve_integral
  end function thermal_curve_integral

  subroutine thermal_get_material(state, mid, material, found)
    type(thermal_state_type), intent(in) :: state
    integer, intent(in) :: mid
    type(thermal_material_type), intent(out) :: material
    logical, intent(out) :: found
    integer :: i

    material = thermal_material_type()
    found = .false.
    if (.not.allocated(state%materials)) return
    do i=1,size(state%materials)
      if (state%materials(i)%mid == mid) then
        material = state%materials(i)
        found = .true.
        return
      end if
    end do
  end subroutine thermal_get_material

  subroutine thermal_free_state(state)
    type(thermal_state_type), intent(inout) :: state
    integer :: i

    if (allocated(state%temperature)) deallocate(state%temperature)
    if (allocated(state%temperature_old)) deallocate(state%temperature_old)

    if (allocated(state%materials)) then
      do i=1,size(state%materials)
        if (allocated(state%materials(i)%k%T)) deallocate(state%materials(i)%k%T)
        if (allocated(state%materials(i)%k%value)) deallocate(state%materials(i)%k%value)
        if (allocated(state%materials(i)%E%T)) deallocate(state%materials(i)%E%T)
        if (allocated(state%materials(i)%E%value)) deallocate(state%materials(i)%E%value)
        if (allocated(state%materials(i)%nu%T)) deallocate(state%materials(i)%nu%T)
        if (allocated(state%materials(i)%nu%value)) deallocate(state%materials(i)%nu%value)
        if (allocated(state%materials(i)%alpha%T)) deallocate(state%materials(i)%alpha%T)
        if (allocated(state%materials(i)%alpha%value)) deallocate(state%materials(i)%alpha%value)
      end do
      deallocate(state%materials)
    end if

    if (allocated(state%tbcs)) deallocate(state%tbcs)
    if (allocated(state%convection)) deallocate(state%convection)
    if (allocated(state%fluxes)) deallocate(state%fluxes)
    if (allocated(state%heat_sources)) deallocate(state%heat_sources)
  end subroutine thermal_free_state

#include "petsc/finclude/petscsys.h"
  subroutine bcast_thermal(state)
    use petscsys
    implicit none
    type(thermal_state_type), intent(inout) :: state
    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    integer :: i,n,ii

#if !defined (PETSC_HAVE_MPIUNI)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

    call MPI_Bcast(state%enabled,1,MPI_LOGICAL,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    call MPI_Bcast(state%coupled_to_structure,1,MPI_LOGICAL,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    call MPI_Bcast(state%thermal_only,1,MPI_LOGICAL,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    call MPI_Bcast(state%max_iterations,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    call MPI_Bcast(state%tolerance,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    call MPI_Bcast(state%relaxation,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    call MPI_Bcast(state%tref,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)

    if (rank == 0) then
      if (allocated(state%temperature)) then
        n = size(state%temperature)
      else
        n = 0
      end if
    end if
    call MPI_Bcast(n,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    if (rank /= 0) then
      if (allocated(state%temperature)) deallocate(state%temperature)
      if (allocated(state%temperature_old)) deallocate(state%temperature_old)
      if (n > 0) then
        allocate(state%temperature(n),state%temperature_old(n))
      else
        allocate(state%temperature(0),state%temperature_old(0))
      end if
    end if
    if (n > 0 .and. rank == 0 .and. .not.allocated(state%temperature_old)) then
      allocate(state%temperature_old(n))
      state%temperature_old = state%temperature
    end if
    if (n > 0) then
      call MPI_Bcast(state%temperature,n,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
      call MPI_Bcast(state%temperature_old,n,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
    end if
#endif
  end subroutine bcast_thermal

  function thermal_is_td(state) result(is_td)
    type(thermal_state_type), intent(in) :: state
    logical :: is_td
    
    is_td = (allocated(state%materials_td) .and. size(state%materials_td) > 0)
    
  end function thermal_is_td

  function thermal_is_conv(state) result(is_conv)
    type(thermal_state_type), intent(in) :: state
    logical :: is_conv
    
    is_conv = (allocated(state%convections) .and. size(state%convections) > 0)
    
  end function thermal_is_conv
  

  end module thermal_typen
