
! =================================================================================================
!
!> @brief Write the finite-element mesh in VTK XML UnstructuredGrid format.
!
!> @details
!> Supported element types:
!>   beam2    -> VTK_LINE                  (3)
!>   quad8    -> VTK_QUADRATIC_QUAD        (23)
!>   lsolid20 -> VTK_QUADRATIC_HEXAHEDRON  (25)
!
!> The element ordering is beam2, quad8, lsolid20, matching the legacy
!> vtkoutmesh routine.
!
! =================================================================================================
SUBROUTINE vtkoutmesh_new(fesim, fname)

  USE fesimulation_typen
  USE, INTRINSIC :: iso_fortran_env, ONLY: int64

  IMPLICIT NONE

  TYPE(fe_simulation), INTENT(IN) :: fesim
  CHARACTER(*), INTENT(IN)        :: fname

  INTEGER :: unit_vtk
  INTEGER :: ios
  INTEGER :: ii, jj
  INTEGER :: num_points
  INTEGER :: num_cells
  INTEGER :: num_beam2, num_quad8, num_lsolid20
  INTEGER(int64) :: num_connectivity
  INTEGER(int64) :: offset

  ! -----------------------------------------------------------------------------------------------
  ! Determine mesh dimensions
  ! -----------------------------------------------------------------------------------------------

  num_points = SIZE(fesim%knoten%nodes)

  num_beam2   = 0
  num_quad8   = 0
  num_lsolid20 = 0

  IF (ALLOCATED(fesim%elemente%beam2s)) THEN
    num_beam2 = SIZE(fesim%elemente%beam2s)
  END IF

  IF (ALLOCATED(fesim%elemente%quad8s)) THEN
    num_quad8 = SIZE(fesim%elemente%quad8s)
  END IF

  IF (ALLOCATED(fesim%elemente%lsolid20s)) THEN
    num_lsolid20 = SIZE(fesim%elemente%lsolid20s)
  END IF

  num_cells = num_beam2 + num_quad8 + num_lsolid20

  num_connectivity = 2_int64  * INT(num_beam2, int64)   &
                   + 8_int64  * INT(num_quad8, int64)   &
                   + 20_int64 * INT(num_lsolid20, int64)

  ! -----------------------------------------------------------------------------------------------
  ! Validate element connectivity before creating the output file
  ! -----------------------------------------------------------------------------------------------

  IF (ALLOCATED(fesim%elemente%beam2s)) THEN
    DO ii = 1, num_beam2
      CALL check_node_ids(fesim%elemente%beam2s(ii)%nids, num_points, &
                          fesim%elemente%beam2s(ii)%eid)
    END DO
  END IF

  IF (ALLOCATED(fesim%elemente%quad8s)) THEN
    DO ii = 1, num_quad8
      CALL check_node_ids(fesim%elemente%quad8s(ii)%nids, num_points, &
                          fesim%elemente%quad8s(ii)%eid)
    END DO
  END IF

  IF (ALLOCATED(fesim%elemente%lsolid20s)) THEN
    DO ii = 1, num_lsolid20
      CALL check_node_ids(fesim%elemente%lsolid20s(ii)%nids, num_points, &
                          fesim%elemente%lsolid20s(ii)%eid)
    END DO
  END IF

  ! -----------------------------------------------------------------------------------------------
  ! Open output file
  ! -----------------------------------------------------------------------------------------------

  OPEN(NEWUNIT=unit_vtk, FILE=TRIM(fname), STATUS='REPLACE', &
       ACTION='WRITE', FORM='FORMATTED', IOSTAT=ios)

  IF (ios /= 0) THEN
    WRITE(*,'(A)') 'ERROR: Could not open VTK XML output file: '//TRIM(fname)
    ERROR STOP
  END IF

  ! -----------------------------------------------------------------------------------------------
  ! XML header
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '<?xml version="1.0"?>'
  WRITE(unit_vtk,'(A)') '<VTKFile type="UnstructuredGrid" version="1.0" '// &
                        'byte_order="LittleEndian">'
  WRITE(unit_vtk,'(A)') '  <UnstructuredGrid>'
  WRITE(unit_vtk,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', num_points, &
                                  '" NumberOfCells="', num_cells, '">'

  ! -----------------------------------------------------------------------------------------------
  ! Points
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '      <Points>'
  WRITE(unit_vtk,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" '// &
                        'format="ascii">'

  DO ii = 1, num_points
    WRITE(unit_vtk,'(3(ES25.16E3,1X))') &
      fesim%knoten%nodes(ii)%coords(1:3)
  END DO

  WRITE(unit_vtk,'(A)') '        </DataArray>'
  WRITE(unit_vtk,'(A)') '      </Points>'

  ! -----------------------------------------------------------------------------------------------
  ! Cells: connectivity
  !
  ! VTK uses zero-based point indices. The FiPPS node IDs are one-based.
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '      <Cells>'

  WRITE(unit_vtk,'(A)') '        <DataArray type="Int64" Name="connectivity" '// &
                        'format="ascii">'

  ! beam2
  IF (ALLOCATED(fesim%elemente%beam2s)) THEN
    DO ii = 1, num_beam2
      DO jj = 1, 2
        WRITE(unit_vtk,'(I0,1X)',ADVANCE='NO') &
          INT(fesim%elemente%beam2s(ii)%nids(jj)-1, int64)
      END DO
      WRITE(unit_vtk,*)
    END DO
  END IF

  ! quad8
  IF (ALLOCATED(fesim%elemente%quad8s)) THEN
    DO ii = 1, num_quad8
      DO jj = 1, 8
        WRITE(unit_vtk,'(I0,1X)',ADVANCE='NO') &
          INT(fesim%elemente%quad8s(ii)%nids(jj)-1, int64)
      END DO
      WRITE(unit_vtk,*)
    END DO
  END IF

  ! lsolid20
  IF (ALLOCATED(fesim%elemente%lsolid20s)) THEN
    DO ii = 1, num_lsolid20
      DO jj = 1, 20
        WRITE(unit_vtk,'(I0,1X)',ADVANCE='NO') &
          INT(fesim%elemente%lsolid20s(ii)%nids(jj)-1, int64)
      END DO
      WRITE(unit_vtk,*)
    END DO
  END IF

  WRITE(unit_vtk,'(A)') '        </DataArray>'

  ! -----------------------------------------------------------------------------------------------
  ! Cells: offsets
  !
  ! Each offset is the cumulative number of connectivity entries.
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '        <DataArray type="Int64" Name="offsets" '// &
                        'format="ascii">'

  offset = 0_int64

  DO ii = 1, num_beam2
    offset = offset + 2_int64
    WRITE(unit_vtk,'(I0)') offset
  END DO

  DO ii = 1, num_quad8
    offset = offset + 8_int64
    WRITE(unit_vtk,'(I0)') offset
  END DO

  DO ii = 1, num_lsolid20
    offset = offset + 20_int64
    WRITE(unit_vtk,'(I0)') offset
  END DO

  WRITE(unit_vtk,'(A)') '        </DataArray>'

  ! -----------------------------------------------------------------------------------------------
  ! Cells: VTK cell types
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '        <DataArray type="UInt8" Name="types" '// &
                        'format="ascii">'

  DO ii = 1, num_beam2
    WRITE(unit_vtk,'(I0)') 3
  END DO

  DO ii = 1, num_quad8
    WRITE(unit_vtk,'(I0)') 23
  END DO

  DO ii = 1, num_lsolid20
    WRITE(unit_vtk,'(I0)') 25
  END DO

  WRITE(unit_vtk,'(A)') '        </DataArray>'
  WRITE(unit_vtk,'(A)') '      </Cells>'

  ! -----------------------------------------------------------------------------------------------
  ! Point data
  !
  ! No nodal fields are written here. They can be added later as DataArrays.
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '      <PointData/>'

  ! -----------------------------------------------------------------------------------------------
  ! Cell data
  !
  ! These arrays follow exactly the same element ordering as connectivity.
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '      <CellData Scalars="element_id">'

  ! Element IDs
  WRITE(unit_vtk,'(A)') '        <DataArray type="Int32" Name="element_id" '// &
                        'format="ascii">'

  IF (ALLOCATED(fesim%elemente%beam2s)) THEN
    DO ii = 1, num_beam2
      WRITE(unit_vtk,'(I0)') fesim%elemente%beam2s(ii)%eid
    END DO
  END IF

  IF (ALLOCATED(fesim%elemente%quad8s)) THEN
    DO ii = 1, num_quad8
      WRITE(unit_vtk,'(I0)') fesim%elemente%quad8s(ii)%eid
    END DO
  END IF

  IF (ALLOCATED(fesim%elemente%lsolid20s)) THEN
    DO ii = 1, num_lsolid20
      WRITE(unit_vtk,'(I0)') fesim%elemente%lsolid20s(ii)%eid
    END DO
  END IF

  WRITE(unit_vtk,'(A)') '        </DataArray>'

  ! Property IDs
  WRITE(unit_vtk,'(A)') '        <DataArray type="Int32" Name="property_id" '// &
                        'format="ascii">'

  IF (ALLOCATED(fesim%elemente%beam2s)) THEN
    DO ii = 1, num_beam2
      WRITE(unit_vtk,'(I0)') fesim%elemente%beam2s(ii)%pid
    END DO
  END IF

  IF (ALLOCATED(fesim%elemente%quad8s)) THEN
    DO ii = 1, num_quad8
      WRITE(unit_vtk,'(I0)') fesim%elemente%quad8s(ii)%pid
    END DO
  END IF

  IF (ALLOCATED(fesim%elemente%lsolid20s)) THEN
    DO ii = 1, num_lsolid20
      WRITE(unit_vtk,'(I0)') fesim%elemente%lsolid20s(ii)%pid
    END DO
  END IF

  WRITE(unit_vtk,'(A)') '        </DataArray>'
  WRITE(unit_vtk,'(A)') '      </CellData>'

  ! -----------------------------------------------------------------------------------------------
  ! Close XML structure
  ! -----------------------------------------------------------------------------------------------

  WRITE(unit_vtk,'(A)') '    </Piece>'
  WRITE(unit_vtk,'(A)') '  </UnstructuredGrid>'
  WRITE(unit_vtk,'(A)') '</VTKFile>'

  CLOSE(unit_vtk)

  WRITE(*,'(A)') 'VTK XML mesh written: '//TRIM(fname)

CONTAINS

  ! -----------------------------------------------------------------------------------------------
  ! Check that all element node IDs refer to existing mesh nodes.
  ! -----------------------------------------------------------------------------------------------

  SUBROUTINE check_node_ids(nids, number_of_nodes, eid)

    INTEGER, INTENT(IN) :: nids(:)
    INTEGER, INTENT(IN) :: number_of_nodes
    INTEGER, INTENT(IN) :: eid

    INTEGER :: kk

    DO kk = 1, SIZE(nids)
      IF (nids(kk) < 1 .OR. nids(kk) > number_of_nodes) THEN
        WRITE(*,'(A,I0,A,I0,A,I0)') &
          'ERROR: Invalid node ID ', nids(kk), &
          ' in element ', eid, &
          '. Number of mesh nodes: ', number_of_nodes
        ERROR STOP
      END IF
    END DO

  END SUBROUTINE check_node_ids

END SUBROUTINE vtkoutmesh_new
