
include ../makeconfig.mk

PROG = PANEL2D

LIBRARY  = panel2d.a

SRC_PATH = ./src/
OBJ_PATH = ./obj/
INC_PATH = -I ./include
MOD_PATH = ./mod/

VPATH = $(SRC_PATH) $(OBJ_PATH) $(INC_PATH)

MAINS = panel2d.o

OBJS  = konst_var.o \
        netz_variablen.o \
        functions.o \
        bl_functions.o \
        runge.o \
        ppval.o \
        spline.o \
        drag_squire_young.o \
        solve_bl_moran.o \
        stagnation_point.o \
        write_vtk.o \
        boundary_layer.o \
        pressure_calculation.o \
        calculate_ABmat.o \
        panel_geometry.o \
        calculate_geometry.o \
        sort_nodes.o \
        add_displacements.o \
        unify_airfoil.o \
        check_airfoil.o \
        read_input.o \
        splineInterpol.o \
        panel2d_main.o \
        xfoil_deps.o \
        refine_mesh.o

OBJ = $(patsubst %,$(OBJ_PATH)%,$(OBJS))

MAIN = $(patsubst %,$(OBJ_PATH)%,$(MAINS))

.PHONY : mkobj

mkobj:
	mkdir -p $(OBJ_PATH)

all: mkobj $(PROG)

$(OBJ_PATH)xfoil_deps.o: $(SRC_PATH)xfoil_deps.f
	$(FC) -c $(FC_FLAGS) -r8 $(SRC_PATH)xfoil_deps.f -o $@

$(OBJ_PATH)%.o: $(SRC_PATH)%.f90
	$(FC) $(FC_FLAGS) -c $< $(INC_PATH) -o $@

$(PROG): $(LIBRARY) $(MAIN)
	$(FC) -o $@ $(MAIN) $(LIBRARY) $(BLASLAPACK_LIB) $(LIBS)

$(LIBRARY): $(OBJ)
	$(AR) $(AR_FLAGS) $(LIBRARY) $^

clean::
	$(RM) -rfv $(PROG).exe $(PROG) $(LIBRARY) $(OBJ_PATH) $(OBJ) $(MAIN) $(MOD_PATH)*.mod
