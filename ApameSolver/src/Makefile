# Makefile
# erzeugt das Programm apame

include ../../makeconfig.mk

BIN_PATH := ../bin

PROGRAM := apame

LIBRARY := apame.a

MAIN := \
	apame.o

OBJECTS :=  \
	module_kind_and_konst.o \
	module_grid.o \
	module_induced_drag.o \
	module_influence.o \
	module_input.o \
	module_modify_wake.o \
	module_pressure.o \
	module_rhs.o \
	module_velocity.o \
	func_IDW.o \
	func_message.o \
	func_new_line.o \
	func_tic.o \
	func_title.o \
	func_toc.o \
	func_vect_ang.o \
	func_vect_inprod.o \
	func_vect_mod.o \
	func_vect_prod.o \
	subr_apame_main.o \
	subr_vtk_results.o \
	subr_results.o \
	subr_speeds.o \
	subr_write_case.o \
	subr_write_coef.o \
	subr_write_key.o \
	subr_write_ptot.o \
	subr_write_time.o

###############################################################################
## rules follow
###############################################################################
# .PHONY: default

default all: $(PROGRAM) move

$(PROGRAM): $(LIBRARY) $(MAIN)
	$(FC) -o $(PROGRAM) $(MAIN) $(BLASLAPACK_LIB) $(LIBRARY) $(LIB)

$(LIBRARY): $(OBJECTS)
	$(AR) $(AR_FLAGS) $(LIBRARY) $^

move:
	mkdir -p $(BIN_PATH)
	$(MV) $(PROGRAM) $(BIN_PATH)

clean::
	rm -rf $(LIBRARY) $(BIN_PATH)
