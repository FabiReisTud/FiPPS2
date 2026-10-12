# Makefile
# erzeugt die Library für alle Module

include ../makeconfig.mk

PROGRAM := libmodule

MODULE := \
 globale_variablen.o \
 eigenschaften_typen.o \
 element_typen.o \
 failurecriteria_typen.o \
 knoten_typen.o \
 koordinatensystem_typen.o \
 laminat_typen.o \
 last_typen.o \
 material_typen.o \
 randbedingung_typen.o \
 residual_typen.o \
 internals_typen.o \
 ausgabe_typen.o \
 ergebnis_typen.o \
 thermal_typen.o \
 fesimulation_typen.o \
 pre_assemble_types.o \
 vtk_variablen.o \
 vec_func.o \
 dof_types.o \
 konstanten.o \
 mat_func.o \
 integration_schemes.o \
 failure_criteria.o

###############################################################################
## rules follow
###############################################################################
.PHONY: default clean

all: $(PROGRAM)

$(PROGRAM): $(MODULE)
	ar rcs $(PROGRAM).a $(MODULE) 

clean2: clean
	$(RM) $(PROGRAM).a
