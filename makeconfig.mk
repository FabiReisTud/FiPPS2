#  This program developed in FORTRAN is a Finite Element solver for linear-
#  static analyses as well as linearized stability analyses. It is inherently
#  coupled to the open-source panel method APAME for providing fluid-structure-
#  interaction capabilites.
#
#  Copyright (C) 2024 TUD Dresden University of Technology
#
#  This file is part of FiPPS².
#
#  FiPPS² is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FiPPS² is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#===============================================================================
# makeconfig.mk:  common configuration file for FiPPS
#===============================================================================

# fast - without debugging options
PETSC_DIR   = /mnt/appl/x86_64/petsc/3.16.1
SLEPC_DIR   = /mnt/appl/x86_64/slepc/3.16.0_mod
APAME_LIB   = ../ApameSolver/src/apame.a
PANEL2D_LIB = ../PANEL2D/panel2d.a
XFOILWRAPPER_LIB = ../XfoilWrapper/xfoilwrapper.a
FIPPS2_LIB  = libfipps.a

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

VERSION=2.0.0
FIPPS_VERSION=2.0.0
REVISION=6969

# Fortran compiler
# LFT
FC  = mpif90
# # HPC Barnard
# FC  = mpiifort

# Standardoptionen
#FC_FLAGS =
#FC_FLAGS = -g -traceback -check all -warn unused -warn uninitialized -warn interfaces -check noarg_temp_created

# Standardoptionen für die modifizierte SLEPC-Version
FC_FLAGS = -Dmod_slepc
#FC_FLAGS = -Dmod_slepc -Dvideo
#FC_FLAGS = -Dmod_slepc -g -traceback -check all -warn unused -warn uninitialized -warn interfaces -check noarg_temp_created

CFLAGS =
FFLAGS = -fpp -fPIC -I../include/ -I../../include ${SLEPC_INCLUDE} ${PETSC_FC_INCLUDES}
CPPFLAGS =
FPPFLAGS =
