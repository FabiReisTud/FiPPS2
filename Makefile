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

#-------------------------------------------------------------------------------
# Makefile for FIPPS
#-------------------------------------------------------------------------------

ALL: all

include makeconfig.mk

# Compile the default rules for each package

all:
	( cd ApameSolver/src ; $(MAKE) all)
	( cd PANEL2D ; $(MAKE) all)
	( cd XfoilWrapper ; $(MAKE) all)
	( echo -n "CHARACTER(len=10), PARAMETER :: version=\"$(VERSION)-R$(REVISION)" > $(FIPPS_VERSION) )
	( svn info | grep -i revision | cut -f2 -d: | tr -d [:space:] >> $(FIPPS_VERSION) )
	( echo -n "\"" >> $(FIPPS_VERSION) ) 
	( cd include ; $(MAKE) all )
	( cd source ; $(MAKE) flib )
	( cd outputLibrary ; $(MAKE) all )
	( cd source ; $(MAKE) all )

# flib:
# 	( cd ApameSolver/src ; $(MAKE) all)
# 	( cd PANEL2D ; $(MAKE) all)
# 	( cd XfoilWrapper ; $(MAKE) all)
# 	( cd include ; $(MAKE) all )
# 	( cd source ; $(MAKE) flib)
# 	ar rcs libfipps.a include/*.o source/*.o

# Remove all files
clean::
	( cd ApameSolver/src ; $(MAKE) clean)
	( cd PANEL2D ; $(MAKE) clean)
	( cd XfoilWrapper ; $(MAKE) clean)
	( cd include ; $(MAKE) clean2 )
	( cd source ; $(MAKE) clean2 )
	( cd outputLibrary ; $(MAKE) clean2 )
	rm -rf $(FIPPS) $(FIPPS_VERSION)
