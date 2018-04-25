# OpenSEMBA
# Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
#                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
#                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
#                    Daniel Mateos Romero            (damarro@semba.guru)
#
# This file is part of OpenSEMBA.
#
# OpenSEMBA is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.

# -- USAGE --------------------------------------------------------------------
# make target     = {debug, release}
#      compiler = {intel, gnu, ...}
# ==================== Default values =========================================
target   = release
compiler = gnu

#===================== GNU Compiler ===========================================
ifeq ($(compiler),gnu)
	CC = gcc
	CXX = g++
	CCFLAGS +=
	CXXFLAGS += -std=c++0x -pthread -fopenmp
endif # endif choosing the GNU compiler.
# ================= Optimization target =======================================
ifeq ($(target),debug)
	CXXFLAGS +=-O0 -g3 -Wall -Wno-write-strings

endif
ifeq ($(target),release)
   	CXXFLAGS +=-O2
endif
ifeq ($(target),optimal)
   	CXXFLAGS +=-O3
endif

# =============================================================================
# -------------------- Paths to directories -----------------------------------
BUILD_DIR = ./build/
OBJ_DIR = ./obj/
SRC_DIR = ./src/
EXTERNAL_DIR = ./external/

BIN_DIR = $(BUILD_DIR)bin/
LIB_DIR = $(BUILD_DIR)lib/

# =============================================================================
.NOTPARALLEL:
# -------------------- RULES --------------------------------------------------
default: all
	@echo "======>>>>> Done <<<<<======"

all: check test vectorfitting

test: check
	$(MAKE) -f ./src/apps/test/test.mk print
	$(MAKE) -f ./src/apps/test/test.mk
#	cp -r testData $(BIN_DIR)test/

vectorfitting: check
	$(MAKE) -f ./src/apps/vectorfitting/vectorfitting.mk print
	$(MAKE) -f ./src/apps/vectorfitting/vectorfitting.mk

clean:
	rm -rf $(OBJ_DIR)

clobber: clean
	rm -rf $(BUILD_DIR)

check:
ifneq ($(target),release)
ifneq ($(target),debug)
ifneq ($(target),optimal)
	@echo "Invalid build target."
	@echo "Please use target=[release|debug|optimal]"
	@exit 1
endif
endif
endif
ifneq ($(compiler),gnu)
	@echo "Invalid build compiler"
	@echo "Please use 'make compiler= intel|gnu|mingw32|mingw64'"
	@exit 2
endif

# Exports current variables when other makefiles are called.
export
