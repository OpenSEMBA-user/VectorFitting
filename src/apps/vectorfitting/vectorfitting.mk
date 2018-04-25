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
# =============================================================================
OUT = vectorfitting
static = yes

# -------------------- Paths to directories -----------------------------------
SRC_CORE_DIRS     := $(shell find $(SRC_DIR)core/ -type d)

SRC_DIRS = $(SRC_CORE_DIRS)

SRCS_CXX := $(shell find $(SRC_DIRS) -maxdepth 1 -type f -name "*.cpp")	
OBJS_CXX := $(addprefix $(OBJ_DIR), $(SRCS_CXX:.cpp=.o))
# =============================================================================
INCLUDES += $(SRC_DIR)core/ 
# =============================================================================
.PHONY: default print

default: $(OUT) 
	@echo "======================================================="
	@echo "           $(OUT) compilation finished                 "
	@echo "======================================================="

$(OBJ_DIR)%.o: %.cpp
	@dirname $@ | xargs mkdir -p
	@echo "Compiling:" $@
	$(CXX) $(CXXFLAGS) $(addprefix -D, $(DEFINES)) $(addprefix -I,$(INCLUDES)) -c -o $@ $<
	
$(LIB_DIR)/$(OUT)/lib/lib$(OUT).a: $(OBJS_CXX)
	@mkdir -p $(LIB_DIR)$(OUT)/lib/ 
	@mkdir -p $(LIB_DIR)$(OUT)/include/
	@echo "Linking:" $@
	-ar rs $@ $^
	-cd $(SRC_DIR); find core/ \( -name "*.h" -o -name "*.hpp" \) -exec cp --parents {} ../$(LIB_DIR)$(OUT)/include/ \;

$(OUT): $(LIB_DIR)/$(OUT)/lib/lib$(OUT).a

print: 
	@echo "======================================================="
	@echo "         ----- Compiling $(OUT) ------                 "
	@echo "Target:           " $(target)
	@echo "Compiler:         " $(compiler)
	@echo "C++ Compiler:     " `which $(CXX)`
	@echo "C++ Flags:        " $(CXXFLAGS)
	@echo "Defines:          " $(DEFINES)
	@echo "======================================================="

# ------------------------------- END ----------------------------------------
