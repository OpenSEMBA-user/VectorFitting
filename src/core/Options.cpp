// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//                    Alejandro García Montoro        (alejandro.garciamontoro@gmail.com)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.

#include "Options.h"

namespace VectorFitting {

Options::Options() {
    relax_                     = true;
    stable_                    = true;
    asymptoticTrend_           = constant;
    polesType_ 				   = lincmplx;
    skipPoleIdentification_    = false;
    skipResidueIdentification_ = false;
//    complexSpaceState_         = true;
}

Options::~Options() {
}

Options::AsymptoticTrend Options::getAsymptoticTrend() const {
    return asymptoticTrend_;
}

Options::PolesType Options::getPolesType() const {
	return polesType_;
}

void Options::setAsymptoticTrend(
        Options::AsymptoticTrend asymptoticTrend) {
    asymptoticTrend_ = asymptoticTrend;
}

void Options::setPolesType(
		Options::PolesType polesType) {
	polesType_ = polesType;
}

bool Options::isRelax() const {
    return relax_;
}

void Options::setRelax(bool relax) {
    relax_ = relax;
}

bool Options::isSkipPoleIdentification() const {
    return skipPoleIdentification_;
}

void Options::setSkipPoleIdentification(
        bool skipPoleIdentification) {
    skipPoleIdentification_ = skipPoleIdentification;
}

bool Options::isSkipResidueIdentification() const {
    return skipResidueIdentification_;
}

void Options::setSkipResidueIdentification(
        bool skipResidueIdentification) {
    skipResidueIdentification_ = skipResidueIdentification;
}

bool Options::isStable() const {
    return stable_;
}

void Options::setStable(bool stable) {
    stable_ = stable;
}

//bool VectorFitting::Options::isComplexSpaceState() const {
//    return complexSpaceState_;
//}
//
//void VectorFitting::Options::setComplexSpaceState(bool complexSpaceState) {
//    complexSpaceState_ = complexSpaceState;
//}

} /* namespace VectorFitting */

