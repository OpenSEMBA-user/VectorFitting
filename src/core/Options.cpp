// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Alejandro García Montoro        (alejandro.garciamontoro@gmail.com)
//					  Alejandra López de Aberasturi Gómez (aloaberasturi@ugr.es)
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
    n_                         = 0;
    nu_                        = 1e-3;
    polesType_ 				   = PolesType::lincmplx;
    iterations_.first 		   = 4;   // Initial poles.
    iterations_.second   	   = 4;

    relax_                     = true;
    stable_                    = true;
    asymptoticTrend_           = AsymptoticTrend::constant;
    weighting_                 = Weighting::one;
    skipPoleIdentification_    = false;
    skipResidueIdentification_ = false;
    complexSpaceState_         = true;
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

bool VectorFitting::Options::isComplexSpaceState() const {
    return complexSpaceState_;
}

void VectorFitting::Options::setComplexSpaceState(bool complexSpaceState) {
    complexSpaceState_ = complexSpaceState;
}

Options::Weighting Options::getWeighting() const {
    return weighting_;
}

std::pair<size_t, size_t> Options::getIterations() const {
    return iterations_;
}

void Options::setIterations(const std::pair<size_t, size_t>& iterations) {
    iterations_ = iterations;
}

size_t Options::getN() const {
    return n_;
}

void Options::setN(size_t n) {
    n_ = n;
}

double Options::getNu() const {
    return nu_;
}

void Options::setNu(double nu) {
    nu_ = nu;
}

} /* namespace VectorFitting */


