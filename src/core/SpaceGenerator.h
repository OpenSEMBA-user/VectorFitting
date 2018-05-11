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


#ifndef SEMBA_MATH_UTIL_SPACEGENERATOR_H_
#define SEMBA_MATH_UTIL_SPACEGENERATOR_H_

#include <utility>
#include <vector>

#include "Types.h"

namespace VectorFitting {

template<class T>
std::vector<Real> logspace(const std::pair<Real, Real>& rangeExponents,
                           const T nPoints);

template<class T>
std::vector<T> linspace(const std::pair<T,T>& range,
                        const std::size_t nPoints);

}

#include "SpaceGenerator.hpp"

#endif /* SEMBA_MATH_UTIL_SPACEGENERATOR_H_ */
