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

#ifndef SEMBA_VECTOR_FITTING_OPTIONS_H_
#define SEMBA_VECTOR_FITTING_OPTIONS_H_

#include <utility>
#include <cstddef>

namespace VectorFitting {

class Options {
public:
    enum class AsymptoticTrend {
        zero,
        constant,
        linear
    };

    enum class PolesType {
    	lincmplx,
		logcmplx
    };

    enum class Weighting {
        one,
        oneOverAbs,
        oneOverSqrtAbs,
        oneOverNorm,
        oneOverSqrtNorm
    };

    Options();
    virtual ~Options();

    AsymptoticTrend getAsymptoticTrend() const;
    PolesType getPolesType() const;
    Weighting getWeighting() const;
    bool isRelax() const;

    bool isSkipPoleIdentification() const;
    bool isSkipResidueIdentification() const;
    bool isStable() const;
    bool isComplexSpaceState() const;

    void setAsymptoticTrend(AsymptoticTrend asymptoticTrend);
    void setPolesType(PolesType polesType);
    void setRelax(bool relax);
    void setSkipPoleIdentification(bool skipPoleIdentification);
    void setSkipResidueIdentification(bool skipResidueIdentification);
    void setStable(bool stable);
    void setComplexSpaceState(bool complexSpaceState);
    std::pair<std::size_t, std::size_t> getIterations() const;
    void setIterations(const std::pair<std::size_t, std::size_t>& iterations);

    size_t getN() const;
    void setN(size_t n);

    double getNu() const;
    void setNu(double nu);

private:

    bool relax_;
    bool stable_;
    AsymptoticTrend asymptoticTrend_;
    Weighting weighting_;
    bool skipPoleIdentification_;
    bool skipResidueIdentification_;
    bool complexSpaceState_;
    double nu_;

    PolesType polesType_;
    size_t N_;
    std::pair<size_t, size_t> iterations_;
};

} /* namespace VectorFitting */

#endif // SEMBA_MATH_FITTING_VECTOR_H_
