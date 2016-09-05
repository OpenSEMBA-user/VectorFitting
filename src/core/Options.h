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

namespace VectorFitting {

class Options {
public:
    enum AsymptoticTrend {
        zero,
        constant,
        linear
    };

    Options();
    virtual ~Options();

    AsymptoticTrend getAsymptoticTrend() const;
    bool isRelax() const;
    bool isSkipPoleIdentification() const;
    bool isSkipResidueIdentification() const;
    bool isStable() const;

    void setAsymptoticTrend(AsymptoticTrend asymptoticTrend);
    void setRelax(bool relax);
    void setSkipPoleIdentification(bool skipPoleIdentification);
    void setSkipResidueIdentification(bool skipResidueIdentification);
    void setStable(bool stable);

private:
    bool relax_;
    bool stable_;
    AsymptoticTrend asymptoticTrend_;
    bool skipPoleIdentification_;
    bool skipResidueIdentification_;
};

} /* namespace VectorFitting */

#endif // SEMBA_MATH_FITTING_VECTOR_H_
