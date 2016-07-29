// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
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

#ifndef SEMBA_MATH_FITTING_VECTOR_H_
#define SEMBA_MATH_FITTING_VECTOR_H_

#include <vector>
#include <complex>
#include <eigen3/Eigen/Dense>

#include "Types.h"

namespace VectorFitting {

using namespace Eigen;
using namespace std;

    typedef complex<Real> Complex;
    typedef pair<Complex, vector<Complex>> Sample;

class VectorFitting {
public:

    // TODO: manage options in constructor.
    /**
     * @param samples   Data to be fitted.
     * @param N         Order of approximation. It shall be an even number.
     */
    VectorFitting(const vector<Sample>& samples, size_t order = 20);

    // This could be called from the constructor, but if an iterative algorithm
    // is preferred, it's a good idea to have it as a public method
    void fit();

    vector<Sample> getFittedSamples(
            const vector<complex<Real>>& frequencies) const;
    vector<complex<Real> > getPoles();
    vector<complex<Real> > getResidues();
    Real getRMSE();

private:
    vector<Sample> samples_;

    vector<Complex> poles_, residues_;

    Real d_, h_;

    Matrix<Complex,Dynamic,Dynamic>  A_, C_;
    vector<Real>  B_, D_, E_;

    size_t order_;

    /**
     * Perfoms the first stage of the algorithm: given a set of starting
     * poles, it returns the new fitted poles.
     * @param startingPoles Vector of complex numbers (Real as its base
     *                      type) containing the starting poles.
     * @return A vector of complex numbers with Real as its base type
     *           containing the fitted poles. This set of poles should
     *           be used to feed the residue identification method.
     */
    vector<Complex> poleIdentification(const vector<Complex>& startingPoles);

    /**
     * Perfoms the second stage of the algorithm: given a set of known
     * poles, it returns the new fitted residues.
     * @param poles Vector of complex numbers (Real as its base type)
     *              containing the known poles.
     * @return A vector of complex numbers with Real as its base type
     *           containing the fitted residues.
     */
    vector<Complex> residueIdentification(const vector<Complex>& poles);
};

} /* namespace VectorFitting */

#endif // SEMBA_MATH_FITTING_VECTOR_H_
