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

#ifndef SEMBA_MATH_FITTING_VECTOR_H_
#define SEMBA_MATH_FITTING_VECTOR_H_

#include <vector>
#include <complex>
#include <eigen3/Eigen/Dense>

#include "Real.h"
#include "Options.h"

namespace VectorFitting {

using namespace Eigen;

typedef std::complex<Real> Complex;

/**
 * Samples are formed by a pair formed by:
 *  - First, the parameter $s = j \omega$ a purely imaginary number.
 *  - Second, a vector with the complex data to be fitted.
 */
typedef std::pair<Complex, std::vector<Complex>> Sample;

class Fitting {
public:

    /**
     * Build a fitter with starting poles computed automatically.
     * @param samples   Data to be fitted.
     * @param order     Order of approximation.
     * @param options   Options.
     */
    Fitting(const std::vector<Sample>& samples,
            const size_t order,
            const Options& options,
            const std::vector<std::vector<Real>>& weights =
                    std::vector<std::vector<Real>>());

    /**
     * Build a fitter with starting poles provided by the user. order_ and
     * poles.size() shall be the same
     * @param samples   Data to be fitted.
     * @param poles     Starting poles.
     * @param options   Options.
     */
    Fitting(const std::vector<Sample>& samples,
            const std::vector<Complex>& poles,
            const Options& options,
            const std::vector<std::vector<Real>>& weight =
                    std::vector<std::vector<Real>>());
    /**
     * Build a fitter with starting weights provided by Driver::initWeightsSum()
     */

    Fitting(const Sample& samples,
            const std::vector<Complex>& poles,
            const Options& options,
            const MatrixXd weightsSum);

    // This could be called from the constructor, but if an iterative algorithm
    // is preferred, it's a good idea to have it as a public method
    void fit();
    MatrixXd initWeights(std::vector<std::vector<Real>>& weights);

    std::vector<Sample>  getFittedSamples() const;
    std::vector<Complex> getPoles();

    /**
     *  Getters and setters to fitting coefficents.
     */
    MatrixXcd getA() {return A_;}    // Size:  N, N.
    MatrixXcd getC() {return C_;}    // Size:  Nc, N.
    std::vector<MatrixXcd> getR() {return R_;} 	 // Size:  Nc,Nc,N.
    RowVectorXi getB() {return B_;}  // Size:  1, N.
    VectorXcd getD() {return D_;}    // Size:  1, Nc.
    VectorXcd getE() {return E_;}    // Size:  1, Nc.
    Real getRMSE() const;
    Real getMaxDeviation() const;
    void setOptions(const Options& options);
    void setR(std::vector<MatrixXcd> R);


private:
    Options options_;

    std::vector<Sample> samples_;
    VectorXcd poles_;

    MatrixXcd A_, C_;
    std::vector<MatrixXcd> R_;
    VectorXcd D_, E_;
    RowVectorXi B_;

    MatrixXd weights_; // Size: Ns, Nc

    static constexpr Real toleranceLow_  = 1e-18;
    static constexpr Real toleranceHigh_ = 1e+18;

    void init(const std::vector<Sample>& samples,
              const std::vector<Complex>& poles,
              const Options& options);

    size_t getSamplesSize() const;
    size_t getResponseSize() const;
    size_t getOrder() const;

    static RowVectorXi getCIndex(const VectorXcd& poles);
};

} /* namespace VectorFitting */

#endif // SEMBA_MATH_FITTING_VECTOR_H_
