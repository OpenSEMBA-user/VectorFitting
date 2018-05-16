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


#ifndef SEMBA_MATH_FITTING_VECTOR_H_
#define SEMBA_MATH_FITTING_VECTOR_H_

#include <vector>
#include <complex>
#include <algorithm>
#include <eigen3/Eigen/Dense>

#include "Real.h"
#include "Options.h"

namespace VectorFitting {

using namespace Eigen;

typedef std::complex<Real> Complex;


class Fitting {

public:
	/**
	 * Samples are formed by a pair formed by:
	 *  - First, the parameter $s = j \omega$ a purely imaginary number.
	 *  - Second, a vector with the complex data to be fitted.
	 */
	typedef std::pair<Complex, VectorXcd> Sample;

	Fitting() {}

	/**
     * Build a fitter with starting poles provided by the user. order_ and
     * poles.size() shall be the same
     * @param samples   Data to be fitted.
     * @param options   Options.
     * @param poles     Starting poles (optional).
     * @param weights   Samples weights (optional).
     */
    Fitting(const std::vector<Sample>& samples,
            const Options& options,
            const std::vector<Complex>& poles = {},
			const std::vector<VectorXd>& weights = {});


    // This could be called from the constructor, but if an iterative algorithm
    // is preferred, it's a good idea to have it as a public method
    void fit();

    std::vector<Sample>  getFittedSamples() const;

    std::vector<Complex> getPoles();

    /**
     *  Getters and setters to fitting coefficents.
     */
    MatrixXcd getA() {return A_;}    // Size:  N, N.
    MatrixXcd getC() {return C_;}    // Size:  Nc, N.
    VectorXi getB()  {return B_;}    // Size:  1, N.
    VectorXcd getD() {return D_;}    // Size:  1, Nc.
    VectorXcd getE() {return E_;}    // Size:  1, Nc.
    Real getRMSE() const;
    Real getMaxDeviation() const;
	const std::vector<Sample>& getSamples() const;

    Options& options() {return options_;};

	void setA(const MatrixXcd& a) {A_ = a;}
	void setB(const MatrixXi&  b) {B_ = b;}
	void setC(const MatrixXcd& c) {C_ = c;}
	void setD(const VectorXcd& d) {D_ = d;}
	void setE(const VectorXcd& e) {E_ = e;}

    template <class T>
    static std::vector<T> toStdVector(
            const Matrix<T, Eigen::Dynamic, 1>& rhs) {
        std::vector<T> res(rhs.size());
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] = rhs(i);
        }
        return res;
    }

    template <class T>
    static Eigen::Matrix<T, Eigen::Dynamic, 1> toEigenVector(
            const std::vector<T>& rhs) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> res(rhs.size());
        for (size_t i = 0; i < rhs.size(); ++i) {
            res(i) = rhs[i];
        }
        return res;
    }

private:
    Options options_;

    std::vector<Sample> samples_;
    std::vector<Complex> poles_;

    MatrixXcd A_, C_;
    VectorXcd D_, E_;
    MatrixXi B_;

    std::vector<VectorXd> weights_; // Size: Ns, Nc

    static constexpr Real toleranceLow_  = 1e-4;
    static constexpr Real toleranceHigh_ = 1e+4;

    size_t getSamplesSize() const;
    size_t getResponseSize() const;
    size_t getOrder() const;

    static RowVectorXi getCIndex(const std::vector<Complex>& poles);

    struct {
        bool operator()(Complex a, Complex b)
        {
            if (lower(a.real(), b.real())) {
                return true;
            }
            if (equal(a.real(), b.real())) {
                return lower(a.imag(), b.imag());
            }
            return false;
        }
    } complexOrdering;

    // Quick check to see if a Complex number is real
    static bool isReal(Complex n){
        return equal(n.imag(), 0.0);
    }

    Real useWeight_(size_t i, VectorXd::Index n) {
        if (weights_[i].size() > 1) {
            return weights_[i](n);
        } else if (weights_[i].size() == 1) {
            return weights_[i](0);
        } else {
            throw std::runtime_error("Invalid weight operation");
        }
    }

};

} /* namespace VectorFitting */

#endif // SEMBA_MATH_FITTING_VECTOR_H_
