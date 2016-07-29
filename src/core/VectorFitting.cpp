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

#include "VectorFitting.h"
#include "SpaceGenerator.h"

#include <iostream>

namespace VectorFitting {

VectorFitting::VectorFitting( const vector<Sample>& samples,
                              size_t order) {

    samples_ = samples;
    order_ = order;

    // Define starting poles as a vector of complex conjugates -a + bi with
    // the imaginary part linearly distributed over the frequency range of
    // interest; i.e., for each pair of complex conjugates (see eqs 9 and 10):
    //      1. imagParts = linspace(range(samples), number_of_poles)
    //      2. realParts = imagParts / 100

    // Generate the imaginary parts of the initial poles from a linear
    // distribution covering the range in the samples.
    // TODO: We are assuming the samples are ordered depending on their
    // imaginary parts; order them before doing this!
    pair<Real,Real> range(  samples_.front().first.imag(),
                            samples_.back().first.imag());

    // This can also be done with a logarithmic distribution (sometimes
    // faster convergence -see Userguide, p.8-)
    vector<Real> imagParts = linspace(range, order_/2);

    // Generate all the starting poles
    for (size_t i = 0; i < order_; i+=2) {
        Real imag = imagParts[i];
        Real real = - imag / (Real) 100.0;
        poles_.push_back(Complex(real, imag));
        poles_.push_back(Complex(real, -imag));
    }
}

vector<Complex> VectorFitting::poleIdentification(const vector<Complex>& startingPoles){
    // Define matrix A following equation (A.3), where:
    //      s_k = sample[k].first()
    //      f(s_k) = sample[k].second()
    //      a_i = starting poles
    // If a_i, a_{i+1} is a complex conjugate (they always come in pairs;
    // otherwise, there is an error), the elements have to follow equation
    // (A.6)
    // Define column vector B following equation (A.4), where:
    //      f(s_k) = sample.second()


    // Number of rows = number of samples
    // Number of columns = 2* order of approximation + 2
    MatrixXcd A(samples_.size(), 2*order_ + 2);
    VectorXcd B(samples_.size());

    // TODO: We are dealing only with the first element in f(s_k), so this
    // is not yet a *vector* fitting.
    for (size_t k = 0; k < A.rows()/2; k++) {
        Complex s_k = samples_[k].first;
        Complex f_k = samples_[k].second[0];

        B(k,0) = f_k;

        // TODO: We are considering all the poles are complex, is it ok?
        for (size_t i = 0; i < order_; i+=2) {
            Complex a_i = startingPoles[i];
            A(k,i) = 1.0 / (s_k - a_i) + 1.0 / (s_k - conj(a_i));
            A(k,i+1) = Complex(0.0,1.0) * A(k,i);

            A(k,2 + 2*i) = -A(k,i) * f_k;
            A(k,2 + 2*i + 1) = -A(k,i+1) * f_k;
        }
        A(k, order_) = Complex(1.0,0.0);
        A(k, order_+1) = s_k;
    }

    // To follow (A.8) note, it is necessary to build A' and B' as a stack
    // of the real and imaginary parts of A and B
    MatrixXd Ap(2*A.rows(), A.cols());
    VectorXd Bp(2*B.size());

    Ap << A.real(),
          A.imag();

    Bp << B.real(),
          B.imag();

    // Solve A'X = B'. The solution X is the following:
    //      X[:N] = f residues
    //      X[N] = d
    //      X[N+1] = h
    //      X[N+2:] = sigma residues <- these values are the important ones
    VectorXd X = Ap.colPivHouseholderQr().solve(Bp);

    // VectorXd fResidues = X.head(order_);
    // Real d = X[order_];
    // Real h = X[order_ + 1];
    VectorXd sigmaResidues = X.tail(X.size() - (order_ + 2));

    // Debug check
    assert(sigmaResidues.size() == startingPoles.size());

    // Define a diagonal matrix containing the starting poles, A_, (see B.2)
    // as follows:
    //      If a_i is real => A_[i,i] = a_i
    //      If a_i, a_{i+1} is a complex conjugate =>
    //          A_[i,i] = A_[i+1, i+1] = real(a_i)
    //          A_[i,i+1] = imag(a_i)
    //          A_[i+1,i] = -imag(a_i)
    // Define a column vector, B_, as follows (see B.2):
    //      If a_i is real => B_[i] = 1
    //      If a_i, a_{i+1} is a complex conjugate => B_[i] = 2; B_[i+1] = 0
    // Define a row vector, C_, as follows:
    //      If a_i is real => C_[i] = sigma_residues[i]
    //      If a_i, a_{i+1} is a complex conjugate =>
    //          C_[i] = real(sigma_residues[i])
    //          C_[i+1] = imag(sigma_residues[i])

    MatrixXd A_(startingPoles.size(), startingPoles.size());
    VectorXd B_(startingPoles.size());
    RowVectorXd C_ = sigmaResidues; // TODO: We're considering only complex poles

    // Populate matrices A_, B_ and C_
    for (size_t i = 0; i < A_.rows(); i+=2) {
        Real poleReal = startingPoles[i].real();
        Real poleImag = startingPoles[i].imag();

        A_(i, i) = A_(i+1, i+1) = poleReal;

        A_(i, i+1) = poleImag;
        A_(i+1, i) = -poleImag;

        B(i) = 2;
        B(i+1) = 0;
    }

    // Compute the eigen values of H = A - BC, which are the new poles
    MatrixXd H = A_ - B_*C_;

    EigenSolver<MatrixXd> eigenSolver(H);

    VectorXcd poles = eigenSolver.eigenvalues();

    // Return the fitted poles as a std::vector of Complex
    return(vector<Complex>(poles.data(),
                           poles.data() + poles.rows() * poles.cols()));
}

vector<Complex> VectorFitting::residueIdentification(const vector<Complex>& poles){
    // Define matrix A following equation (A.3) without the negative terms, where:
    //      s_k = sample[k].first()
    //      f(s_k) = sample[k].second()
    //      a_i = starting poles
    // If a_i, a_{i+1} is a complex conjugate (they always come in pairs;
    // otherwise, there is an error), the elements have to follow equation
    // (A.6)
    // Define column vector B following equation (A.4), where:
    //      f(s_k) = sample.second()


    // Number of rows = number of samples
    // Number of columns = 2* order of approximation + 2
    MatrixXcd A(samples_.size(), order_ + 2);
    VectorXcd B(samples_.size());

    // TODO: We are dealing only with the first element in f(s_k), so this
    // is not yet a *vector* fitting.
    for (size_t k = 0; k < A.rows()/2; k++) {
        Complex s_k = samples_[k].first;
        Complex f_k = samples_[k].second[0];

        B(k) = f_k;

        // TODO: We are considering all the poles are complex, is it ok?
        for (size_t i = 0; i < order_; i+=2) {
            Complex a_i = poles[i];
            A(k,i) = 1.0 / (s_k - a_i) + 1.0 / (s_k - conj(a_i));
            A(k,i+1) = Complex(0.0,1.0) * A(k,i);
        }

        A(k, order_) = Complex(1.0,0.0);
        A(k, order_+1) = s_k;
    }

    // To follow (A.8) note, it is necessary to build A' and B' as a stack
    // of the real and imaginary parts of A and B
    MatrixXd Ap(2*A.rows(), A.cols());
    VectorXd Bp(2*B.size());

    Ap << A.real(),
          A.imag();

    Bp << B.real(),
          B.imag();

    // Solve A'X = B' (¿by singular value decomposition? Eigen!), whose elements
    // are:
    //      X[:N] = f residues
    //      X[N] = d
    //      X[N+1] = h
    VectorXd X = Ap.colPivHouseholderQr().solve(Bp);

    vector<Complex> residues;

    for (size_t i = 0; i < order_; i+=2) {
        Real real = X(i);
        Real imag = X(i+1);

        residues.push_back(Complex(real, -imag));
        residues.push_back(Complex(real, imag));
    }

    d_ = X(order_);
    h_ = X(order_ + 1);

    // Return the fitted residues
    return(residues);
}

void VectorFitting::fit(){
    poles_ = poleIdentification(poles_);
    residues_ = residueIdentification(poles_);
}

// Return the fitted samples: a vector of pairs s <-> f(s), where f(s) is
// computed with the model in (2)
// TODO: This is not yet a *vector* fitting, the [0] has to be removed :)
vector<Sample> VectorFitting::getFittedSamples(
        const vector<Complex >& frequencies) const {

    // Vector to store the fitted samples and loop variable
    vector<Sample> fittedSamples;
    Sample sample;

    for (size_t i = 0; i < samples_.size(); i++) {
        // Independent variable s
        sample.first = samples_[i].first;

        // Computation of the dependent variable f(s) with the model (see (2))
        sample.second = vector<Complex>();
        sample.second.push_back(Complex(0,0));

        for (size_t j = 0; j < order_; j++) {
            sample.second[0] += residues_[j] / (sample.first - poles_[j]);
        }

        sample.second[0] += d_ + sample.first * h_;

        fittedSamples.push_back(sample);
    }

    return fittedSamples;
}

vector<complex<Real>> VectorFitting::getPoles() {
    // assert(computed_ == true);
    return poles_;
}

vector<complex<Real>> VectorFitting::getResidues() {
    // assert(computed_ == true);
    return residues_;
}

/**
 * Returns the error of the model, measured as the root mean
 * square of the estimated data with respect to the samples.
 * @return Real - Root mean square error of the model.
 */
Real VectorFitting::getRMSE() {
    vector<Sample> fittedSamples = getFittedSamples(vector<Complex>());

    Real error = 0.0;
    Complex actual, fitted;

    for (size_t i = 0; i < samples_.size(); i++) {
        actual = samples_[i].second[0];
        fitted = fittedSamples[i].second[0];

        cout << actual << "\t" << fitted << endl;

        error += abs((actual - fitted) * (actual - fitted));
    }

    return sqrt(error/samples_.size());
}

} /* namespace VectorFitting */
