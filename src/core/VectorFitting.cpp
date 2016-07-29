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

#include "VectorFitting.h"
#include "SpaceGenerator.h"

#include <iostream>

namespace VectorFitting {

// Custom ordering for the samples, depending on the imaginary parts of the
// frequencies
struct {
    bool operator()(Sample a, Sample b)
    {
        Complex aS = a.first;
        Complex bS = b.first;

        if(aS.imag() == bS.imag())
            return(aS.real() < bS.real());
        else
            return(aS.imag() < bS.imag());
    }
} sampleOrdering;

bool isReal(Complex n){ return(n.imag() == 0); }

// TODO: Sanity check on the poles: if there are complex ones, they shall come
// in pairs
VectorFitting::VectorFitting(const vector<Sample>& samples,
                             const vector<Complex>& poles,
                             size_t order) {
    samples_ = samples;

    // The samples should be ordered depending on their imaginary parts.
    sort(samples_.begin(), samples_.end(), sampleOrdering);

    order_ = order;
    poles_ = poles;
}


VectorFitting::VectorFitting( const vector<Sample>& samples,
                              size_t order) {
    samples_ = samples;

    // The samples should be ordered depending on their imaginary parts.
    sort(samples_.begin(), samples_.end(), sampleOrdering);

    order_ = order;

    // Define starting poles as a vector of complex conjugates -a + bi with
    // the imaginary part linearly distributed over the frequency range of
    // interest; i.e., for each pair of complex conjugates (see eqs 9 and 10):
    //      1. imagParts = linspace(range(samples), number_of_poles)
    //      2. realParts = imagParts / 100

    // Generate the imaginary parts of the initial poles from a linear
    // distribution covering the range in the samples.
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

    int order = startingPoles.size();

    // Number of rows = number of samples
    // Number of columns = 2* order of approximation + 2
    MatrixXcd A(samples_.size(), 2*order + 2);
    VectorXcd B(samples_.size());

    // TODO: We are dealing only with the first element in f(s_k), so this
    // is not yet a *vector* fitting.
    for (size_t k = 0; k < samples_.size(); k++) {
        Complex s_k = samples_[k].first;
        Complex f_k = samples_[k].second[0];

        B(k) = f_k;

        for (size_t i = 0; i < order; i++) {
            Complex a_i = startingPoles[i];

            // Real pole
            if(isReal(a_i)){
                A(k, i) = 1.0 / (s_k - a_i);
                A(k, 2 + 2*i) = -A(k,i) * f_k;
            }
            // Complex pair
            else{
                A(k,i) = 1.0 / (s_k - a_i) + 1.0 / (s_k - conj(a_i));
                A(k,i+1) = Complex(0.0,1.0) * (1.0 / (s_k - a_i) - 1.0 / (s_k - conj(a_i)));

                A(k,2 + 2*i)        = -A(k,i) * f_k;
                A(k,2 + 2*i + 1)    = -A(k,i+1) * f_k;

                // FIXME: Ugly, beugh!
                i++;
            }
        }

        A(k, order) = Complex(1.0,0.0);
        A(k, order+1) = s_k;
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
    VectorXd sigmaResidues = X.tail(X.size() - (order + 2));

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

    MatrixXd A_ = MatrixXd::Zero(order, order);
    VectorXd B_(order);
    RowVectorXd C_(order);

    // Populate matrices A_, B_ and C_
    for (size_t i = 0; i < order; i++) {
        Complex pole = startingPoles[i];
        Complex residue = sigmaResidues(i);

        // Real pole
        if(isReal(pole)){
            A_(i,i) = pole.real();
            B_(i) = 1.0;
            C_(i) = residue.real();
        }
        // Complex pole
        else{
            A_(i, i) = A_(i+1, i+1) = pole.real();

            A_(i, i+1) = pole.imag();
            A_(i+1, i) = -pole.imag();

            B_(i) = 2.0;
            B_(i+1) = 0.0;

            C_(i) = residue.real();
            C_(i+1) = residue.imag();

            // FIXME: Ugly, beugh!
            i++;
        }
    }

    // Define matrix H = A - BC, whose eigenvalues are the fitted poles
    MatrixXd H = A_ - B_*C_;

    // Compute the eigenvalues!
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

    int order = poles.size();

    // Number of rows = number of samples
    // Number of columns = 2* order of approximation + 2
    MatrixXcd A(samples_.size(), order + 2);
    VectorXcd B(samples_.size());

    // TODO: We are dealing only with the first element in f(s_k), so this
    // is not yet a *vector* fitting.
    for (size_t k = 0; k < samples_.size(); k++) {
        Complex s_k = samples_[k].first;
        Complex f_k = samples_[k].second[0];

        B(k) = f_k;

        for (size_t i = 0; i < order; i++) {
            Complex a_i = poles[i];

            // Real pole
            if(isReal(a_i)){
                A(k, i) = 1.0 / (s_k - a_i);
            }
            // Complex pair
            else{
                A(k,i) = 1.0 / (s_k - a_i) + 1.0 / (s_k - conj(a_i));
                A(k,i+1) = Complex(0.0,1.0) * (1.0 / (s_k - a_i) - 1.0 / (s_k - conj(a_i)));

                // FIXME: Ugly, beugh!
                i++;
            }
        }

        A(k, order) = Complex(1.0,0.0);
        A(k, order+1) = s_k;
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

    vector<Complex> residues(order);

    for (size_t i = 0; i < order; i++) {
        Complex pole = poles[i];

        // Real pole
        if(isReal(pole)){
            residues[i] = Complex(X(i),0.0);
        }
        // Complex pole
        else{
            Real real = X(i);
            Real imag = X(i+1);

            residues[i] = Complex(real, imag);
            residues[i+1] = Complex(real, -imag);

            i++;
        }

    }

    d_ = X(order);
    h_ = X(order + 1);

    // Return the fitted residues
    return(residues);
}

void VectorFitting::fit(){
    for (size_t i = 0; i < 30; i++) {
        poles_ = poleIdentification(poles_);
        residues_ = residueIdentification(poles_);
    }
}

// Return the fitted samples: a vector of pairs s <-> f(s), where f(s) is
// computed with the model in (2)
// TODO: This is not yet a *vector* fitting, the [0] has to be removed :)
vector<Sample> VectorFitting::getFittedSamples() const {

    // Vector to store the fitted samples and loop variable
    vector<Sample> fittedSamples(samples_.size());

    for (size_t k = 0; k < samples_.size(); k++) {
        // Independent variable s
        Complex sk = samples_[k].first;

        // Computation of the dependent variable f(s) with the model (see (2))
        vector<Complex> f(1);
        f[0] = Complex(0,0);

        for (size_t n = 0; n < order_; n++) {
            Complex cn = residues_[n];
            Complex an = poles_[n];

            f[0] += cn / (sk - an);
        }

        f[0] += d_ + sk * h_;

        fittedSamples[k] = Sample(sk, f);
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
    vector<Sample> actualSamples = samples_;
    vector<Sample> fittedSamples = getFittedSamples();

    // The samples should be ordered in the same way
    sort(actualSamples.begin(), actualSamples.end(), sampleOrdering);
    sort(fittedSamples.begin(), fittedSamples.end(), sampleOrdering);

    Real error = 0.0;
    Complex actual, fitted, diff;

    for (size_t i = 0; i < samples_.size(); i++) {
        actual = actualSamples[i].second[0];
        fitted = fittedSamples[i].second[0];

        diff = actual - fitted;

        error += abs(diff * diff);
    }

    return sqrt(error/samples_.size());
}

} /* namespace VectorFitting */
