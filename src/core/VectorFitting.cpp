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
        return a.first.imag() < b.first.imag();
    }
} sampleOrdering;

// Quick check to see if a Complex number is real
bool isReal(Complex n, Real tol = 1e-20){ return(n.imag() < tol); }


void VectorFitting::init(const vector<Sample>& samples,
    const vector<Complex>& poles,
    size_t order) {
        // Sanity check: the complex poles should come in pairs; otherwise, there
        // is an error
        Complex currentPole, conjugate;
        for (size_t i = 0; i < poles.size(); i++) {
            currentPole = poles[i];

            if(!isReal(currentPole)){
                assert(conj(currentPole) == poles[i+1]);
                i++;
            }
        }


        samples_ = samples;
        poles_ = poles;
        order_ = order;
    }


VectorFitting::VectorFitting(const vector<Sample>& samples,
                             const vector<Complex>& poles,
                             size_t order) {
     init(samples, poles, order);
}



VectorFitting::VectorFitting( const vector<Sample>& samples,
                              size_t order) {
    // The startin poles are all complex, so the order has to be even.
    assert(order % 2 == 0);

    // Define starting poles as a vector of complex conjugates -a + bi with
    // the imaginary part linearly distributed over the frequency range of
    // interest; i.e., for each pair of complex conjugates (see eqs 9 and 10):
    //      1. imagParts = linspace(range(samples), number_of_poles)
    //      2. realParts = imagParts / 100

    // Get range of the samples frequencies:
    Real min = (*min_element(samples.begin(), samples.end(),
                sampleOrdering)).first.imag();
    Real max = (*max_element(samples.begin(), samples.end(),
                sampleOrdering)).first.imag();
    pair<Real, Real> range(min, max);

    // Generate the imaginary parts of the initial poles from a linear
    // distribution covering the range in the samples.
    // This can also be done with a logarithmic distribution (sometimes
    // faster convergence -see Userguide, p.8-)
    vector<Real> imagParts = linspace(range, order/2);

    cout << "HOLA" << endl;
    // Generate all the starting poles
    vector<Complex> poles(order);
    for (size_t i = 0; i < order; i+=2) {
        Real imag = imagParts[i];
        Real real = - imag / (Real) 100.0;
        poles_[i] = Complex(real, imag);
        poles_[i+1] = conj(poles_[i]);
    }

    init(samples, poles, order);
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
                A(k, i) = Complex(1.0, 0.0) / (s_k - a_i);
                A(k, i + order + 2) = -A(k,i) * f_k;
            }
            // Complex pair
            else{
                // Sanity check
                assert(conj(a_i) == startingPoles[i+1]);

                A(k,i) =   Complex(1.0, 0.0) / (s_k - a_i) + Complex(1.0, 0.0) / (s_k - conj(a_i));
                A(k,i+1) = Complex(0.0, 1.0) / (s_k - a_i) - Complex(0.0, 1.0) / (s_k - conj(a_i));

                A(k, i + order + 2)     = -A(k,i) * f_k;
                A(k, i + order + 2 + 1) = -A(k,i+1) * f_k;

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

    cout << "A:" << endl << Ap << endl;
    cout << "B:" << endl << Bp << endl;


    // Solve A'X = B'. The solution X is the following:
    //      X[:N] = f residues
    //      X[N] = d
    //      X[N+1] = h
    //      X[N+2:] = sigma residues <- these values are the important ones
    VectorXd X = Ap.fullPivHouseholderQr().solve(Bp);

    VectorXd sigmaResidues = X.tail(order);

    cout << "residues: " << endl << sigmaResidues << endl;

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
    //      If a_i is real => C_[i] = sigma_residues[i] (which is real)
    //      If a_i, a_{i+1} is a complex conjugate =>
    //          C_[i] = sigma_residues[i] (which is the real part)
    //          C_[i+1] = sigma_residues[i+1] (which is the imaginary part)

    MatrixXd A_ = MatrixXd::Zero(order, order);
    VectorXd B_(order);
    RowVectorXd C_(sigmaResidues);

    // Populate matrices A_ and B_
    for (size_t i = 0; i < order; i++) {
        Complex a_i = startingPoles[i];

        // Real pole
        if(isReal(a_i)){
            A_(i,i) = a_i.real();
            B_(i) = 1.0;
        }
        // Complex pole
        else{
            // Sanity check
            assert(conj(a_i) == startingPoles[i+1]);

            A_(i, i) = A_(i+1, i+1) = a_i.real();

            A_(i, i+1) = a_i.imag();
            A_(i+1, i) = -a_i.imag();

            B_(i) = 2.0;
            B_(i+1) = 0.0;

            // FIXME: Ugly, beugh!
            i++;
        }
    }

    // Define matrix H = A_ - B_C_, whose eigenvalues are the fitted poles
    MatrixXd H = A_ - B_*C_;

    // Compute the eigenvalues!
    EigenSolver<MatrixXd> eigenSolver(H);
    VectorXcd poles = eigenSolver.eigenvalues();

    // Force unstable poles to be stable
    for (size_t i = 0; i < poles.size(); i++) {
        if(poles(i).real() > 0)
            poles(i) = Complex(-poles(i).real(), poles(i).imag());
    }

    cout << "POLES:" << endl << poles << endl;

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
    // Number of columns = order of approximation + 2
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
                A(k, i) = Complex(1.0, 0.0) / (s_k - a_i);
            }
            // Complex pair
            else{
                // Sanity check: complex poles must come in pairs
                assert(conj(a_i) == poles[i+1]);

                A(k,i) =   Complex(1.0, 0.0) / (s_k - a_i) + Complex(1.0, 0.0) / (s_k - conj(a_i));
                A(k,i+1) = Complex(0.0, 1.0) / (s_k - a_i) - Complex(0.0, 1.0) / (s_k - conj(a_i));

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
    VectorXd X = Ap.fullPivHouseholderQr().solve(Bp);

    vector<Complex> residues(order);

    cout << "RESIDUES:" << endl;

    for (size_t i = 0; i < order; i++) {
        Complex a_i = poles[i];

        // Real pole
        if(isReal(a_i)){
            residues[i] = Complex(X(i), 0.0);
            cout << residues[i] << endl;
        }
        // Complex pole
        else{
            // Sanity check: complex poles must come in pairs
            assert(conj(a_i) == poles[i+1]);

            Real real = X(i);
            Real imag = X(i+1);

            residues[i] = Complex(real, imag);
            residues[i+1] = Complex(real, -imag);
            cout << residues[i] << endl;
            cout << residues[i+1] << endl;

            i++;
        }

    }

    d_ = X(order);
    h_ = X(order + 1);

    cout << "PARAMETERS:" << endl;
    cout << d_ << endl;
    cout << h_ << endl;

    // Return the fitted residues
    return(residues);
}

void VectorFitting::fit(){
    // for (size_t i = 0; i < 1; i++) {
        poles_ = poleIdentification(poles_);
        residues_ = residueIdentification(poles_);
    // }
}

// Return the fitted samples: a vector of pairs s <-> f(s), where f(s) is
// computed with the model in (2)
vector<Sample> VectorFitting::getFittedSamples() const {
    // Vector to store the fitted samples
    vector<Sample> fittedSamples(samples_.size());

    for (size_t k = 0; k < samples_.size(); k++) {
        // Independent variable s
        Complex sk = samples_[k].first;

        // Response of the model
        vector<Complex> response = predictResponse(sk);

        // Building of the sample
        fittedSamples[k] = Sample(sk, response);
    }

    return fittedSamples;
}

// TODO: This is not yet a vector fitting
vector<Complex> VectorFitting::predictResponse(Complex freq) const {
    // Computation of the response with the fitted model (see (2))
    vector<Complex> response(1);
    response[0] = Complex(0,0);

    for (size_t n = 0; n < order_; n++) {
        Complex an = poles_[n];
        Complex cn = residues_[n];

        response[0] += cn / (freq - an);
    }

    response[0] += d_ + freq * h_;

    return response;
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

    Real error = 0.0;
    Complex actual, fitted, diff;

    for (size_t i = 0; i < samples_.size(); i++) {
        // Sanity check: the response should be on the *same* frequency
        assert(actualSamples[i].first == fittedSamples[i].first);

        // Retrieve the actual and fitted responses
        actual = actualSamples[i].second[0];
        fitted = fittedSamples[i].second[0];

        diff = actual - fitted;

        error += abs(diff * diff);
    }

    return sqrt(error/samples_.size());
}

} /* namespace VectorFitting */
