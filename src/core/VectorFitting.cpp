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
bool isReal(Complex n){
    return equal(n.imag(), 0.0);
}

void VectorFitting::init(const std::vector<Sample>& samples,
                         const std::vector<Complex>& poles,
                         const Options& options) {
    options_ = options;

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
    weights_ = MatrixXd::Ones(getSamplesSize(), getResponseSize());
}

VectorFitting::VectorFitting(const std::vector<Sample>& samples,
        const std::vector<Complex>& poles,
        const Options& options) {
    if (samples.size() == 0) {
        throw std::runtime_error("Samples size cannot be zero");
    }
    init(samples, poles, options);
}

VectorFitting::VectorFitting(const std::vector<Sample>& samples,
        const size_t order,
        const Options& options) {
    if (samples.size() == 0) {
        throw std::runtime_error("Samples size cannot be zero");
    }
    if (order % 2 == 0) {
        throw std::runtime_error("Default starting poles are complex, order must be even");
    }

    // Define starting poles as a vector of complex conjugates -a + bi with
    // the imaginary part linearly distributed over the frequency range of
    // interest; i.e., for each pair of complex conjugates (see eqs 9 and 10):
    //      1. imagParts = linspace(range(samples), number_of_poles)
    //      2. realParts = imagParts / 100

    // Get range of the samples frequencies:
    Sample minSample = *min_element(samples.begin(), samples.end(), sampleOrdering);
    Sample maxSample = *max_element(samples.begin(), samples.end(), sampleOrdering);
    std::pair<Real, Real> range(minSample.first.imag(), maxSample.first.imag());

    // Generate the imaginary parts of the initial poles from a linear
    // distribution covering the range in the samples.
    // This can also be done with a logarithmic distribution (sometimes
    // faster convergence -see Userguide, p.8-)
    std::vector<Real> imagParts = linspace(range, order/2);

    // Generate all the starting poles
    std::vector<Complex> poles(order);

    for (size_t i = 0; i < order; i+=2) {
        Real imag = imagParts[i/2];
        Real real = - imag / (Real) 100.0;
        poles[i] = Complex(real, imag);
        poles[i+1] = conj(poles[i]);
    }

    init(samples, poles, options);
}

void VectorFitting::fit(){
    // Following Gustavssen notation in vectfit3.m .
    const size_t Ns = getSamplesSize();
    const size_t N  = getOrder();
    const size_t Nc = getResponseSize();

    MatrixXcd LAMBD(N, N), LAMBDprime(N,N);
    for (size_t i = 0; i < N; ++i) {
        LAMBD(i,i) = poles_[i];
    }
    LAMBDprime = LAMBD.transpose().conjugate();

    // --- Pole identification ---
    if (!options_.isSkipPoleIdentification()) {
        // Finds out which starting poles are complex.
        RowVectorXi cindex = RowVectorXi::Zero(N);
        for (size_t m = 0; m < N; ++m) {
            if (!equal(std::imag(LAMBD(m,m)), 0.0)) {
                if (m == 0) {
                    cindex(m) = 1;
                } else {
                    if (cindex(m-1) == 0 || cindex(m-1) == 2) {
                        cindex(m) = 1;
                        cindex(m+1) = 2;
                    } else {
                        cindex(m) = 2;
                    }
                }
            }
        }

        // Builds system - matrix.
        MatrixXcd Dk(Ns,N+2);
        for (size_t m = 0; m < N; ++m) {
            if (cindex(m) == 0) { // Real pole.
                for (size_t i = 0; i < Ns; ++i) {
                    Dk(i,m) = Complex(1,0) / (samples_[i].first - LAMBD(m,m));
                }
            } else if (cindex(m) == 1) { // Complex pole, first part.
                for (size_t i = 0; i < Ns; ++i) {
                    Dk(i,m)   = Complex(1,0) / (samples_[i].first - LAMBD(m,m))
                               + Complex(1,0) / (samples_[i].first - LAMBDprime(m,m));
                    Dk(i,m+1) = Complex(0,1) / (samples_[i].first - LAMBD(m,m))
                               - Complex(0,1) / (samples_[i].first - LAMBDprime(m,m));
                }
            }
        }
        for (size_t i = 0; i < Ns; ++i) {
            Dk(i,N) = (Real) 1.0;
            if (options_.getAsymptoticTrend() == Options::linear) {
                Dk(i,N+1) = samples_[i].first;
            }
        }
        // Scaling for last row of LS-problem (pole identification).
        Real scale = 0.0;
        for (size_t m = 0; m < Nc; ++m) {
            for (size_t i = 0; i < Ns; ++i) {
                const Real weight = weights_(i,m);
                const Complex sample = samples_[i].second[m];
                scale += std::pow(std::abs(weight * std::conj(sample)), 2);
            }
        }
        scale = std::sqrt(scale) / (Real) Ns;

        if (options_.isRelax()) {
            size_t offs;
            switch (options_.getAsymptoticTrend()) {
            case Options::zero:
                offs = 0;
                break;
            case Options::constant:
                offs = 1;
                break;
            case Options::linear:
                offs = 2;
                break;
            }

            // Computes AA and bb.
            MatrixXd AA(Nc*(N+1), N+1);
            VectorXd bb(Nc*(N+1));
            VectorXd x(Nc*(N+1));
            for (size_t n = 0; n < Nc; ++n) {
                MatrixXd A(2*Ns+1, (N+offs)+N+1);
                VectorXd weig(Ns);
                for (size_t i = 0; i < Ns; ++i) {
                    weig(i) = weights_(i,n);
                }
                // Left block.
                for (size_t m = 0; m < N + offs; ++m) {
                    for (size_t i = 0; i < Ns; ++i) {
                        const Complex entry = weig(i) * Dk(i,m);
                        A(i   ,m) = std::real(entry);
                        A(i+Ns,m) = std::imag(entry);
                    }
                }
                // Right block.
                const size_t inda = N + offs;
                for (size_t m = 0; m < N+1; ++m) {
                    for (size_t i = 0; i < Ns; ++i) {
                        const Complex entry =
                         - weig(i) * Dk(i,m) * samples_[i].second[n];
                        A(i   ,inda+m) = std::real(entry);
                        A(i+Ns,inda+m) = std::imag(entry);
                    }
                }

                // Integral criterion for sigma.
                const size_t offset = N + offs;
                if (n == Nc-1) {
                    for (size_t mm = 0; mm < N+1; ++mm) {
                        A(2*Ns, offset+mm) = std::real(scale*Dk.col(mm).sum());
                    }
                }

                // Performs QR decomposition.
                MatrixXd Q, R;
                HouseholderQR<MatrixXd> qr(A.rows(), A.cols());
                qr.compute(A);
                Q = qr.householderQ() * MatrixXd::Identity(A.rows(),A.cols());
                R = Q.transpose() * A;

                const size_t ind = N + offs;
                MatrixXd R22 = R.block(ind,ind, N+1,N+1);
                AA.block(n*(N+1), 0, N+1, N+1) = R22;
                if (n == Nc-1) {
                    for (size_t i = 0; i < N+1; ++i) {
                        // FIXME: This is buggy.
                        bb(i + n*(N+1)) = Q(2*Ns + i, N+offs)
                                * (Real) Ns * (Real) scale;
                    }
                }
            }  // End of for loop n=1:Nc

            // Computes scaling factor.
            VectorXd Escale(Nc*(N+1));
            for (size_t col = 0; col < N+1; ++col) {
                Escale(col) = 1.0 / AA.col(col).norm();
                for (size_t i = 0; i < Nc*(N+1); ++i) {
                    AA(i,col) = Escale(col) * AA(i,col);
                }
            }

            x = AA.inverse() * bb;
            for (size_t i = 0; i < Nc*(N+1); ++i) {
                x(i) *= Escale(i);
            }

        } // End of if for "relax" flag.

        // TODO: If relax is zero or sigma is out of tolerance range.

    } // End of if for "skip pole identification" flag.

    VectorXcd B(N), SERD(Nc), SERE(Nc);
    RowVectorXcd SERA(1,N);
    MatrixXcd  SERC(Nc, N);
    for (size_t i = 0; i < N; ++i) {
        SERA(0,i) = poles_[i];
    }
}

// Return the fitted samples: a vector of pairs s <-> f(s), where f(s) is
// computed with the model in (2)
std::vector<Sample> VectorFitting::getFittedSamples(
        std::vector<Complex> freqs) const {
    // Vector to store the fitted samples
    std::vector<Sample> fittedSamples(freqs.size());

    for (size_t k = 0; k < freqs.size(); k++) {
        // Independent variable s
        Complex sk = freqs[k];

        // Response of the model
        std::vector<Complex> response = predictResponse(sk);

        // Building of the sample
        fittedSamples[k] = Sample(sk, response);
    }

    return fittedSamples;
}

std::vector<Complex> VectorFitting::predictResponse(Complex freq) const {
//    // Computation of the response with the fitted model (see (2))
//    vector<Complex> response(getResponseSize());
//
//    for (size_t i = 0; i < getResponseSize(); i++) {
//        response[i] = Complex(0,0);
//
//        for (size_t n = 0; n < getOrder(); n++) {
//            Complex an = poles_[n];
//            //TODO: residues should be different for every element in the
//            //response; e.g. residues_[i][n].
//            Complex cn = residues_[n];
//
//            response[i] += cn / (freq - an);
//        }
//
//        // d_ and h_ should have responseSize_ elements; e.g., d_[i], h_[i];
//        response[i] += d_ + freq * h_;
//    }
//
//    return response;
}

std::vector<Complex> VectorFitting::getPoles() {
    return poles_;
}

/**
 * Returns the error of the model, measured as the root mean
 * square of the estimated data with respect to the samples.
 * @return Real - Root mean square error of the model.
 */
Real VectorFitting::getRMSE() {
    std::vector<Sample> actualSamples = samples_;

    // Retrieve the frequencies from the samples
    std::vector<Complex> frequencies(actualSamples.size());
    for (size_t k = 0; k < actualSamples.size(); k++) {
        frequencies[k] = actualSamples[k].first;
    }

    // Get the fitted samples
    std::vector<Sample> fittedSamples = getFittedSamples(frequencies);

    Real error = 0.0;
    Complex actual, fitted, diff;

    // Compute the error between the real responses and the fitted ones
    for (size_t i = 0; i < getSamplesSize(); i++) {
        // Sanity check: the response should be on the *same* frequency
        assert(actualSamples[i].first == fittedSamples[i].first);

        // Iterate through all the responses in the vector of each sample
        for (size_t j = 0; j < actualSamples[i].second.size(); j++) {
            // Retrieve the actual and fitted responses
            actual = actualSamples[i].second[j];
            fitted = fittedSamples[i].second[j];

            diff = actual - fitted;

            error += abs(diff * diff);
        }
    }

    return sqrt(error/(getSamplesSize()*getResponseSize()));
}

size_t VectorFitting::getSamplesSize() const {
    return samples_.size();
}

size_t VectorFitting::getResponseSize() const {
    assert(samples_.size() > 0);
    return samples_.front().second.size();
}

size_t VectorFitting::getOrder() const {
    return poles_.size();
}

} /* namespace VectorFitting */

