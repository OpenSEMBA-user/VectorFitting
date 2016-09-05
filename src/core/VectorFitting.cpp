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
    poles_ = VectorXcd::Zero(poles.size());
    for (size_t i = 0; i < poles.size(); ++i) {
        poles_(i) = poles[i];
    }
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

    VectorXcd SERD(Nc), SERE(Nc);
    VectorXi SERB(N);
    RowVectorXcd SERA(1,N);
    MatrixXcd  SERC(Nc, N);
    for (size_t i = 0; i < N; ++i) {
        SERA(0,i) = poles_[i];
    }
    VectorXcd roetter;

    // --- Pole identification ---
    if (!options_.isSkipPoleIdentification()) {

        // Finds out which starting poles are complex.
        RowVectorXi cindex = getCIndex(poles_);

        // Builds system - matrix.
        MatrixXcd LAMBD(N, N);
        for (size_t i = 0; i < N; ++i) {
            LAMBD(i,i) = poles_[i];
        }

        MatrixXcd Dk(Ns,N+2);
        MatrixXcd LAMBDprime = LAMBD.transpose().conjugate();
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

        VectorXd x(Nc*(N+1));

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
                        bb(i + n*(N+1)) = Q(2*Ns, N + offs + i)
                                * (Real) Ns * (Real) scale;
                    }
                }
            }  // End of for loop n=1:Nc

            // Computes scaling factor.
            VectorXd Escale(N+1);
            for (size_t col = 0; col < N+1; ++col) {
                Escale(col) = 1.0 / AA.col(col).norm();
                for (size_t i = 0; i < Nc*(N+1); ++i) {
                    AA(i,col) = Escale(col) * AA(i,col);
                }
            }

            x = (AA.transpose() * AA).inverse() * AA.transpose() * bb; // FIXME: Buggy.
            for (size_t i = 0; i < N+1; ++i) {
                x(i) *= Escale(i);
            }

        } // End of if for "relax" flag.

        if (!options_.isRelax()
                || lower  (std::abs(x(0)), toleranceLow_)
                || greater(std::abs(x(N)), toleranceHigh_) ) {
            throw std::runtime_error("Option to do not relax is not implemented");
        }

        VectorXcd C(N);
        for (int i = 0; i < x.rows()-1; ++i) {
            C(i) = x(i);
        }
        for (size_t m = 0; m < N; ++m) {
            if (cindex(m) == 1) {
                const Real r1 = std::real(C(m  ));
                const Real r2 = std::real(C(m+1));
                C(m)   = Complex(r1,  r2);
                C(m+1) = Complex(r1, -r2);
            }
        }
        Real D = x(x.rows()-1);

        // Calculates the zeros for sigma.
        VectorXi B = VectorXi::Ones(N);
        size_t m = 0;
        for (size_t n = 0; n < N; ++n) {
            if (m < N) {
                if (greater(std::abs(LAMBD(m,m)),
                            std::abs(std::real(LAMBD(m,m))))) {
                    LAMBD(m+1,m  ) = - std::imag(LAMBD(m,m));
                    LAMBD(m  ,m+1) =   std::imag(LAMBD(m,m));
                    LAMBD(m  ,m  ) =   std::real(LAMBD(m,m));
                    LAMBD(m+1,m+1) =             LAMBD(m,m);
                    B(m  ) = 2;
                    B(m+1) = 0;
                    const Real aux = std::real(C(m));
                    C(m  ) = std::real(aux);
                    C(m+1) = std::imag(aux);
                    m++;
                }
            }
            m++;
        }

        // Checks LAMBD and C are purely real.
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (!equal(std::imag(LAMBD(i,j)), 0.0)) {
                    throw std::runtime_error("LAMBD is not purely real");
                }
            }
        }
        for (size_t i = 0; i < N; ++i) {
            if (!equal(std::imag(C(i)), 0.0)) {
                throw std::runtime_error("LAMBD is not purely real");
            }
        }

        MatrixXd ZER(N,N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                ZER(i,j) = std::real(LAMBD(i,j)) -  B(i) * std::real(C(i)) / D;
            }
        }

        // Stores roetter
        roetter = EigenSolver<MatrixXd>(ZER, false).eigenvalues();
        if (options_.isStable()) {
            for (size_t i = 0; i < N; ++i) {
                const Real realPart = std::real(roetter(i));
                if (greater(realPart, 0.0)) {
                    roetter(i) = roetter(i) - 2.0 * realPart;
                }
            }
        }

        // TODO Sorterer polene s.a. de reelle kommer first .

        // Stores results for poles.
        SERA = roetter;

    } // End of if for "skip pole identification" flag.

    // --- Residue identification ---
    if (!options_.isSkipResidueIdentification()) {
        // We now calculate SER for f, using the modified zeros of sigma
        // as new poles.
        VectorXcd LAMBD = roetter;
        RowVectorXi cindex = getCIndex(LAMBD);

        // We now calculate the SER for f (new fitting), using the above
        // calculated zeros as known poles.
        MatrixXcd Dk(Ns,N);
        for (size_t m = 0; m < N; ++m) {
            for (size_t i = 0; i < Ns; ++i) {
                if (cindex(m) == 0) {
                    Dk(i,m) = Complex(1,0) / (samples_[i].first - LAMBD(m));
                } else if (cindex(m) == 1) {
                    Dk(i,m)   = Complex(1,0) / (samples_[i].first - LAMBD(m))
                                      + Complex(1,0) / (samples_[i].first - std::conj(LAMBD(m)));
                    Dk(i,m+1) = Complex(0,1) / (samples_[i].first - LAMBD(m))
                                      - Complex(0,1) / (samples_[i].first - std::conj(LAMBD(m)));
                }
            }
        }

        MatrixXcd C(Nc, N);
        VectorXcd BB = VectorXcd::Zero(2*Ns);
        MatrixXcd A;
        switch (options_.getAsymptoticTrend()) {
        case Options::zero:
            A = MatrixXcd::Zero(2*Ns, N);
            break;
        case Options::constant:
            A = MatrixXcd::Zero(2*Ns, N+1);
            break;
        case Options::linear:
            A = MatrixXcd::Zero(2*Ns, N+2);
            break;
        }
        for (size_t n = 0; n < Nc; ++n) {
            for (size_t i = 0; i < Ns; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    A (i    ,j) =   std::real(Dk(i,j)) * weights_(i,n);
                    A (i+Ns ,j) =   std::imag(Dk(i,j)) * weights_(i,n);
                    BB(i)    = std::real(samples_[i].second[n]) * weights_(i,n);
                    BB(i+Ns) = std::imag(samples_[i].second[n]) * weights_(i,n);
                }
            }
            switch (options_.getAsymptoticTrend()) {
            case Options::zero:
                break;
            case Options::constant:
                A.block( 0, N, Ns, 1) = MatrixXcd::Ones(Ns, 1);
                A.block(Ns, N, Ns, 1) = MatrixXcd::Zero(Ns, 1);
                break;
            case Options::linear:
                A.block( 0, N, Ns, 1) = MatrixXcd::Ones(Ns, 1);
                A.block(Ns, N, Ns, 1) = MatrixXcd::Zero(Ns, 1);
                for (size_t i = 0; i < Ns; ++i) {
                    A(i,    N+1) = std::real(samples_[i].first);
                    A(i+Ns, N+1) = std::imag(samples_[i].first);
                }
                break;
            }

            // Computes scaling factor.
            VectorXd Escale(A.cols());
            for (int col = 0; col < A.cols(); ++col) {
                Escale(col) = A.col(col).norm();
                for (int i = 0; i < A.rows(); ++i) {
                    A(i,col) = A(i,col) / Escale(col);
                }
            }

            VectorXcd x = (A.transpose() * A).inverse() * A.transpose() * BB;
            for (int i = 0; i < A.cols(); ++i) {
                x(i) /= Escale(i);
            }

            // Stores results for response;
            for (size_t i = 0; i < N; ++i) {
                C(n,i) = x(i);
            }
            switch (options_.getAsymptoticTrend()) {
            case Options::zero:
                break;
            case Options::constant:
                SERD(n) = x(N);
                break;
            case Options::linear:
                SERD(n) = x(N);
                SERE(n) = x(N+1);
                break;
            }
        } // End of loop over Nc responses.

        for (size_t m = 0; m < N; ++m) {
            if (cindex(m) == 1) {
                for (size_t n = 0; n < Nc; ++n) {
                    const Real r1 = std::real(C(n, m  ));
                    const Real r2 = std::real(C(n, m+1));
                    C(n, m  ) = Complex(r1,  r2);
                    C(n, m+1) = Complex(r1, -r2);
                }
            }
        }

        SERA = LAMBD;
        SERB = VectorXi::Ones(N);
        SERC = C;
    } // End of if for "skip residue identification" flag.

    A_ = MatrixXcd::Zero(N,N);
    for (size_t i = 0; i < N; ++i) {
        A_(i,i) = SERA(i);
    }
    poles_ = SERA.transpose();
    if (!options_.isSkipResidueIdentification()) {
        B_ = SERB;
        C_ = SERC;
        D_ = SERD;
        E_ = SERE;
    } else {
        B_ = VectorXi::Ones(N);
        C_ = MatrixXcd::Zero(Nc, N);
        D_ = VectorXcd::Zero(Nc);
        E_ = VectorXcd::Zero(Nc);
    }
}

// Return the fitted samples: a vector of pairs s <-> f(s), where f(s) is
// computed with the model in (2)
std::vector<Sample> VectorFitting::getFittedSamples() const {
    const size_t N  = getOrder();
    const size_t Ns = getSamplesSize();
    const size_t Nc = getResponseSize();

    MatrixXcd Dk(Ns,N);
    for (size_t m = 0; m < N; ++m) {
        for (size_t i = 0; i < Ns; ++i) {
            Dk(i,m) = Complex(1.0, 0) / (samples_[i].first - A_(m,m));
        }
    }

    std::vector<Sample> res(
            Ns, Sample(Complex(0.0,0.0), std::vector<Complex>(Nc)));
    for (size_t n = 0; n < Nc; ++n) {
        VectorXcd fit(Ns);
        fit = Dk * C_;
        switch (options_.getAsymptoticTrend()) {
        case Options::zero:
            break;
        case Options::constant:
            for (size_t i = 0; i < Ns; ++i) {
                fit(i) += D_(n);
            }
            break;
        case Options::linear:
            for (size_t i = 0; i < Ns; ++i) {
                fit(i) += D_(n) + samples_[i].first * E_(n);
            }
        }
        for (size_t i = 0; i < Ns; ++i) {
            res[i].first = samples_[i].first;
            res[i].second[n] = fit[i];
        }
    }
    return res;
}

std::vector<Complex> VectorFitting::getPoles() {
    std::vector<Complex> res(poles_.rows());
    for (int i = 0; i < poles_.rows(); ++i) {
        res[i] = poles_[i];
    }
    return res;
}

/**
 * Returns the error of the model, measured as the root mean
 * square of the estimated data with respect to the samples.
 * @return Real - Root mean square error of the model.
 */
Real VectorFitting::getRMSE() const {
    std::vector<Sample> fittedSamples = getFittedSamples();

    Real error = 0.0;
    Complex actual, fitted, diff;

    // Compute the error between the real responses and the fitted ones
    for (size_t i = 0; i < getSamplesSize(); i++) {
        // Sanity check: the response should be on the *same* frequency
        assert(samples_[i].first == fittedSamples[i].first);

        // Iterate through all the responses in the vector of each sample
        for (size_t j = 0; j < samples_[i].second.size(); j++) {
            // Retrieve the actual and fitted responses
            actual = samples_[i].second[j];
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
    return (size_t) poles_.rows();
}

RowVectorXi VectorFitting::getCIndex(const VectorXcd& poles) {
    const size_t N = poles.rows();
    RowVectorXi cindex = RowVectorXi::Zero(N);
    for (size_t m = 0; m < N; ++m) {
        if (!equal(std::imag(poles(m)), 0.0)) {
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
    return cindex;
}

} /* namespace VectorFitting */


