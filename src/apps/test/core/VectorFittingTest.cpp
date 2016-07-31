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

#include "gtest/gtest.h"
#include "VectorFitting.h"
#include "SpaceGenerator.h"

using namespace VectorFitting;
using namespace std;

class MathFittingVectorFittingTest : public ::testing::Test {

};

// Test first example of Bjorn Gustavsen's code
TEST_F(MathFittingVectorFittingTest, ex1) {
    // Define samples frequencies
    const size_t nS = 101;
    vector<VectorFitting::Sample> samples(nS);

    // Compute distribution of the frequencies
    vector<Real> sImag = logspace(pair<Real,Real>(0.0,4.0), nS);


    cout << "SAMPLES" << endl;

    // Populate frequencies and responses
    for (size_t k = 0; k < nS; k++) {
        // Frequency
        const Complex sk = Complex(0.0, 2.0 * M_PI * sImag[k]);

        // Response
        vector<Complex> f(1);
        f[0] =  2.0 /(sk+5.0)
                + Complex(30.0,40.0)  / (sk-Complex(-100.0,500.0))
                + Complex(30.0,-40.0) / (sk-Complex(-100.0,-500.0))
                + 0.5;

        // Build sample
        samples[k].first = sk;
        samples[k].second = f;

        cout << samples[k].first << endl;
    }

    // Define starting poles
    const size_t N = 3;
    vector<Complex> poles(N);

    // Compute distribution of the poles
    vector<Real> pReal = logspace(pair<Real,Real>(0.0,4.0), N);

    cout << "STARTING" << endl;
    // Populate starting poles
    for (size_t i = 0; i < N; i++) {
        poles[i] = Complex(-2 * M_PI * pReal[i], 0.0);
        cout << poles[i] << endl;
    }

    // Model fitting
    VectorFitting::VectorFitting fitting(samples, poles, N);
    fitting.fit();

    // Error check
    EXPECT_NEAR(0.0, fitting.getRMSE(), 1e-3);
}


TEST_F(MathFittingVectorFittingTest, ex2){
    // Order of approximation
    const int N = 18;

    // Number of samples
    const int Ns = 100;

    Real D = 0.2;
    Real E = 2e-5;

    vector<Complex> p;
    p.push_back(Complex(-4500.0, +0.0));
    p.push_back(Complex(-41000.0, +0.0));
    p.push_back(Complex(-100.0, -5000.0));
    p.push_back(Complex(-100.0, +5000.0));
    p.push_back(Complex(-120.0, -15000.0));
    p.push_back(Complex(-120.0, +15000.0));
    p.push_back(Complex(-3000.0, -35000.0));
    p.push_back(Complex(-3000.0, +35000.0));
    p.push_back(Complex(-200.0, -45000.0));
    p.push_back(Complex(-200.0, +45000.0));
    p.push_back(Complex(-1500.0, -45000.0));
    p.push_back(Complex(-1500.0, +45000.0));
    p.push_back(Complex(-500.0, -70000.0));
    p.push_back(Complex(-500.0, +70000.0));
    p.push_back(Complex(-1000.0, -73000.0));
    p.push_back(Complex(-1000.0, +73000.0));
    p.push_back(Complex(-2000.0, -90000.0));
    p.push_back(Complex(-2000.0, +90000.0));

    vector<Complex> r;
    r.push_back(Complex(-3000.0, +0.0));
    r.push_back(Complex(-83000.0, +0.0));
    r.push_back(Complex(-5.0, -7000.0));
    r.push_back(Complex(-5.0, +7000.0));
    r.push_back(Complex(-20.0, -18000.0));
    r.push_back(Complex(-20.0, +18000.0));
    r.push_back(Complex(6000.0, -45000.0));
    r.push_back(Complex(6000.0, +45000.0));
    r.push_back(Complex(40.0, -60000.0));
    r.push_back(Complex(40.0, +60000.0));
    r.push_back(Complex(90.0, -10000.0));
    r.push_back(Complex(90.0, +10000.0));
    r.push_back(Complex(50000.0, -80000.0));
    r.push_back(Complex(50000.0, +80000.0));
    r.push_back(Complex(1000.0, -45000.0));
    r.push_back(Complex(1000.0, +45000.0));
    r.push_back(Complex(-5000.0, -92000.0));
    r.push_back(Complex(-5000.0, +92000.0));

    for (size_t i = 0; i < N; i++) {
        p[i] *= 2.0 * M_PI;
        r[i] *= 2.0 * M_PI;
    }

    vector<Real> w = linspace(pair<Real, Real>(1, 1e5), Ns);

    // Frequencies
    vector<Complex> s(Ns);

    // Samples
    vector<Sample> samples(Ns);

    // Build samples
    for (size_t k = 0; k < Ns; k++) {
        // Build frequencies
        w[k] *= 2.0 * M_PI;
        s[k] = Complex(0.0, w[k]);

        // Initialize responses to zero
        vector<Complex> f(2, 0.0);

        // Compute responses with the model
        for (size_t n = 0; n < N; n++) {
            if(n < 10)
                f[0] += r[n] / (s[k] - p[n]);
            if(n >= 8)
                f[1] += r[n] / (s[k] - p[n]);
        }
        f[0] += s[k] * E + D;
        f[1] += s[k] * 3.0*E + 2.0*D;

        // Add pair frequency-response
        samples[k] = Sample(s[k], f);
    }

    // Starting poles
    vector<Complex> startingPoles;

    // Imaginary parts
    vector<Real> beta = linspace(pair<Real,Real>(w.front(), w.back()), N/2);
    for (size_t n = 0; n < N/2; n++) {
        // Real part
        Real alpha = -beta[n]*1e-2;

        startingPoles.push_back(Complex(alpha, -beta[n]));
        startingPoles.push_back(Complex(alpha, beta[n]));
    }

    // Model fitting
    VectorFitting::VectorFitting fitting(samples, startingPoles, N);
    fitting.fit();

    // Error check
    EXPECT_NEAR(0.0, fitting.getRMSE(), 1e-3);


}

// Test paper's example described in section 4
TEST_F(MathFittingVectorFittingTest, test) {
    // Known poles
    vector<Complex> knownPoles;
    knownPoles.push_back(Complex(-4500, 0.0));
    knownPoles.push_back(Complex(-41000, 0.0));
    knownPoles.push_back(Complex(-100, +5000));
    knownPoles.push_back(Complex(-100, -5000));
    knownPoles.push_back(Complex(-120, +15000));
    knownPoles.push_back(Complex(-120, -15000));
    knownPoles.push_back(Complex(-3000, +35000));
    knownPoles.push_back(Complex(-3000, -35000));
    knownPoles.push_back(Complex(-200, +45000));
    knownPoles.push_back(Complex(-200, -45000));
    knownPoles.push_back(Complex(-1500, +45000));
    knownPoles.push_back(Complex(-1500, -45000));
    knownPoles.push_back(Complex(-500, +70000));
    knownPoles.push_back(Complex(-500, -70000));
    knownPoles.push_back(Complex(-1000, +73000));
    knownPoles.push_back(Complex(-1000, -73000));
    knownPoles.push_back(Complex(-2000, +90000));
    knownPoles.push_back(Complex(-2000,-90000));

    // Known residues
    vector<Complex> knownResidues;
    knownResidues.push_back(Complex(-3000, 0.0));
    knownResidues.push_back(Complex(-83000, 0.0));
    knownResidues.push_back(Complex(-5, +7000));
    knownResidues.push_back(Complex(-5, -7000));
    knownResidues.push_back(Complex(-20, +18000));
    knownResidues.push_back(Complex(-20, -18000));
    knownResidues.push_back(Complex(6000, +45000));
    knownResidues.push_back(Complex(6000, -45000));
    knownPoles.push_back(Complex(40, +60000));
    knownPoles.push_back(Complex(40, -60000));
    knownPoles.push_back(Complex(90, +10000));
    knownPoles.push_back(Complex(90, -10000));
    knownPoles.push_back(Complex(50000, +80000));
    knownPoles.push_back(Complex(50000, -80000));
    knownPoles.push_back(Complex(1000, +45000));
    knownPoles.push_back(Complex(1000, -45000));
    knownPoles.push_back(Complex(-5000, +92000));
    knownPoles.push_back(Complex(-5000, -92000));

    // Known parameters
    Real knownD = 0.2;
    Real knownH = 2e-5;

    // Define samples frequencies
    const size_t nS = 100;
    vector<VectorFitting::Sample> samples(nS);

    // Compute distribution of the frequencies
    vector<Real> sImag = linspace(pair<Real,Real>(1,1e5), nS);

    // Vector to store the fitted samples and loop variable
    vector<Sample> knownResponses(nS);

    for (size_t k = 0; k < nS; k++) {
        // Independent variable s
        Complex sk = Complex(0.0, 2 * M_PI * sImag[k]);

        // Computation of the dependent variable f(s) with the model (see (2))
        vector<Complex> f(1);
        f[0] = Complex(0,0);

        for (size_t n = 0; n < knownPoles.size(); n++) {
            Complex cn = knownResidues[n];
            Complex an = knownPoles[n];

            f[0] += cn / (sk - an);
        }

        f[0] += knownD + sk * knownH;

        knownResponses[k] = Sample(sk, f);
    }


    // Define starting poles
    const size_t N = 20;
    vector<Complex> poles;

    poles.push_back(Complex(-1e-2, +1.0));
    poles.push_back(Complex(-1e-2, -1.0));
    poles.push_back(Complex(-1.11e2, +1.11e4));
    poles.push_back(Complex(-1.11e2, -1.11e4));
    poles.push_back(Complex(-2.22e2, +2.22e4));
    poles.push_back(Complex(-2.22e2, -2.22e4));
    poles.push_back(Complex(-3.33e2, +3.33e4));
    poles.push_back(Complex(-3.33e2, -3.33e4));
    poles.push_back(Complex(-4.44e2, +4.44e4));
    poles.push_back(Complex(-4.44e2, -4.44e4));
    poles.push_back(Complex(-5.55e2, +5.55e4));
    poles.push_back(Complex(-5.55e2, -5.55e4));
    poles.push_back(Complex(-6.66e2, +6.66e4));
    poles.push_back(Complex(-6.66e2, -6.66e4));
    poles.push_back(Complex(-7.77e2, +7.77e4));
    poles.push_back(Complex(-7.77e2, -7.77e4));
    poles.push_back(Complex(-8.88e2, +8.88e4));
    poles.push_back(Complex(-8.88e2, -8.88e4));
    poles.push_back(Complex(-1e3, +1e5));
    poles.push_back(Complex(-1e3, -1e5));

    // Model fitting
    VectorFitting::VectorFitting fitting(knownResponses, poles, N);
    fitting.fit();

    cout << "RMSE: " << fitting.getRMSE() << endl;

    // Error check
    EXPECT_NEAR(0.0, fitting.getRMSE(), 1e-3);
}
