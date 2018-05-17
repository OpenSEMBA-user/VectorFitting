// OpenSEMBAPtions
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


#include <fstream>

#include "gtest/gtest.h"

#include "../../../core/Fitting.h"
#include "SpaceGenerator.h"

using namespace VectorFitting;
using namespace std;

class FittingTest : public ::testing::Test {
protected:
    Real tol_ = 1e-12;
};

TEST_F(FittingTest, ctor) {
    vector<Fitting::Sample> noSamples;

    Options opts;
    opts.setN(3);

    EXPECT_THROW(Fitting(noSamples, opts), runtime_error);
}

// Test first example of Bjorn Gustavsen's code
TEST_F(FittingTest, ex1) {
    // Define samples frequencies
    const size_t nS = 101;
    vector<Fitting::Sample> samples(nS);

    // Compute distribution of the frequencies
    vector<Real> sImag = logspace(pair<Real,Real>(0.0,4.0), nS);
    vector<Complex> s(nS);
    for (size_t k = 0; k < nS; k++) {
        s[k] = Complex(0.0, 2.0 * M_PI * sImag[k]);
    }

    // Populate frequencies and responses
    for (size_t k = 0; k < nS; k++) {
        // Response
        vector<Complex> f(1);
        f[0] =  2.0 /(s[k] + 5.0)
                        + Complex(30.0,40.0)  / (s[k] - Complex(-100.0,500.0))
                        + Complex(30.0,-40.0) / (s[k] - Complex(-100.0,-500.0))
                        + 0.5;

        // Build sample
        samples[k].first = s[k];
        samples[k].second = Fitting::toEigenVector(f);
    }

    // Define starting poles
    const size_t N = 3;
    vector<Complex> poles(N);

    // Compute distribution of the poles
    vector<Real> pReal = logspace(pair<Real,Real>(0.0,4.0), N);

    // Populate starting poles
    for (size_t i = 0; i < N; i++) {
        poles[i] = Complex(-2 * M_PI * pReal[i], 0.0);
    }

    // Model fitting
    Options opts;
    opts.setRelax(true);
    opts.setStable(true);
    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);
    opts.setSkipResidueIdentification(false);

    Fitting fitting(samples, opts, poles);
    fitting.fit();

    // Compare poles.
    vector<Complex> obtainedPoles = fitting.getPoles();
    vector<Complex> gustavssenPoles = {
            Complex(-5.00000000000118,    0.0           ),
            Complex(-100.000000000017, +499.999999999981),
            Complex(-100.000000000017, -499.999999999981)};
    EXPECT_EQ(gustavssenPoles.size(), obtainedPoles.size());
    for (size_t i = 0; i < gustavssenPoles.size(); ++i) {
        EXPECT_NEAR(gustavssenPoles[i].real(), obtainedPoles[i].real(), 1e-6);
        EXPECT_NEAR(gustavssenPoles[i].imag(), obtainedPoles[i].imag(), 1e-6);
    }

    // Compare residues.
    MatrixXcd obtainedResidues = fitting.getC();
    MatrixXcd gustavssenResidues(1,3);
    gustavssenResidues(0,0)= Complex( 2.0000,   0.0   );
    gustavssenResidues(0,1)= Complex(30.0000, +40.0000);
    gustavssenResidues(0,2)= Complex(30.0000, -40.0000);

    EXPECT_EQ(gustavssenResidues.rows(), obtainedResidues.rows());
    EXPECT_EQ(gustavssenResidues.cols(), obtainedResidues.cols());
    for (int i = 0; i < gustavssenResidues.rows(); ++i) {
        for (int j = 0; j < gustavssenResidues.cols(); ++j) {
            Complex gus = gustavssenResidues(i,j);
            Complex obt = obtainedResidues(i,j);
            EXPECT_NEAR(gus.real(), obt.real(), 1e-6);
            EXPECT_NEAR(gus.imag(), obt.imag(), 1e-6);
        }
    }

    // Compare fitted samples.
    vector<Fitting::Sample> obtained = fitting.getFittedSamples();

    vector<Complex> gustavssen = {
            Complex(0.524311630961292, -0.192810298822856),
            Complex(0.507219511076778, -0.187864130127161),
            Complex(0.491073209318934, -0.181598419003598),
            Complex(0.476044935959655, -0.174227918434680),
            Complex(0.462249052829261, -0.165984432707701),
            Complex(0.449743969666221, -0.157100741291502),
            Complex(0.438538305678014, -0.147797320123519),
            Complex(0.428599696241759, -0.138272709027386),
            Complex(0.419864674392760, -0.128697644501629),
            Complex(0.412248359171777, -0.119212544032106),
            Complex(0.405653092416093, -0.109927633001059),
            Complex(0.399975565740413, -0.100924918944686),
            Complex(0.395112301831748, -0.092261276759377),
            Complex(0.390963574797126, -0.083972045091533),
            Complex(0.387435979372693, -0.076074694896263),
            Complex(0.384443909417894, -0.068572281595080),
            Complex(0.381910206794167, -0.061456515289742),
            Complex(0.379766213923313, -0.054710374300441),
            Complex(0.377951422946968, -0.048310248259993),
            Complex(0.376412871775550, -0.042227633671923),
            Complex(0.375104398206323, -0.036430423447639),
            Complex(0.373985830313125, -0.030883838012634),
            Complex(0.373022165096636, -0.025551043434821),
            Complex(0.372182767492424, -0.020393494707570),
            Complex(0.371440607356632, -0.015371031679651),
            Complex(0.370771541985841, -0.010441742041492),
            Complex(0.370153645203925, -0.005561590280341),
            Complex(0.369566580450718, -0.000683792827819),
            Complex(0.368991014397864, +0.004242103938614),
            Complex(0.368408069703299, +0.009271517491601),
            Complex(0.367798821765265, +0.014466613136855),
            Complex(0.367143857416534, +0.019898760309373),
            Complex(0.366422938882318, +0.025651939489784),
            Complex(0.365614865002824, +0.031827527539850),
            Complex(0.364697715963778, +0.038551127804586),
            Complex(0.363649853714853, +0.045982501468674),
            Complex(0.362452426298378, +0.054330309938759),
            Complex(0.361094907869365, +0.063874496523684),
            Complex(0.359586897036525, +0.075001087651863),
            Complex(0.357983195252852, +0.088257639311398),
            Complex(0.356438123537790, +0.104443585510194),
            Complex(0.355327148674715, +0.124759644251374),
            Complex(0.355531558548490, +0.151052583218850),
            Complex(0.359139279135928, +0.186181751872123),
            Complex(0.371247309574021, +0.234342118302268),
            Complex(0.404558809648146, +0.299904526702350),
            Complex(0.488308626271165, +0.376885484580625),
            Complex(0.658441585990537, +0.407020593805587),
            Complex(0.835589930126997, +0.279329873321834),
            Complex(0.851817486962841, +0.083322770582632),
            Complex(0.774482992787240, -0.026996854456600),
            Complex(0.699751975480063, -0.067271822470660),
            Complex(0.645920927746225, -0.077593749471503),
            Complex(0.608733540325068, -0.076787456930570),
            Complex(0.582683923300524, -0.072092104766679),
            Complex(0.563970385934307, -0.066260807152894),
            Complex(0.550188943628670, -0.060352870576803),
            Complex(0.539818368248686, -0.054766437816503),
            Complex(0.531872916338377, -0.049632682142883),
            Complex(0.525695003756959, -0.044973482933641),
            Complex(0.520833207499591, -0.040767094575789),
            Complex(0.516969354871025, -0.036976211184243),
            Complex(0.513873852038642, -0.033560011196527),
            Complex(0.511377561292997, -0.030479169071932),
            Complex(0.509353616355513, -0.027697731143911),
            Complex(0.507705361758002, -0.025183599196943),
            Complex(0.506358162536538, -0.022908418122227),
            Complex(0.505253719298392, -0.020847234119190),
            Complex(0.504346042026707, -0.018978090259822),
            Complex(0.503598545440102, -0.017281632245087),
            Complex(0.502981917735740, -0.015740752869338),
            Complex(0.502472532506207, -0.014340283146353),
            Complex(0.502051248743605, -0.013066728939179),
            Complex(0.501702492599007, -0.011908048329312),
            Complex(0.501413546772093, -0.010853463948596),
            Complex(0.501173995038109, -0.009893304648812),
            Complex(0.500975284179314, -0.009018871500378),
            Complex(0.500810375810801, -0.008222323858795),
            Complex(0.500673467772828, -0.007496581964229),
            Complex(0.500559769877646, -0.006835243184580),
            Complex(0.500465322491454, -0.006232509558518),
            Complex(0.500386849131643, -0.005683124744720),
            Complex(0.500321636257415, -0.005182318847994),
            Complex(0.500267434927955, -0.004725759885255),
            Complex(0.500222380135250, -0.004309510887435),
            Complex(0.500184924485590, -0.003929991818713),
            Complex(0.500153783573737, -0.003583945641766),
            Complex(0.500127890916233, -0.003268407974898),
            Complex(0.500106360721102, -0.002980679880368),
            Complex(0.500088457096683, -0.002718303398059),
            Complex(0.500073568561910, -0.002479039498846),
            Complex(0.500061186928619, -0.002260848180706),
            Complex(0.500050889794465, -0.002061870470320),
            Complex(0.500042326021126, -0.001880412125441),
            Complex(0.500035203683232, -0.001714928860201),
            Complex(0.500029280063802, -0.001564012937968),
            Complex(0.500024353346026, -0.001426380995136),
            Complex(0.500020255711901, -0.001300862975208),
            Complex(0.500016847608251, -0.001186392066092),
            Complex(0.500014012981769, -0.001081995545217),
            Complex(0.500011655318711, -0.000986786447191)
    };

    EXPECT_EQ(gustavssen.size(), obtained.size());
    for (size_t i = 0; i < gustavssen.size(); ++i) {
        EXPECT_NEAR(gustavssen[i].real(), obtained[i].second[0].real(), 1e-6);
        EXPECT_NEAR(gustavssen[i].imag(), obtained[i].second[0].imag(), 1e-6);
    }

    // Root Mean Square Error check.
    EXPECT_NEAR(0.0, fitting.getRMSE(), 1e-8);

    // Get maximum deviation.
    EXPECT_NEAR(0.0, fitting.getMaxDeviation(), 1e-10);
}

TEST_F(FittingTest, ex2){
    // Order of approximation
    const int N = 18;

    // Number of samples
    const int Ns = 100;

    Real D = 0.2;
    Real E = 2e-5;

    vector<Complex> p = {
            Complex(-4.5e+03,   0.0    ),
            Complex(-4.1e+04,   0.0    ),
            Complex(-1.0e+02, - 5.0e+03),
            Complex(-1.0e+02, + 5.0e+03),
            Complex(-1.2e+02, - 1.5e+04),
            Complex(-1.2e+02, + 1.5e+04),
            Complex(-3.0e+03, - 3.5e+04),
            Complex(-3.0e+03, + 3.5e+04),
            Complex(-2.0e+02, - 4.5e+04),
            Complex(-2.0e+02, + 4.5e+04),
            Complex(-1.5e+03, - 4.5e+04),
            Complex(-1.5e+03, + 4.5e+04),
            Complex(-5.0e+02, - 7.0e+04),
            Complex(-5.0e+02, + 7.0e+04),
            Complex(-1.0e+03, - 7.3e+04),
            Complex(-1.0e+03, + 7.3e+04),
            Complex(-2.0e+03, - 9.0e+04),
            Complex(-2.0e+03, + 9.0e+04)};


    vector<Complex> r = {
            Complex(-3.0e+03,   0.0    ),
            Complex(-8.3e+04,   0.0    ),
            Complex(-5.0e+00, - 7.0e+03),
            Complex(-5.0e+00, + 7.0e+03),
            Complex(-2.0e+01, - 1.8e+04),
            Complex(-2.0e+01, + 1.8e+04),
            Complex( 6.0e+03, - 4.5e+04),
            Complex( 6.0e+03, + 4.5e+04),
            Complex( 4.0e+01, - 6.0e+04),
            Complex( 4.0e+01, + 6.0e+04),
            Complex( 9.0e+01, - 1.0e+04),
            Complex( 9.0e+01, + 1.0e+04),
            Complex( 5.0e+04, - 8.0e+04),
            Complex( 5.0e+04, + 8.0e+04),
            Complex( 1.0e+03, - 4.5e+04),
            Complex( 1.0e+03, + 4.5e+04),
            Complex(-5.0e+03, - 9.2e+04),
            Complex(-5.0e+03, + 9.2e+04)};

    for (size_t i = 0; i < N; i++) {
        p[i] *= 2.0 * M_PI;
        r[i] *= 2.0 * M_PI;
    }

    vector<Real> w = linspace(pair<Real, Real>(1, 1e5), Ns);

    // Frequencies
    vector<Complex> s(Ns);

    // Samples
    vector<Fitting::Sample> samples(Ns);

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
        f[0] += s[k] * E + D + 2.0*D;
        f[1] += s[k] * 3.0*E; // Misleading in Ex2 of Gustavssen.

        // Add pair frequency-response
        samples[k] = Fitting::Sample(s[k], Fitting::toEigenVector(f));
    }

    // Starting poles
    vector<Complex> startingPoles;

    // Imaginary parts
    vector<Real> beta = linspace(pair<Real,Real>(w.front(), w.back()), N/2);
    for (size_t n = 0; n < N/2; n++) {
        // Real part
        Real alpha = -beta[n]*1e-2;

        startingPoles.push_back(Complex(alpha, -beta[n]));
        startingPoles.push_back(Complex(alpha,  beta[n]));
    }

    // Model fitting
    Options opts;
    opts.setRelax(true);
    opts.setStable(true);
    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);
    opts.setSkipResidueIdentification(true);

    const size_t nIter = 3;
    Fitting fitting(samples, opts, startingPoles);
    for (size_t iter = 0; iter < nIter; ++iter) {
        if (iter == nIter-1) {
            fitting.options().setSkipResidueIdentification(false);
        }
        fitting.fit();
    }

    // Error check
    EXPECT_NEAR(0.0, fitting.getRMSE(), 1e-8);

    // Get maximum deviation.
    EXPECT_NEAR(0.0, fitting.getMaxDeviation(), 1e-10);

    // Compares with Gustavssen's poles.
    vector<Complex> gusPoles = {
        Complex( - 2.827433388231069e+04, + 0.0                  ),
        Complex( - 2.576105975943579e+05, + 0.0                  ),
        Complex( - 6.283185307179983e+02, + 3.141592653589790e+04),
        Complex( - 6.283185307179983e+02, - 3.141592653589790e+04),
        Complex( - 7.539822368615460e+02, + 9.424777960769372e+04),
        Complex( - 7.539822368615460e+02, - 9.424777960769372e+04),
        Complex( - 1.884955592153879e+04, + 2.199114857512855e+05),
        Complex( - 1.884955592153879e+04, - 2.199114857512855e+05),
        Complex( - 1.256637061435897e+03, + 2.827433388230813e+05),
        Complex( - 1.256637061435897e+03, - 2.827433388230813e+05),
        Complex( - 9.424777960769381e+03, + 2.827433388230815e+05),
        Complex( - 9.424777960769381e+03, - 2.827433388230815e+05),
        Complex( - 3.141592653589785e+03, + 4.398229715025704e+05),
        Complex( - 3.141592653589785e+03, - 4.398229715025704e+05),
        Complex( - 6.283185307179650e+03, + 4.586725274241097e+05),
        Complex( - 6.283185307179650e+03, - 4.586725274241097e+05),
        Complex( - 1.256637061435921e+04, + 5.654866776461634e+05),
        Complex( - 1.256637061435921e+04, - 5.654866776461634e+05)};

    vector<Complex> obtPoles = fitting.getPoles();
    EXPECT_EQ(gusPoles.size(), obtPoles.size());
    for (size_t i = 0; i < gusPoles.size(); ++i) {
        EXPECT_NEAR(gusPoles[i].real(), obtPoles[i].real(), 1e-6);
        EXPECT_NEAR(gusPoles[i].imag(), obtPoles[i].imag(), 1e-6);
    }

}

TEST_F(FittingTest, ex4a){

    // Reads raw data from file.
    ifstream file("testData/fdne.txt");
    EXPECT_TRUE(file.is_open());
    size_t Nc, Ns;
    file >> Nc >> Ns;
    vector<pair<Complex,MatrixXcd>> bigY(Ns);
    for (size_t k = 0; k < Ns; ++k) {
        Real readS;
        file >> readS;
        bigY[k].first = Complex(0.0, readS);
        bigY[k].second = MatrixXcd::Zero(Nc,Nc);
        for (size_t row = 0; row < Nc; ++row) {
            for (size_t col = 0; col < Nc; ++col) {
                Real re, im;
                file >> re >> im;
                bigY[k].second(row,col) = Complex(re,im);
            }
        }
    }
    file.close();

    // Prepares samples. Only first column of bigY is used.
    vector<Fitting::Sample> f;
    for (size_t k = 0; k < Ns; ++k) {
        VectorXcd aux = bigY[k].second.row(0).transpose();
        f.push_back({bigY[k].first, aux});
    }

    // Prepares fitting.
    const size_t N = 50;
    pair<Real,Real> range(f.front().first.imag(),  f.back().first.imag());
    vector<Real> bet = linspace(range, N/2);
    vector<Complex> poles(N); // Starting poles.
    for (size_t n = 0; n < N/2; ++n) {
        poles[2*n  ] = Complex( - bet[n]*1e-2, - bet[n]);
        poles[2*n+1] = Complex( - bet[n]*1e-2, + bet[n]);
    }
    vector<vector<Real>> weights(Ns, vector<Real>(Nc));
    for (size_t i = 0; i < Ns; ++i) {
        for (size_t j = 0; j < Nc; ++j) {
            weights[i][j] = 1.0 / sqrt(std::abs(f[i].second[j]));
        }
    }

    Options opts;
    opts.setRelax(true);
    opts.setStable(true);
    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);
    opts.setSkipPoleIdentification(false);
    opts.setSkipResidueIdentification(false);

    const size_t Niter = 4; // RMS is not reduced after this...
    Fitting fitting(f, opts, poles);
    Real rmse = numeric_limits<Real>::max();
    for (size_t iter = 0; iter < Niter; ++iter) {
        fitting.fit();

        Real newRmse = fitting.getRMSE();
        EXPECT_TRUE(newRmse < rmse);
        rmse = newRmse;

        EXPECT_NEAR(0.0, rmse, 1e-3);
    }


}

// Test Gustavsen's 1999 paper example described in section 4
TEST_F(FittingTest, paperSection4) {
    // Known poles
    vector<Complex> knownPoles;
    knownPoles.push_back(Complex(-  4500,       0));
    knownPoles.push_back(Complex(- 41000,       0));
    knownPoles.push_back(Complex(-   100, +  5000));
    knownPoles.push_back(Complex(-   100, -  5000));
    knownPoles.push_back(Complex(-   120, + 15000));
    knownPoles.push_back(Complex(-   120, - 15000));
    knownPoles.push_back(Complex(-  3000, + 35000));
    knownPoles.push_back(Complex(-  3000, - 35000));
    knownPoles.push_back(Complex(-   200, + 45000));
    knownPoles.push_back(Complex(-   200, - 45000));
    knownPoles.push_back(Complex(-  1500, + 45000));
    knownPoles.push_back(Complex(-  1500, - 45000));
    knownPoles.push_back(Complex(-   500, + 70000));
    knownPoles.push_back(Complex(-   500, - 70000));
    knownPoles.push_back(Complex(-  1000, + 73000));
    knownPoles.push_back(Complex(-  1000, - 73000));
    knownPoles.push_back(Complex(-  2000, + 90000));
    knownPoles.push_back(Complex(-  2000, - 90000));

    // Known residues
    vector<Complex> knownResidues;
    knownResidues.push_back(Complex(-  3000,       0));
    knownResidues.push_back(Complex(- 83000,       0));
    knownResidues.push_back(Complex(-     5, +  7000));
    knownResidues.push_back(Complex(-     5, -  7000));
    knownResidues.push_back(Complex(-    20, + 18000));
    knownResidues.push_back(Complex(-    20, - 18000));
    knownResidues.push_back(Complex(   6000, + 45000));
    knownResidues.push_back(Complex(   6000, - 45000));
    knownResidues.push_back(Complex(     40, + 60000));
    knownResidues.push_back(Complex(     40, - 60000));
    knownResidues.push_back(Complex(     90, + 10000));
    knownResidues.push_back(Complex(     90, - 10000));
    knownResidues.push_back(Complex(  50000, + 80000));
    knownResidues.push_back(Complex(  50000, - 80000));
    knownResidues.push_back(Complex(   1000, + 45000));
    knownResidues.push_back(Complex(   1000, - 45000));
    knownResidues.push_back(Complex(-  5000, + 92000));
    knownResidues.push_back(Complex(-  5000, - 92000));

    // Known parameters
    Real knownD = 0.2;
    Real knownH = 2e-5;

    // Define samples frequencies
    const size_t nS = 100;
    vector<Fitting::Sample> samples(nS);

    // Compute distribution of the frequencies
    vector<Real> sImag = linspace(pair<Real,Real>(1,1e5), nS);

    // Vector to store the fitted samples and loop variable
    vector<Fitting::Sample> knownResponses(nS);

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

        VectorXcd aux = Fitting::toEigenVector(f);
        knownResponses[k] = Fitting::Sample(sk, aux);
    }

    // Define starting poles
    vector<Complex> poles;
    poles.push_back(Complex(-    1e-2, +  1.0));
    poles.push_back(Complex(-    1e-2, -  1.0));
    poles.push_back(Complex(- 1.11e+2, +  1.11e+4));
    poles.push_back(Complex(- 1.11e+2, -  1.11e+4));
    poles.push_back(Complex(- 2.22e+2, +  2.22e+4));
    poles.push_back(Complex(- 2.22e+2, -  2.22e+4));
    poles.push_back(Complex(- 3.33e+2, +  3.33e+4));
    poles.push_back(Complex(- 3.33e+2, -  3.33e+4));
    poles.push_back(Complex(- 4.44e+2, +  4.44e+4));
    poles.push_back(Complex(- 4.44e+2, -  4.44e+4));
    poles.push_back(Complex(- 5.55e+2, +  5.55e+4));
    poles.push_back(Complex(- 5.55e+2, -  5.55e+4));
    poles.push_back(Complex(- 6.66e+2, +  6.66e+4));
    poles.push_back(Complex(- 6.66e+2, -  6.66e+4));
    poles.push_back(Complex(- 7.77e+2, +  7.77e+4));
    poles.push_back(Complex(- 7.77e+2, -  7.77e+4));
    poles.push_back(Complex(- 8.88e+2, +  8.88e+4));
    poles.push_back(Complex(- 8.88e+2, -  8.88e+4));
    poles.push_back(Complex(-    1e+3, +     1e+5));
    poles.push_back(Complex(-    1e+3, -     1e+5));

    // Model fitting
    Options opts;
    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);

    Fitting fitting(knownResponses, opts, poles);
    fitting.fit();

    // Error check
    EXPECT_NEAR(0.0, fitting.getRMSE(), 1e-6);

    // Get maximum deviation.
    EXPECT_NEAR(0.0, fitting.getMaxDeviation(), 1e-8);
}


TEST_F(FittingTest, constant) {
    const size_t Ns = 20;

    std::vector<Fitting::Sample> samples;
    std::vector<VectorXd> weights;
    std::vector<Real> freq = linspace(std::make_pair(1.0, 1000.0), Ns);
    for (size_t i = 0; i < freq.size(); i++) {
        VectorXcd f(3);
        f << 1.0, 2.0, 3.0;
        samples.push_back(Fitting::Sample(Complex(0.0, 2.0*M_PI*freq[i]), f));
        VectorXd weight(1);
        weight << 1.0;
        weights.push_back(weight);
    }



    Options opts;
    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);
    opts.setStable(true);
    opts.setRelax(true);
    opts.setSkipResidueIdentification(false);
    opts.setSkipPoleIdentification(false);

    std::vector<Complex> poles = {
            Complex(-0.000006283185307e3, - 0.006283185307180e3),
            Complex(-0.000006283185307e3, + 0.006283185307180e3),
            Complex(-0.006283185307180e3, - 6.283185307179586e3),
            Complex(-0.006283185307180e3, + 6.283185307179586e3)
    };

    Fitting fitting(samples, opts, poles, weights);
    fitting.fit();

    auto B = fitting.getB();
    for (auto i = 0; i < B.size(); ++i) {
        EXPECT_EQ(1, B(i));
    }

    auto C = fitting.getC();
    for (auto i = 0; i < C.size(); ++i) {
        EXPECT_NEAR(0.0, C(i).real(), tol_);
        EXPECT_NEAR(0.0, C(i).imag(), tol_);
    }

    auto D = fitting.getD();
    EXPECT_EQ(3, D.size());
    EXPECT_NEAR(1.0, D(0).real(), tol_);
    EXPECT_NEAR(2.0, D(1).real(), tol_);
    EXPECT_NEAR(3.0, D(2).real(), tol_);
    EXPECT_NEAR(0.0, D(0).imag(), tol_);
    EXPECT_NEAR(0.0, D(1).imag(), tol_);
    EXPECT_NEAR(0.0, D(2).imag(), tol_);

    auto E = fitting.getE();
    for (auto i = 0; i < E.size(); ++i) {
        EXPECT_NEAR(0.0, E(i).real(), tol_);
        EXPECT_NEAR(0.0, E(i).imag(), tol_);
    }

    EXPECT_NEAR(0.0, fitting.getRMSE(), tol_);

}
