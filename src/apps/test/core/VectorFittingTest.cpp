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

TEST_F(MathFittingVectorFittingTest, ctor) {
    Options defaultOptions;
    vector<Sample> noSamples;
    EXPECT_THROW(VectorFitting::VectorFitting(noSamples, 3, defaultOptions),
            runtime_error);
}

// Test first example of Bjorn Gustavsen's code
TEST_F(MathFittingVectorFittingTest, ex1) {
    // Define samples frequencies
    const size_t nS = 101;
    vector<Sample> samples(nS);

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
        samples[k].second = f;
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
    opts.setAsymptoticTrend(Options::linear);
    opts.setSkipPoleIdentification(false);
    opts.setSkipResidueIdentification(false);
    opts.setComplexSpaceState(true);

    VectorFitting::VectorFitting fitting(samples, poles, opts);
    fitting.fit();

    // Compare poles.
    vector<Complex> obtainedPoles = fitting.getPoles();
    vector<Complex> gustavssenPoles = {
            Complex(-5.00000000000118,    0.0           ),
            Complex(-100.000000000017, +499.999999999981),
            Complex(-100.000000000017, -499.999999999981)};
    EXPECT_EQ(gustavssenPoles.size(), obtainedPoles.size());
    for (size_t i = 0; i < gustavssenPoles.size(); ++i) {
        ASSERT_NEAR(gustavssenPoles[i].real(), obtainedPoles[i].real(), 1e-8);
        ASSERT_NEAR(gustavssenPoles[i].imag(), obtainedPoles[i].imag(), 1e-8);
    }

    // Compare residues.
    Matrix<Complex,Dynamic,Dynamic> obtainedResidues = fitting.getC();
    Matrix<Complex,Dynamic,Dynamic> gustavssenResidues;
    gustavssenResidues(0,0)= Complex( 2.0000,   0.0   );
    gustavssenResidues(0,1)= Complex(30.0000, +40.0000);
    gustavssenResidues(0,2)= Complex(30.0000, -40.0000);

    EXPECT_EQ(gustavssenResidues.rows(), obtainedResidues.rows());
    EXPECT_EQ(gustavssenResidues.cols(), obtainedResidues.cols());
    for (int i = 0; i < gustavssenResidues.rows(); ++i) {
        for (int j = 0; j < gustavssenResidues.cols(); ++j) {
            Complex gus = gustavssenResidues(i,j);
            Complex obt = obtainedResidues(i,j);
            ASSERT_NEAR(gus.real(), obt.real(), 1e-8);
            ASSERT_NEAR(gus.imag(), obt.imag(), 1e-8);
        }
    }

    // Compare fitted samples.
    vector<Sample> obtained;
    obtained = fitting.getFittedSamples(s);

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
        EXPECT_FLOAT_EQ(gustavssen[i].real(),
                obtained[i].second[0].real());
        EXPECT_FLOAT_EQ(gustavssen[i].real(),
                obtained[i].second[0].imag());
    }

    // Root Mean Square Error check.
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
    Options defaults;
    VectorFitting::VectorFitting fitting(samples, startingPoles, defaults);
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
    vector<Sample> samples(nS);

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
    Options defaults;
    VectorFitting::VectorFitting fitting(knownResponses, poles, defaults);
    fitting.fit();

    cout << "RMSE: " << fitting.getRMSE() << endl;

    // Error check
    EXPECT_NEAR(0.0, fitting.getRMSE(), 1e-3);
}
