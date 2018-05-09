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

#include <fstream>

#include "gtest/gtest.h"

#include "Driver.h"
#include "SpaceGenerator.h"

using namespace VectorFitting;
using namespace std;

class DriverTest : public ::testing::Test {};

TEST_F(DriverTest, ctor) {
    Options defaultOptions;
    defaultOptions.setN(3);
    defaultOptions.setIterations({4,1});

    vector<Driver::Sample> noSamples;

    EXPECT_THROW(Driver(noSamples, defaultOptions), runtime_error);
}

TEST_F(DriverTest, simple_case){

    const size_t Ns = 20;
    const size_t N = 4;
    const size_t Nc = 2;

    Options opts;
    opts.setPolesType(Options::PolesType::lincmplx);
    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);
    opts.setN(4);
    opts.setIterations({4, 1});

    std::vector<Driver::Sample> samples;
    {
        std::vector<Real> freq =
                linspace(std::make_pair(2*M_PI*1.0, 2*M_PI*1000.0), Ns);
        for (size_t i = 0; i < freq.size(); ++i) {
            Complex s = {0, freq[i]};
            MatrixXcd data(Nc,Nc);
            data << Complex(1.0, 0.0), Complex(2.0, 0.0),
                    Complex(2.0, 0.0), Complex(3.0, 0.0);
            samples.push_back({s, data});
        }
    }

    Driver driver(samples, opts);

    // Checks sizes of returned fitting parameters.
    EXPECT_EQ(Nc*N, driver.getA().rows());
    EXPECT_EQ(Nc*N, driver.getA().cols());

    EXPECT_EQ(Nc*N, driver.getB().rows());
    EXPECT_EQ(Nc,   driver.getB().cols());

    EXPECT_EQ(Nc,   driver.getC().rows());
    EXPECT_EQ(Nc*N, driver.getC().cols());

    EXPECT_EQ(Nc, driver.getE().rows());
    EXPECT_EQ(Nc, driver.getE().cols());

    EXPECT_EQ(N, driver.getR().size());
    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(Nc, driver.getR()[i].rows());
        EXPECT_EQ(Nc, driver.getR()[i].cols());
    }

    // Checks values of returned parameters.
    // TODO
}
