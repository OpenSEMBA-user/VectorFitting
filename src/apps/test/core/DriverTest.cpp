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
    vector<Sample> noSamples;
    EXPECT_THROW(Driver(noSamples, 3, defaultOptions), runtime_error);
}

TEST_F(DriverTest, simple_case){
    Options opts;
    opts.setPolesType(Options::lincmplx);
    opts.setAsymptoticTrend(Options::linear);

    const size_t Ns = 100;
    const size_t N = 4;
    const std::pair<size_t,size_t> iterations(4,1);

    std::vector<Real> freq = linspace(1,1000, Ns);
    std::vector<Complex> s;
    for (size_t i = 0; i < freq.size(); ++i) {
        s.push_back({0, freq[i]});
    }



    Driver driver(samples, N, opts, {}, iterations);

}