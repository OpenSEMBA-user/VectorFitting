// OpenSEMBA
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

#include "Driver.h"
#include "SpaceGenerator.h"



using namespace VectorFitting;
using namespace std;

class DriverTest : public ::testing::Test {
protected:
	const double tol_ = 1.5e-5;
};

TEST_F(DriverTest, ctor) {
	Options defaultOptions;
	defaultOptions.setN(3);
	defaultOptions.setIterations({4,1});

	vector<Driver::Sample> noSamples;

	EXPECT_THROW(Driver(noSamples, defaultOptions), runtime_error);
}

TEST_F(DriverTest, initial_poles_lincmplx) {
    std::pair<Real,Real> range(2*M_PI*1.0, 2*M_PI*1000.0);
    Options opts;
    opts.setN(4);
    opts.setPolesType(Options::PolesType::lincmplx);

    std::vector<Complex> poles = Driver::buildPoles(range, opts);

    EXPECT_FLOAT_EQ(-0.000006283185307e3, poles[0].real());
    EXPECT_FLOAT_EQ(-0.000006283185307e3, poles[1].real());
    EXPECT_FLOAT_EQ(-0.006283185307180e3, poles[2].real());
    EXPECT_FLOAT_EQ(-0.006283185307180e3, poles[3].real());

    EXPECT_FLOAT_EQ(-0.006283185307180e3, poles[0].imag());
    EXPECT_FLOAT_EQ(+0.006283185307180e3, poles[1].imag());
    EXPECT_FLOAT_EQ(-6.283185307179586e3, poles[2].imag());
    EXPECT_FLOAT_EQ(+6.283185307179586e3, poles[3].imag());
}

TEST_F(DriverTest, constant){
    size_t Ns = 20;
    size_t N = 4;
    size_t Nc = 2;

    Options opts;
    opts.setPolesType(Options::PolesType::lincmplx);
    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);
    opts.setN(N);
    opts.setIterations({0, 1});

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

	EXPECT_EQ(N, driver.ss2pr().second.size());
	for (size_t i = 0; i < N; ++i) {
		EXPECT_EQ(Nc, driver.ss2pr().second[i].rows());
		EXPECT_EQ(Nc, driver.ss2pr().second[i].cols());
	}


	EXPECT_EQ(1, driver.getB()(0,0));
	EXPECT_EQ(1, driver.getB()(1,0));
	EXPECT_EQ(1, driver.getB()(2,0));
	EXPECT_EQ(1, driver.getB()(3,0));
	EXPECT_EQ(0, driver.getB()(4,0));
	EXPECT_EQ(0, driver.getB()(5,0));
	EXPECT_EQ(0, driver.getB()(6,0));
	EXPECT_EQ(0, driver.getB()(7,0));

	EXPECT_EQ(0, driver.getB()(0,1));
	EXPECT_EQ(0, driver.getB()(1,1));
	EXPECT_EQ(0, driver.getB()(2,1));
	EXPECT_EQ(0, driver.getB()(3,1));
	EXPECT_EQ(1, driver.getB()(4,1));
	EXPECT_EQ(1, driver.getB()(5,1));
	EXPECT_EQ(1, driver.getB()(6,1));
	EXPECT_EQ(1, driver.getB()(7,1));

	for (MatrixXcd::Index i = 0; i < driver.getC().size(); ++i) {
	    EXPECT_NEAR(0.0, driver.getC()(i).real(), tol_);
	    EXPECT_NEAR(0.0, driver.getC()(i).imag(), tol_);
	}

	EXPECT_FLOAT_EQ(1.0,driver.getD()(0,0).real());
	EXPECT_FLOAT_EQ(2.0,driver.getD()(0,1).real());
	EXPECT_FLOAT_EQ(2.0,driver.getD()(1,0).real());
	EXPECT_FLOAT_EQ(3.0,driver.getD()(1,1).real());

	EXPECT_FLOAT_EQ(0.0,driver.getD()(0,0).imag());
	EXPECT_FLOAT_EQ(0.0,driver.getD()(0,1).imag());
	EXPECT_FLOAT_EQ(0.0,driver.getD()(1,0).imag());
	EXPECT_FLOAT_EQ(0.0,driver.getD()(1,1).imag());


    for (MatrixXcd::Index i = 0; i < driver.getE().size(); ++i) {
        EXPECT_NEAR(0.0, driver.getE()(i).real(), tol_);
        EXPECT_NEAR(0.0, driver.getE()(i).imag(), tol_);
    }

	auto R = driver.ss2pr().second;
	for (size_t p = 0; p < R.size(); ++p) {
	    EXPECT_EQ(Nc, R[p].rows());
	    EXPECT_EQ(Nc, R[p].cols());
	    for (size_t i = 0; i < Nc; ++i) {
	        for (size_t j = 0; j < Nc; ++j) {
	            EXPECT_NEAR(0.0, R[p](i,j).real(), tol_);
	            EXPECT_NEAR(0.0, R[p](i,j).imag(), tol_);
	        }
	    }
	}


	EXPECT_NEAR(0.0, driver.getRMSE(), tol_);

}



