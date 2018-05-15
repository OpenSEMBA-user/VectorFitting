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

	size_t Ns, N, Nc;
	Options opts;
    std::vector<Driver::Sample> samples;

	void SetUp() {
	    Ns = 20;
	    N = 4;
	    Nc = 2;

	    opts.setPolesType(Options::PolesType::lincmplx);
	    opts.setAsymptoticTrend(Options::AsymptoticTrend::linear);
	    opts.setN(4);
	    opts.setIterations({0, 1});

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
	}
};

TEST_F(DriverTest, ctor) {
	Options defaultOptions;
	defaultOptions.setN(3);
	defaultOptions.setIterations({4,1});

	vector<Driver::Sample> noSamples;

	EXPECT_THROW(Driver(noSamples, defaultOptions), runtime_error);
}

TEST_F(DriverTest, simple_case){
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

}

TEST_F(DriverTest, Initial_poles_lincmplx) {
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

TEST_F(DriverTest, Matrix_A) {
	Driver driver(samples, opts);
	// Checks values of returned parameters.
	EXPECT_FLOAT_EQ(-0.2223e3, driver.getA()(0,0).real());
	EXPECT_FLOAT_EQ(-0.4390e3, driver.getA()(1,1).real());
	EXPECT_FLOAT_EQ(-0.0420e3, driver.getA()(2,2).real());
	EXPECT_FLOAT_EQ(-0.0420e3, driver.getA()(3,3).real());
	EXPECT_FLOAT_EQ(-0.2223e3, driver.getA()(4,4).real());
	EXPECT_FLOAT_EQ(-0.4390e3, driver.getA()(5,5).real());
	EXPECT_FLOAT_EQ(-0.0420e3, driver.getA()(6,6).real());
	EXPECT_FLOAT_EQ(-0.0420e3, driver.getA()(7,7).real());

	EXPECT_FLOAT_EQ( 0.0000e3, driver.getA()(0,0).imag());
	EXPECT_FLOAT_EQ( 0.0000e3, driver.getA()(1,1).imag());
	EXPECT_FLOAT_EQ( 6.2730e3, driver.getA()(2,2).imag());
	EXPECT_FLOAT_EQ(-6.2730e3, driver.getA()(3,3).imag());
	EXPECT_FLOAT_EQ( 0.0000e3, driver.getA()(4,4).imag());
	EXPECT_FLOAT_EQ( 0.0000e3, driver.getA()(5,5).imag());
	EXPECT_FLOAT_EQ( 6.2730e3, driver.getA()(6,6).imag());
	EXPECT_FLOAT_EQ(-6.2730e3, driver.getA()(7,7).imag());

	for (auto i = 0; i < driver.getA().rows(); ++i){
		for (auto j = 0; j < driver.getA().cols(); ++j){
			if (i != j){
				EXPECT_FLOAT_EQ(0.0, driver.getA()(i,j).real());
				EXPECT_FLOAT_EQ(0.0, driver.getA()(i,j).imag());
			}
		}
	}
}

TEST_F(DriverTest,Matrix_B){

	Driver driver(samples, opts);

	EXPECT_FLOAT_EQ(1,driver.getB()(0,0));
	EXPECT_FLOAT_EQ(1,driver.getB()(1,0));
	EXPECT_FLOAT_EQ(1,driver.getB()(2,0));
	EXPECT_FLOAT_EQ(1,driver.getB()(3,0));
	EXPECT_FLOAT_EQ(0,driver.getB()(4,0));
	EXPECT_FLOAT_EQ(0,driver.getB()(5,0));
	EXPECT_FLOAT_EQ(0,driver.getB()(6,0));
	EXPECT_FLOAT_EQ(0,driver.getB()(7,0));

	EXPECT_FLOAT_EQ(0,driver.getB()(0,1));
	EXPECT_FLOAT_EQ(0,driver.getB()(1,1));
	EXPECT_FLOAT_EQ(0,driver.getB()(2,1));
	EXPECT_FLOAT_EQ(0,driver.getB()(3,1));
	EXPECT_FLOAT_EQ(1,driver.getB()(4,1));
	EXPECT_FLOAT_EQ(1,driver.getB()(5,1));
	EXPECT_FLOAT_EQ(1,driver.getB()(6,1));
	EXPECT_FLOAT_EQ(1,driver.getB()(7,1));
}

TEST_F(DriverTest,Matrix_C){

	Driver driver(samples, opts);

	EXPECT_NEAR(0.0, driver.getC()(0,0).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,1).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,2).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,3).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,4).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,5).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,6).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,7).real(), tol_);

	EXPECT_NEAR(0.0, driver.getC()(0,0).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,1).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,2).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,3).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,4).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,5).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,6).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(0,7).imag(), tol_);

	EXPECT_NEAR(0.0, driver.getC()(1,0).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,1).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,2).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,3).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,4).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,5).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,6).real(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,7).real(), tol_);

	EXPECT_NEAR(0.0, driver.getC()(1,0).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,1).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,2).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,3).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,4).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,5).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,6).imag(), tol_);
	EXPECT_NEAR(0.0, driver.getC()(1,7).imag(), tol_);

}

TEST_F(DriverTest,Matrix_D){

	Driver driver(samples, opts);

	EXPECT_FLOAT_EQ(1.0,driver.getD()(0,0).real());
	EXPECT_FLOAT_EQ(2.0,driver.getD()(0,1).real());
	EXPECT_FLOAT_EQ(2.0,driver.getD()(1,0).real());
	EXPECT_FLOAT_EQ(3.0,driver.getD()(1,1).real());

	EXPECT_FLOAT_EQ(0.0,driver.getD()(0,0).imag());
	EXPECT_FLOAT_EQ(0.0,driver.getD()(0,1).imag());
	EXPECT_FLOAT_EQ(0.0,driver.getD()(1,0).imag());
	EXPECT_FLOAT_EQ(0.0,driver.getD()(1,1).imag());

}



TEST_F(DriverTest,Matrix_E){

	Driver driver(samples, opts);

	EXPECT_NEAR(0.0, driver.getE()(0,0).real(), tol_);
	EXPECT_NEAR(0.0, driver.getE()(0,1).real(), tol_);
	EXPECT_NEAR(0.0, driver.getE()(1,0).real(), tol_);
	EXPECT_NEAR(0.0, driver.getE()(1,1).real(), tol_);
}



TEST_F(DriverTest,Matrix_R){

	Driver driver(samples, opts);
	auto R = driver.ss2pr().second;

	EXPECT_NEAR(0.0, R[0](0,0).real(), tol_);
	EXPECT_NEAR(0.0, R[0](0,1).real(), tol_);
	EXPECT_NEAR(0.0, R[0](1,0).real(), tol_);
	EXPECT_NEAR(0.0, R[0](1,1).real(), tol_);
	EXPECT_NEAR(0.0, R[1](0,0).real(), tol_);
	EXPECT_NEAR(0.0, R[1](0,1).real(), tol_);
	EXPECT_NEAR(0.0, R[1](1,0).real(), tol_);
	EXPECT_NEAR(0.0, R[1](1,1).real(), tol_);
	EXPECT_NEAR(0.0, R[2](0,0).real(), tol_);
	EXPECT_NEAR(0.0, R[2](0,1).real(), tol_);
	EXPECT_NEAR(0.0, R[2](1,0).real(), tol_);
	EXPECT_NEAR(0.0, R[2](1,1).real(), tol_);
	EXPECT_NEAR(0.0, R[3](0,0).real(), tol_);
	EXPECT_NEAR(0.0, R[3](0,1).real(), tol_);
	EXPECT_NEAR(0.0, R[3](1,0).real(), tol_);
	EXPECT_NEAR(0.0, R[3](1,1).real(), tol_);
	EXPECT_NEAR(0.0, R[0](0,0).imag(), tol_);
	EXPECT_NEAR(0.0, R[0](0,1).imag(), tol_);
	EXPECT_NEAR(0.0, R[0](1,0).imag(), tol_);
	EXPECT_NEAR(0.0, R[0](1,1).imag(), tol_);
	EXPECT_NEAR(0.0, R[1](0,0).imag(), tol_);
	EXPECT_NEAR(0.0, R[1](0,1).imag(), tol_);
	EXPECT_NEAR(0.0, R[1](1,0).imag(), tol_);
	EXPECT_NEAR(0.0, R[1](1,1).imag(), tol_);
	EXPECT_NEAR(0.0, R[2](0,0).imag(), tol_);
	EXPECT_NEAR(0.0, R[2](0,1).real(), tol_);
	EXPECT_NEAR(0.0, R[2](1,0).real(), tol_);
	EXPECT_NEAR(0.0, R[2](1,1).real(), tol_);
	EXPECT_NEAR(0.0, R[3](0,0).real(), tol_);
	EXPECT_NEAR(0.0, R[3](0,1).real(), tol_);
	EXPECT_NEAR(0.0, R[3](1,1).real(), tol_);
	EXPECT_NEAR(0.0, R[3](1,0).real(), tol_);

}

TEST_F(DriverTest,Matrix_poles){

	Driver driver(samples, opts);

	auto poles = driver.ss2pr().first;

	EXPECT_FLOAT_EQ(-2.222719257260298e+02, poles[0].real());
	EXPECT_FLOAT_EQ(-4.390094803597557e+02, poles[1].real());
	EXPECT_FLOAT_EQ(-4.201729969020380e+01, poles[2].real());
	EXPECT_FLOAT_EQ(-4.201729969020380e+01, poles[3].real());

	EXPECT_FLOAT_EQ(0.0,                   poles[0].imag());
	EXPECT_FLOAT_EQ(0.0,                   poles[1].imag());
	EXPECT_FLOAT_EQ( 6.273049560187480e+03, poles[2].imag());
	EXPECT_FLOAT_EQ(-6.273049560187480e+03,poles[3].imag());
}



