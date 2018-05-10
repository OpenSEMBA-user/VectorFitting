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

}


TEST_F(DriverTest,Matrix_A){

	  	const size_t Ns = 20;
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
	    // Checks values of returned parameters.
	    EXPECT_FLOAT_EQ(1.0e+02 * 0.0291,driver.getA()(0,0).real());
	    EXPECT_FLOAT_EQ(1.0e+02 * 0.1855,driver.getA()(1,1).real());
        EXPECT_FLOAT_EQ(1.0e+02 * 2.3410,driver.getA()(2,2).real());
        EXPECT_FLOAT_EQ(1.0e+02 * 2.3410,driver.getA()(3,3).real());
        EXPECT_FLOAT_EQ(1.0e+02 * 0.0291,driver.getA()(4,4).real());
        EXPECT_FLOAT_EQ(1.0e+02 * 0.1855,driver.getA()(5,5).real());
        EXPECT_FLOAT_EQ(1.0e+02 * 2.3410,driver.getA()(6,6).real());
        EXPECT_FLOAT_EQ(1.0e+02 * 2.3410,driver.getA()(7,7).real());

        EXPECT_FLOAT_EQ(1.0e+02 * .0,driver.getA()(0,0).imag());
        EXPECT_FLOAT_EQ(1.0e+02 * .0,driver.getA()(1,1).imag());
        EXPECT_FLOAT_EQ(1.0e+02 * .3406,driver.getA()(2,2).imag());
        EXPECT_FLOAT_EQ(1.0e+02 * .3406,driver.getA()(3,3).imag());
        EXPECT_FLOAT_EQ(1.0e+02 * .0,driver.getA()(4,4).imag());
        EXPECT_FLOAT_EQ(1.0e+02 * .0,driver.getA()(5,5).imag());
        EXPECT_FLOAT_EQ(1.0e+02 * .3406,driver.getA()(6,6).imag());
        EXPECT_FLOAT_EQ(1.0e+02 * .3406,driver.getA()(7,7).imag());

        for (size_t i = 0; i < driver.getA().rows(); ++i){
        	for (size_t j = 0; j < driver.getA().cols(); ++j){
        		if (i != j){
        			EXPECT_FLOAT_EQ(0.0, driver.getA()(i,j).real());
        			EXPECT_FLOAT_EQ(0.0, driver.getA()(i,j).imag());
        		}
        	}
        }

  }

TEST_F(DriverTest,Matrix_B){

	  	const size_t Ns = 20;
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

	  	const size_t Ns = 20;
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



	    EXPECT_FLOAT_EQ(2.03407406382732e-17,driver.getC()(0,0).real());
	    EXPECT_FLOAT_EQ(7.461084540623050e-17,driver.getC()(0,1).real());
	    EXPECT_FLOAT_EQ(-4.141748673186497e-14,driver.getC()(0,2).real());
	    EXPECT_FLOAT_EQ(-4.141748673186497e-14,driver.getC()(0,3).real());
	    EXPECT_FLOAT_EQ(4.068148127654645e-17,driver.getC()(0,4).real());
	    EXPECT_FLOAT_EQ(1.492216908124610e-16,driver.getC()(0,5).real());
	    EXPECT_FLOAT_EQ(-8.283497346372995e-14,driver.getC()(0,6).real());
	    EXPECT_FLOAT_EQ(-8.283497346372995e-14,driver.getC()(0,7).real());

	    EXPECT_FLOAT_EQ(0.0,driver.getC()(0,0).imag());
	    EXPECT_FLOAT_EQ(0.0,driver.getC()(0,1).imag());
	    EXPECT_FLOAT_EQ(1.113383197391269e-13,driver.getC()(0,2).imag());
	    EXPECT_FLOAT_EQ(-1.113383197391269e-13,driver.getC()(0,3).imag());
	    EXPECT_FLOAT_EQ(0.0,driver.getC()(0,4).imag());
	    EXPECT_FLOAT_EQ(0.0,driver.getC()(0,5).imag());
	    EXPECT_FLOAT_EQ(2.226766394782538e-13,driver.getC()(0,6).imag());
	    EXPECT_FLOAT_EQ(-2.226766394782538e-13,driver.getC()(0,7).imag());

	    EXPECT_FLOAT_EQ(4.068148127654645e-17,driver.getC()(1,0).real());
	    EXPECT_FLOAT_EQ(1.492216908124610e-16,driver.getC()(1,1).real());
        EXPECT_FLOAT_EQ(-8.283497346372995e-14,driver.getC()(1,2).real());
        EXPECT_FLOAT_EQ(-8.283497346372995e-14,driver.getC()(1,3).real());
        EXPECT_FLOAT_EQ(2.771204309737382e-16,driver.getC()(1,4).real());
        EXPECT_FLOAT_EQ(-1.007803076036114e-14,driver.getC()(1,5).real());
        EXPECT_FLOAT_EQ(-3.099535157216429e-13,driver.getC()(1,6).real());
        EXPECT_FLOAT_EQ(-3.099535157216429e-13,driver.getC()(1,7).real());

        EXPECT_FLOAT_EQ(0.0,driver.getC()(1,0).imag());
        EXPECT_FLOAT_EQ(0.0,driver.getC()(1,1).imag());
        EXPECT_FLOAT_EQ(2.226766394782538e-13,driver.getC()(1,2).imag());
        EXPECT_FLOAT_EQ(-2.226766394782538e-13,driver.getC()(1,3).imag());
        EXPECT_FLOAT_EQ(0.0,driver.getC()(1,4).imag());
        EXPECT_FLOAT_EQ(0.0,driver.getC()(1,5).imag());
        EXPECT_FLOAT_EQ(4.228396699597442e-13,driver.getC()(1,6).imag());
        EXPECT_FLOAT_EQ(-4.228396699597442e-13,driver.getC()(1,7).imag());

}

TEST_F(DriverTest,Matrix_D){

	const size_t Ns = 20;
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

	const size_t Ns = 20;
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

   EXPECT_FLOAT_EQ(-3.794285178830511e-19,driver.getE()(0,0).real());
   EXPECT_FLOAT_EQ(-7.588570357661021e-19,driver.getE()(0,1).real());
   EXPECT_FLOAT_EQ(-7.588570357661021e-19,driver.getE()(1,0).real());
   EXPECT_FLOAT_EQ(-1.641368479809876e-18,driver.getE()(1,1).real());

}



TEST_F(DriverTest,Matrix_R){

	const size_t Ns = 20;
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

	   EXPECT_FLOAT_EQ(0.0203,driver.getR()[0](0,0).real());
	   EXPECT_FLOAT_EQ(0.0407,driver.getR()[0](0,1).real());
	   EXPECT_FLOAT_EQ(0.0407,driver.getR()[0](1,0).real());
	   EXPECT_FLOAT_EQ(0.2771,driver.getR()[0](1,1).real());
	   EXPECT_FLOAT_EQ(0.0007,driver.getR()[1](0,0).real());
	   EXPECT_FLOAT_EQ(0.0015,driver.getR()[1](0,1).real());
	   EXPECT_FLOAT_EQ(0.0015,driver.getR()[1](1,0).real());
	   EXPECT_FLOAT_EQ(-0.1008,driver.getR()[1](1,1).real());
       EXPECT_FLOAT_EQ(-0.0414,driver.getR()[2](0,0).real());
       EXPECT_FLOAT_EQ(-0.0828,driver.getR()[2](0,1).real());
       EXPECT_FLOAT_EQ(-0.0828,driver.getR()[2](1,0).real());
       EXPECT_FLOAT_EQ(-0.3100,driver.getR()[2](1,1).real());
       EXPECT_FLOAT_EQ(-0.0414,driver.getR()[3](0,0).real());
       EXPECT_FLOAT_EQ(-0.0828,driver.getR()[3](0,1).real());
       EXPECT_FLOAT_EQ(-0.0828,driver.getR()[3](1,0).real());
       EXPECT_FLOAT_EQ(-0.3100,driver.getR()[3](1,1).real());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[0](0,0).imag());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[0](0,1).imag());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[0](1,0).imag());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[0](1,1).imag());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[1](0,0).imag());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[1](0,1).imag());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[1](1,0).imag());
       EXPECT_FLOAT_EQ(0.0,driver.getR()[1](1,1).imag());
       EXPECT_FLOAT_EQ(0.1113,driver.getR()[2](0,0).imag());
       EXPECT_FLOAT_EQ(0.2227,driver.getR()[2](0,1).real());
       EXPECT_FLOAT_EQ(0.2227,driver.getR()[2](1,0).real());
       EXPECT_FLOAT_EQ(0.4228,driver.getR()[2](1,1).real());
       EXPECT_FLOAT_EQ(-0.1113,driver.getR()[3](0,0).real());
       EXPECT_FLOAT_EQ(-0.2227,driver.getR()[3](0,1).real());
       EXPECT_FLOAT_EQ(-0.2227,driver.getR()[3](1,0).real());
       EXPECT_FLOAT_EQ(-0.4228,driver.getR()[3](1,1).real());

}

TEST_F(DriverTest,Matrix_poles){

	const size_t Ns = 20;
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

	  EXPECT_FLOAT_EQ(-2.912231273155982,driver.ss2pr()[0].real());
      EXPECT_FLOAT_EQ(-18.548225577281970 ,driver.ss2pr()[0].real());
      EXPECT_FLOAT_EQ(-2.340975638498256e+02,driver.ss2pr()[0].real());
      EXPECT_FLOAT_EQ(-2.340975638498256e+02,driver.ss2pr()[0].real());
      EXPECT_FLOAT_EQ(0.0,driver.ss2pr()[0].imag());
      EXPECT_FLOAT_EQ(0.0,driver.ss2pr()[0].imag());
      EXPECT_FLOAT_EQ(8.340621081056502e+02,driver.ss2pr()[0].imag());
      EXPECT_FLOAT_EQ(-8.340621081056502e+02,driver.ss2pr()[0].imag());
}



