// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
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


#include "Driver.h"

namespace VectorFitting {

void Driver::init_(
        const std::vector<Sample>& samples,
        const std::vector<Complex>& poles,
        Options options,
        const std::vector<std::vector<Real> >& weights,
        const std::pair<size_t, size_t> iterations) {
    Fitting fitting1(
            { calcFsum(squeeze(samples)) },
            poles,
            options,
            weights);
    for (size_t i = 0; i < iterations.first; ++i) {
        //lines 382-392
        fitting1.fit();
    }
    Fitting fitting2(Driver::squeeze(samples), poles, options, weights);
    for (size_t i = 0; i < iterations.second; ++i) {
        //lines 394-410
        if (i == iterations.second - 1) {
            bool skipResidueIdentification = false;
            options.setSkipResidueIdentification(skipResidueIdentification);
        }
        fitting2.fit();
    }
    tri2full(fitting2);
    ss2pr(fitting2);

    if (!fitting2.getFittedSamples().empty()) {
        fitting_ = fitting2;
    } else {
        fitting_ = fitting1;
    }
}

Fitting Driver::getFitting() const {
    return fitting_;
}

Driver::Driver(const std::vector<Sample>& samples,
               const std::vector<Complex>& poles,
               Options options,
			   const std::vector<std::vector<Real>>& weights,
			   std::pair <size_t, size_t> iterations) {

    init_(samples, poles, options, weights, iterations);
};

std::vector<Fitting::Sample> Driver::squeeze(
        const std::vector<Driver::Sample>& samples){
	for (size_t i = 0; i < samples.size(); ++i){
		for (auto j = 0; j < samples[i].second.cols(); ++j){
			for (auto k = j; k < samples[i].second.rows(); ++k) {
				const MatrixXcd& data = samples[i].second;
				if (!equal(((Complex) data(j,k)).real(),
				           ((Complex) data(k,j)).real()) ||
					!equal(((Complex) data(j,k)).imag(),
						   ((Complex) data(k,j)).imag())) {
					throw std::runtime_error(
						  "Matrices in samples must be symmetric");
				}
			}
		}
	}

	std::vector<Fitting::Sample> squeezedSample;
	for (size_t i = 0; i < samples.size(); ++i) {
		Complex aux1 = samples[i].first;
		std::vector<Complex> aux2;
		for (auto j = 0; j < samples[i].second.cols(); ++j){
			for (auto k = j; k < samples[i].second.rows(); ++k){
				aux2.push_back(samples[i].second(k,j));
			}
		}
		squeezedSample.push_back(std::make_pair(aux1,aux2));
	}
	return squeezedSample;
}


Fitting::Sample Driver::calcFsum(const std::vector<Fitting::Sample>& f) {

	Fitting::Sample fSum;
	std::vector<Complex> sum(f.size());
	for (size_t i; i < f.size(); ++i){
		Complex aux = f[i].first;
		for (size_t j; j < f[i].second.size(); ++j){
			sum[i] += f[i].second[j];
		}

		fSum.first = aux;
		fSum.second = sum;
	}

	return fSum;
}

Driver::Driver(
        const std::vector<Sample>& samples,
        const size_t order,
        Options options,
        const std::vector<std::vector<Real> >& weights,
        std::pair<size_t, size_t> iterations) {
    Fitting fittingPoles(squeeze(samples), order, options, weights);
    init_(samples, fittingPoles.getPoles(), options, weights, iterations);
}

void Driver::tri2full(Fitting fitting){

	MatrixXcd A = fitting.getA();
	MatrixXcd AA(1,1);
	MatrixXcd C = fitting.getC();
	MatrixXcd CC;
	VectorXcd D = fitting.getD();
	MatrixXcd DD(1,1);
	VectorXcd E = fitting.getE();
	MatrixXcd EE(1,1);
	MatrixXi B = fitting.getB();
	MatrixXi BB(1,1);

	size_t tell = 0;
	size_t Nc;

	for (size_t k = 0; k < 10000; ++k){
		tell += (k+1);
		if (tell == (size_t) D.size()){
			Nc = (k+1);
			break;
		}
	}

	tell = 0;
	size_t N = A.cols();
	AA = {};
	BB = {};
	CC = MatrixXcd::Zero(Nc,Nc*N);
	for (size_t i = 0; i < Nc; ++i){
		DD.conservativeResize(Eigen::NoChange, i+1);
		EE.conservativeResize(Eigen::NoChange, i+1);
		AA = blkdiag(AA, A);
		BB = blkdiag(BB, B);
		for (size_t j = i; j < Nc; ++j){
			DD.conservativeResize(j+1,Eigen::NoChange);
			DD(i,j) = D(tell);
			EE.conservativeResize(j+1,Eigen::NoChange);
			EE(i,j) = E(tell);
			for (size_t k = 0; k < i*N; ++k){
				for (auto m = 0; m < C.cols(); ++m){
					CC(i,(j-1)*N+ k) = C(tell,m);
				}
			}
			for (size_t l = 0; l < j*N; ++l){
				for (auto m = 0; m < C.cols(); ++m){
					CC(i,(i-1)*N + l) = C(tell,m);
				}
			}
			tell++;
		}
	}

	MatrixXcd aux1 = DD.eigenvalues().diagonal().transpose();
	DD += aux1;
	MatrixXcd aux2 = EE.eigenvalues().diagonal().transpose();
	EE += aux2;

	fitting.setA(AA);
	fitting.setB(BB);
	fitting.setC(CC);
	fitting.setD(DD);
	fitting.setE(EE);

}

void Driver::ss2pr(Fitting fitting) { //assumption: A & C are complex

	size_t Nc = fitting.getC().cols();
	size_t N = (fitting.getA().rows()) / Nc;
	std::vector<MatrixXcd> R;
	MatrixXcd C = fitting.getC();
	MatrixXi B = fitting.getB();

	for (size_t i = 0; i < N; ++i){
		MatrixXcd Raux = MatrixXcd::Zero(Nc,Nc);
		for (size_t j = 0; j < Nc; ++j){
			size_t ind = j*N + i;
			for (size_t k = 0; k < Nc; ++k){
			    Raux(i,k) += ((Complex) C(k,ind)) * ((double) B(ind,k));
			}
		}
		R.push_back(Raux);
	}
	fitting.setR(R);
	MatrixXcd A = fitting.getA();
	std::vector<Complex> poles;

	for (size_t i = 0; i < N; ++i){
		for (size_t j = 0; j < N; ++j){
			if (j == i){
				poles.push_back(A(i,j));
			}
		}
	}
}

}/* namespace VectorFitting */
