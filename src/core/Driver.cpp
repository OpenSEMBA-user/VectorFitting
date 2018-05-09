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

Fitting Driver::getFitting() const {
    return fitting_;
}

std::vector<Fitting::Sample> Driver::squeeze(
        const std::vector<Driver::Sample>& samples){
    std::vector<MatrixXcd> response;
    for (size_t i = 0; i < samples.size(); ++i) {
        response.push_back(samples[i].second);
    }
    std::vector<VectorXcd> squeezedResponse = squeeze(response);

    std::vector<Fitting::Sample> res;
    for (size_t i = 0; i < samples.size(); ++i) {
        res.push_back(Fitting::Sample(samples[i].first, squeezedResponse[i]));
    }
    return res;
}

std::vector<Fitting::Sample> Driver::calcFsum(
        const std::vector<Fitting::Sample>& f,
        const Options& options) {
    switch (options.getWeighting()) {
    case Options::Weighting::one:
    {
        std::vector<Fitting::Sample> fSum;
        for (size_t i; i < f.size(); ++i){
            VectorXcd sum(1,1);
            sum << f[i].second.sum();
            fSum.push_back(std::make_pair(f[i].first, sum));
        }
        return fSum;
    }
    default:
        throw std::runtime_error("Weighting parameter not implemented");
    }


}

Driver::Driver(
        std::vector<Sample> samples,
        const Options& opts,
        const std::vector<Complex>& inputPoles,
        const std::vector<MatrixXd>& weights) {

    std::sort(samples.begin(), samples.end(), [](Sample a, Sample b) {
        return lower(a.first.imag(), b.first.imag());
    });

    std::vector<Complex> poles = inputPoles;
    if (poles.empty()) {
        std::pair<Real,Real> range(samples.front().first.imag(),
                                   samples.back().first.imag());
        poles = Fitting::buildPoles(range, opts);
    } else {
        poles = inputPoles;
    }

    std::vector<Fitting::Sample> squeezedSum = calcFsum(squeeze(samples), opts);
    Fitting fitting1(squeezedSum, opts, poles, squeeze(weights));
    for (size_t i = 0; i < opts.getIterations().first; ++i) {
        //lines 382-392
        fitting1.fit();
    }

    Fitting fitting2(squeeze(samples), opts, poles, squeeze(weights));
    for (size_t i = 0; i < opts.getIterations().second; ++i) {
        //lines 394-410
        if (i == opts.getIterations().second - 1) {
            fitting2.options().setSkipResidueIdentification(false);
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

	size_t Nc;
	{
        size_t tell = 0;
        for (size_t k = 0; k < 10000; ++k){
            tell += (k+1);
            if (tell == (size_t) D.size()){
                Nc = (k+1);
                break;
            }
        }
	}

	size_t tell = 0;
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
