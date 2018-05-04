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

Driver::Driver(const std::vector<Sample>& samples,
               const std::vector<Complex>& poles,
               Options options,
			   std::vector<std::vector<Real>>& weights,
			   const std::pair <size_t, size_t> iterations) :
			   iterations_(iterations){

	MatrixXd weightsSum = Driver::initWeightsSUm(weights,samples);
	std::vector<Sample> f = Driver::squeeze(samples);
	Sample fSum = Driver::calcFsum(f);

	Fitting fitting1(fSum, poles, options, weightsSum);


	for (size_t i = 0; i < iterations_.first; ++i){ //lines 382-392
		fitting1.fit();
	}

	Fitting fitting2(f,poles,options,weights);
	fitting2.initWeights(weights);

	for (size_t i = 0; i < iterations_.second; ++i){ //lines 394-410
		if (i == iterations_.second - 1){
			bool skipResidueIdentification = false;
			options.setSkipResidueIdentification(skipResidueIdentification);
		}

		fitting2.fit();
	}

	Fitting fittingFull = tri2full(fitting2);
	//set poles and set residues (lines 497-499)

	ss2pr(fitting2);

	std::vector<Sample> fit;

	if (fitting2.getFittedSamples().empty()) {
		fit = fitting2.getFittedSamples();

	} else {

		fit =fitting1.getFittedSamples();
	}


	std::vector<std::vector<Complex>> diff;
	for (size_t i = 0; i < fit.size(); ++i){
		std::vector<Complex> aux;
		for (size_t j = 0; j < fit[i].second.size(); ++j){
			aux.push_back(samples[i].second[j]-fit[i].second[j]);
		}

		diff.push_back(aux);

	}



};

MatrixXd Driver::initWeightsSUm(
		std::vector<std::vector<Real>>& weights,
		const std::vector<Sample>& samples){

	MatrixXd weightsSum = MatrixXd::Ones(samples.size(),
		   	   	   	   	   	   	   	   	 samples[0].second.size());
	return weightsSum;
}

std::vector<Sample> Driver::squeeze(const std::vector<Sample>& samples){



	for (size_t i = 0; i < samples.size(); ++i){
		for (size_t j = 0; j < samples[i].second.size(); ++j){
			assert(samples[1].second[j] == samples[2].second[j]);
		}
	}

	std::vector<Sample> f;
	for (size_t i = 0; i < samples.size(); ++i){
		Complex aux1 = samples[i].first;
		std::vector<Complex> aux2;
		for (size_t j = 0; j < samples[i].second.size(); ++j){
			if(j!= 1){
				aux2.push_back(samples[i].second[j]);
			}
		}
		f.push_back(std::make_pair(aux1,aux2));

	}


	return f;

}


Sample calcFsum(std::vector<Sample> f){

	Sample fSum;
	for (size_t i; i < f.size(); ++i){
		Complex aux = f[i].first;
//		size_t N = f.size();
		std::vector<Complex> sum(f.size());
		for (size_t j; j < f[i].second.size(); ++j){
			if (j != 1){
				sum[i] += f[i].second[j];
			}
		}

		fSum.first = aux;
		fSum.second = sum;
	}

	return fSum;


}


Fitting tri2full(Fitting fitting){

	MatrixXcd A = fitting.getA();
	MatrixXcd AA(1,1);
	MatrixXcd C = fitting.getC();
	MatrixXcd CC;
	VectorXcd D = fitting.getD();
	MatrixXcd DD(1,1);
	VectorXcd E = fitting.getE();
	MatrixXcd EE(1,1);
	RowVectorXi B = fitting.getB();
	MatrixXi BB(1,1);

	size_t tell = 0;
	size_t Nc;

	for (size_t k = 0; k < 10000; ++k){
		tell += k;
		if (tell == D.size()){
			Nc = k;
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
		AA.topLeftCorner(AA.rows(), AA.cols()) = AA;
		AA.bottomRightCorner(A.rows(), A.cols()) = A;
		BB.topLeftCorner(BB.rows(), BB.cols()) = BB;
		BB.bottomRightCorner(B.rows(), B.cols()) = B;
		for (size_t j = i; j < Nc; ++j){
			tell += 1;
			DD.conservativeResize(j+1,Eigen::NoChange);
			DD(i, j) = D(tell);
			EE.conservativeResize(j+1,Eigen::NoChange);
			EE(i,j) = E(tell);
			for (size_t k = 0; k < i*N; ++k){
				for (size_t m = 0; m < C.cols(); ++m){
					CC(i,(j-1)*N+ k) = C(tell,m);
				}
			}
			for (size_t l = 0; l < j*N; ++l){
				for (size_t m = 0; m < C.cols(); ++m){
					CC(i,(i-1)*N + l) = C(tell,m);
				}
			}
		}
	}

	DD += DD - DD.eigenvalues().asDiagonal();

//	EE = EE - EE.diagonal() << EE.eigenvalues();

	fitting.setA(AA);
	fitting.setB(BB);
	fitting.setC(CC);
	fitting.setD(DD);
	fitting.setE(EE);

	return fitting;


}

void ss2pr(Fitting fitting) { //assumption: A,B,C are complex

	size_t Nc = fitting.getC().cols();
	size_t N = fitting.getA().rows() / Nc;
	std::vector<MatrixXcd> R;
	MatrixXcd C = fitting.getC();
	MatrixXcd B = fitting.getB();

	for (size_t i = 0; i < N; ++i){
		MatrixXcd Raux = MatrixXcd::Zero(Nc,Nc);
		for (size_t j = 0; j < Nc; ++j){
			size_t ind = j*N + i;
			for (size_t k = 0; k < Nc; ++k){
				Raux(i,k) += C(k,ind) * B(ind,k);
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
