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

	std::vector<std::vector<Complex>> f = VFdriver::squeeze(samples);
	std::vector<Complex> fSum = VFdriver::calcFsum(f);
	MatrixXd weightsSum = VFdriver::initWeights(weights);
	Fitting::Fitting fitting1(fSum, poles, options, weightsSum);

	for (size_t i = 0; i < iterations_.first; ++i){ //lines 382-392
		fitting1.fit;
	}

	Fitting::Fitting fitting2(f,poles,options,weights);

	for (size_t i = 0; i < iterations_.second; ++i){ //lines 394-410
		if (i == iterations_.second - 1){
			bool skipResidueIdentification = false;
			options.setSkipResidueIdentification(skipResidueIdentification);
		}

		fitting2.fit;
	}

	std::vector<std::vector<Complex>> fFull = tri2full(f);
	//set poles and set residues (lines 497-499)





	std::vector<Sample> fit;

	if (!fitting2.getFittedSamples().empty) {
		fit = fitting2.getFittedSamples();

	} else {

		fit =fitting1.getFittedSamples();
	}


	std::vector<std::vector<Complex>> diff;
	for (size_t i = 0; i < fit.size; ++i){
		std::vector<Complex> aux;
		for (size_t j = 0; j < fit[0].second.size; ++j){
			aux.push_back(samples[i].second[j]-fit[i].second[j]);
		}

		diff.push_back(aux);

	}



};

std::vector<std::vector<Complex>> Driver::squeeze(
				const std::vector<Sample>& samples){



	for (size_t i = 0; i < Fitting::getSamplesSize; ++i){
		for (size_t j = 0; j < samples[0].second.size(); ++j){
			assert(samples[1].second[j] == samples[2].second[j]);
		}
	}

	std::vector<std::vector<Complex>> f;
	for (size_t i = 0; i < Fitting::getSamplesSize; ++i){
		std::vector<Complex> aux;
		for (size_t j = 0; j < samples[0].second.size(); ++j){
			if(j!= 1){
				aux.push_back(samples[i].second[j]);
			}
		}
		f.push_back(aux);
	}


	return f;

}


std::vector<Complex> calcFsum(std::vector<std::vector<Complex>> f){

	std::vector<Complex> fSum;
	for (size_t i; i < f.size(); ++i){
		Complex sum = 0;
		for (size_t j; j < f[0].size(); ++j){
			if (j != 1){
				sum += f[i][j];
			}
		}

		fSum.push_back(sum);
	}

	return fSum;


}

MatrixXd Driver::initWeights(std::vector<std::vector<Real>>& weights){

	if (weights.size() != 0 && weights.size() != Fitting::samples_.size()) {
		throw std::runtime_error("Weights and samples must have same size.");
	}
	if (weights.size() == 0) {
		Fitting::weights_ = MatrixXd::Ones(Fitting::getSamplesSize,
				   	   	   	   	   	   	   Fitting::getResponseSize);

	} else {
		Fitting::weights_ = MatrixXd::Zero(Fitting::getSamplesSize,
										   Fitting::getResponseSize);
	    for (size_t i = 0; i < Fitting::getSamplesSize; ++i) {
	    	if (weights[i].size() != Fitting::getResponseSize) {
	    		throw std::runtime_error(
	    		 "All weights must have the same size as the samples");
	        }
            for (size_t j = 0; j < Fitting::getResponseSize; ++j) {
                Fitting::weights_(i,j) = weights[i][j];
            }
        }
    }
	MatrixXd weightsSum = MatrixXd::Ones(Fitting::getSamplesSize,
		   	   	   	   	   	   	   	   	 Fitting::getResponseSize);

	return weightsSum;
}


std::vector<std::vector<Complex>> tri2full(
			std::vector<std::vector<Complex>> f){

	std::vector<std::vector<Complex>> fFull;
	for (size_t i = 0; i < f.size(); ++i){
		std::vector<Complex> aux;
		for (size_t j = 0; j < f[0].size(); ++j){
			aux.push_back(f[i][j]);
			if (j == 1){
				aux.push_back(f[i][j]);
			}
		}
		fFull.push_back(aux);
	}
	return fFull;
}

void ss2pr() { //assumption: A,B,C are complex

	size_t Nc = Fitting::getC.cols();
	size_t N = Fitting::getA.rows() / Nc;
	std::vector<MatrixXd> R;

	for (size_t i = 0; i < N; ++i){
		MatrixXd Rdum = MatrixXd::Zero(Nc);
		for (size_t j = 0; j < Nc; ++j){
			size_t ind = j*N + i;
			for (size_t k = 0; k < Nc; ++k){
				Rdum[j][k] += Fitting::getC[k][ind] * Fitting::getB[ind][k];
			}
		}

		R.push_back(Rdum);
	}
	//a??
	Fitting::setR(R);


}


}/* namespace VectorFitting */
