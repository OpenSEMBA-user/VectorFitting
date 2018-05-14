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


#include "Driver.h"

namespace VectorFitting {


Driver::Driver(
        std::vector<Sample> samples2,
        const Options& opts,
        const std::vector<Complex>& inputPoles,
        const std::vector<MatrixXd>& weights) {

    std::vector<Sample>& samples = samples2; // GDB bug.

    std::sort(samples.begin(), samples.end(), [](Sample a, Sample b) {
        return lower(a.first.imag(), b.first.imag());
    });

    std::vector<Complex> poles = inputPoles;
    if (poles.empty() && !samples.empty()) {
        std::pair<Real,Real> range(samples.front().first.imag(),
                                   samples.back().first.imag());
        poles = Fitting::buildPoles(range, opts);
    } else {
        poles = inputPoles;
    }

    std::vector<Fitting::Sample> squeezedSum = calcFsum(squeeze(samples), opts);
    Fitting fitting1(squeezedSum, opts, poles, squeeze(weights));
    for (size_t i = 0; i < opts.getIterations().first; ++i) {

        fitting1.fit();
    }

    Fitting fitting2(squeeze(samples), opts, poles, squeeze(weights));
    for (size_t i = 0; i < opts.getIterations().second; ++i) {

        if (i == opts.getIterations().second - 1) {
            fitting2.options().setSkipResidueIdentification(false);
        }
        fitting2.fit();
    }
    tri2full(fitting2);
    std::pair<std::vector<Complex>, std::vector<MatrixXcd>> res = ss2pr();

    if (!fitting2.getFittedSamples().empty()) {
        samples_ = fitting2.getFittedSamples();
    } else {
        samples_ = fitting1.getFittedSamples();
    }
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
        for (size_t i = 0; i < f.size(); ++i){
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



void Driver::tri2full(Fitting fitting){

	MatrixXcd A = fitting.getA();
	VectorXcd E = fitting.getE();
	MatrixXi B =  fitting.getB();
	MatrixXcd C = fitting.getC();
	VectorXcd D = fitting.getD();
	const size_t N = A.cols();

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

	MatrixXcd AA;
	MatrixXi BB;
	MatrixXcd CC(Nc,Nc*N);
	MatrixXcd DD(Nc,Nc);
	MatrixXcd EE(Nc,Nc);

	size_t tell = 0;
	for (size_t i = 0; i < Nc; ++i){
		AA = blkdiag(AA, A);
		BB = blkdiag(BB, B);
		for (size_t j = i; j < Nc; ++j){
			DD(i,j) = D(tell);
			EE(i,j) = E(tell);
			if (i != j){
				DD(j,i) = D(tell);
				E(j,i) = E(tell);
			}
			for (size_t k = 0; k < i*N; ++k){
				for (MatrixXcd::Index m = 0; m < C.cols(); ++m){
					CC(i,j*N + k) = C(tell,m);
				}
			}
			for (size_t l = 0; l < j*N; ++l){
				for (MatrixXcd::Index m = 0; m < C.cols(); ++m){
					CC(i,i*N + l) = C(tell,m);
				}
			}
			tell++;
		}
	}

	A_ = AA;
	B_ = BB;
	C_ = CC;
	D_ = DD;
	E_ = EE;

}

std::pair<std::vector<Complex>, std::vector<MatrixXcd>> Driver::ss2pr() const {

	size_t Nr = C_.rows();
	size_t N = A_.rows() / Nr;

	std::vector<MatrixXcd> R;
	for (size_t i = 0; i < N; ++i){
		MatrixXcd Raux = MatrixXcd::Zero(Nr,Nr);
		std::vector<Complex> aux(N);
		for (size_t j = 0; j < Nr; ++j){
			const size_t ind = j*N + i;
			for (size_t k = 0; k < Nr; ++k){
				aux[i] += ((Complex) C_(k,ind)) * ((double) B_(ind,k));
			}
		}
		Raux(0,0) = aux[0];
		Raux(0,1) = aux[1];
		Raux(1,0) = aux[2];
		Raux(1,1) = aux[3];
		R.push_back(Raux);
 	}

	std::vector<Complex> poles;
	for (size_t i = 0; i < N; ++i){
		for (size_t j = 0; j < N; ++j){
			if (j == i){
				poles.push_back(A_(i,j));
			}
		}
	}

	return {poles, R};
}

const MatrixXcd& Driver::getA() const {
	return A_;
}

const MatrixXi& Driver::getB() const {
	return B_;
}

const MatrixXcd& Driver::getC() const {
	return C_;
}

const MatrixXcd& Driver::getD() const {
	return D_;
}

const MatrixXcd& Driver::getE() const {
	return E_;
}


const std::vector<Fitting::Sample>& Driver::getSamples() const {
	return samples_;
}
}/* namespace VectorFitting */
