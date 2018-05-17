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
        poles = buildPoles(range, opts);
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

std::vector<Complex> Driver::buildPoles(
        const std::pair<Real, Real>& range,
        const Options& options) {
    if (options.getPolesType() == Options::PolesType::lincmplx) {
        std::vector<Real> imagParts = linspace(range, options.getN()/2);
        std::vector<Complex> poles(options.getN());
        for (size_t i = 0; i < options.getN(); i+=2) {
            Real imag = - imagParts[i/2];
            Real real = imag *  options.getNu();
            poles[i] = Complex(real, imag);
            poles[i+1] = conj(poles[i]);
        }

        if (options.getN() % 2 != 0) {
            std::complex<Real> extraPole;
            extraPole = -(range.first + range.second)/2.0;
            poles.push_back(extraPole);
        }
        return poles;
    } else {
        throw std::runtime_error(
                "log distributed initial poles hasn't been implemented yet");
    }
}


std::vector<Fitting::Sample> Driver::calcFsum(
        const std::vector<Fitting::Sample>& f,
        const Options& options) {
    switch (options.getWeighting()) {
    case Options::Weighting::one:
    {
        std::vector<Fitting::Sample> fSum;
        for (size_t i = 0; i < f.size(); ++i){
            VectorXcd sum(1);
            sum << f[i].second.sum();
            fSum.push_back(std::make_pair(f[i].first, sum));
        }
        return fSum;
    }
    default:
        throw std::runtime_error("Weighting parameter not implemented");
    }


}



void Driver::tri2full(const Fitting& fitting) {

	const size_t N = fitting.getOrder();

	size_t Nc;
	{
        size_t tell = 0;
        for (size_t k = 0; k < 10000; ++k){
            tell += (k+1);
            if (tell == (size_t) fitting.getD().size()){
                Nc = (k+1);
                break;
            }
        }
	}

    C_ = MatrixXcd::Zero(Nc,Nc*N);
    D_ = MatrixXcd::Zero(Nc,Nc);
    E_ = MatrixXcd::Zero(Nc,Nc);

    size_t tell = 0;
	for (size_t i = 0; i < Nc; ++i){
		A_ = blkdiag(A_, fitting.getA());
		MatrixXi auxB(fitting.getB());
		B_ = blkdiag(B_, auxB);
		for (size_t j = i; j < Nc; ++j){
			D_(i,j) = fitting.getD()(tell);
			E_(i,j) = fitting.getE()(tell);
			if (i != j){
				D_(j,i) = fitting.getD()(tell);
				E_(j,i) = fitting.getE()(tell);
			}
			for (size_t k = 0; k < i*N; ++k){
				for (MatrixXcd::Index m = 0; m < N; ++m){
					C_(i,j*N + k) = fitting.getC()(tell,m);
				}
			}
			for (size_t l = 0; l < j*N; ++l){
				for (MatrixXcd::Index m = 0; m < N; ++m){
					C_(i,i*N + l) = fitting.getC()(tell,m);
				}
			}
			tell++;
		}
	}

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
