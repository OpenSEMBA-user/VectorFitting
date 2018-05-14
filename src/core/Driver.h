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


#include "Fitting.h"
#include "SpaceGenerator.h"

#include <cmath>

namespace VectorFitting {

using namespace Eigen;

class Driver {
public:
    typedef std::pair<Complex, MatrixXcd> Sample;

	/**
	 * A fitter with starting poles computed automatically will be called
	 * from VFdriver
	 * @param samples   Data to be fitted.
	 * @param order     Order of approximation.
	 * @param options   Options.
     */
	Driver(std::vector<Sample> samples,
           const Options& options,
           const std::vector<Complex>& poles = {},
           const std::vector<MatrixXd>& weights = {});
	const MatrixXcd& getA() const;
	const MatrixXi& getB() const;
	const MatrixXcd& getC() const;
	const MatrixXcd& getD() const;
	const MatrixXcd& getE() const;
	const std::vector<Fitting::Sample>& getSamples() const;

	std::pair<std::vector<Complex>, std::vector<MatrixXcd>> ss2pr() const;
private:

	MatrixXcd A_;
	MatrixXi B_;
	MatrixXcd C_;
	MatrixXcd D_;
	MatrixXcd E_;
	std::vector<Fitting::Sample> samples_;

	template <class T>
	static T blkdiag(const T& a, const T& b) {
	    T res = T::Zero(a.rows() + b.rows(), a.cols() + b.cols());
	    for (typename T::Index i = 0; i < a.rows(); ++i) {
	        for (typename T::Index j = 0; j < a.cols(); ++j) {
	            res(i,j) = a(i,j);
	        }
	    }
	    for (auto i = 0; i < b.rows(); ++i) {
	        for (auto j = 0; j < b.cols(); ++j) {
	            res(i+a.rows(),j+a.cols()) = b(i,j);
	        }
	    }
	    return res;
	}



	static std::vector<Fitting::Sample>
	                squeeze(const std::vector<Sample>& samples);


	template <class T>
	static std::vector<Eigen::Matrix<T,-1,1>> squeeze(
	        const std::vector<Eigen::Matrix<T,-1,-1>>& rhs) {
	    for (size_t i = 0; i < rhs.size(); ++i){
	        for (auto j = 0; j < rhs[i].cols(); ++j){
	            for (auto k = j; k < rhs[i].rows(); ++k) {
	                if (!equal(((Complex) rhs[i](j,k)).real(),
	                        ((Complex) rhs[i](k,j)).real()) ||
	                        !equal(((Complex) rhs[i](j,k)).imag(),
	                                ((Complex) rhs[i](k,j)).imag())) {
	                    throw std::runtime_error(
	                            "Matrices must be symmetric to be squeezed");
	                }
	            }
	        }
	    }

	    std::vector<Eigen::Matrix<T,-1,1>> res;
	    for (size_t i = 0; i < rhs.size(); ++i) {
	        std::vector<T> aux;
	        for (auto j = 0; j < rhs[i].cols(); ++j){
	            for (auto k = j; k < rhs[i].rows(); ++k){
	                aux.push_back(rhs[i](k,j));
	            }
	        }
	        Eigen::Matrix<T,-1,1> eigAux = Fitting::toEigenVector(aux);;
	        res.push_back(eigAux);
	    }
	    return res;
	}


	static std::vector<Fitting::Sample> calcFsum(
	        const std::vector<Fitting::Sample>& f,
	        const Options& options);
	void tri2full(Fitting fitting);

};

} /* namespaceVectorFitting */
