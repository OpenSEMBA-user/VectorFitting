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
	Driver(const std::vector<Sample>& samples,
           const Options& options,
           const std::vector<Complex>& poles = {},
           const std::vector<std::vector<Real>>& weights = {});

	Fitting getFitting() const;

private:
	Fitting fitting_;

	template <class T>
	static T blkdiag(const T& a, const T& b) {
	    T res = T::Zero(a.rows() + b.rows(), a.cols() + b.cols());
	    for (MatrixXcd::Index i = 0; i < a.rows(); ++i) {
	        for (MatrixXcd::Index j = 0; j < a.cols(); ++j) {
	            res(i,j) = a(i,j);
	        }
	    }
	    for (auto i = a.rows(); i < a.rows() + b.rows(); ++i) {
	        for (auto j = a.cols(); j < a.cols() + b.cols(); ++j) {
	            res(i,j) = b(i,j);
	        }
	    }
	    return res;
	}

	static std::vector<Fitting::Sample>
	                squeeze(const std::vector<Sample>& samples);
	static std::vector<Fitting::Sample> calcFsum(
	        const std::vector<Fitting::Sample>& f,
	        const Options& options);
	static void ss2pr(Fitting fitting);
	static void tri2full(Fitting fitting);
    void init_(
            const std::vector<Sample>& samples,
            const std::vector<Complex>& poles,
            Options options,
            const std::vector<std::vector<Real> >& weights,
            std::pair<size_t, size_t> iterations);
};

} /* namespaceVectorFitting */
