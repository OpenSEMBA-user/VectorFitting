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

//#include </usr/include/eigen3/Eigen/Dense>

namespace VectorFitting {

using namespace Eigen;



class Driver {


public:

	/*
	 * A fitter with starting poles computed automatically will be called
	 * from VFdriver
	 * @param samples   Data to be fitted.
	 * @param order     Order of approximation.
	 * @param options   Options.
     */

	Driver(const std::vector<Sample>& samples,
             const size_t order,
             Options options,
			 std::vector<std::vector<Real>>& weights,
		     const std::pair <size_t, size_t> iterations);

	/**
	 * A fitter with starting poles provided by the user will be called
	 * from VFdriver. order_ and poles.size() shall be the same
     * @param samples   Data to be fitted.
     * @param poles     Starting poles.
     * @param options   Options.
     */

	Driver(const std::vector<Sample>& samples,
             const std::vector<Complex>& poles,
             Options options,
			 std::vector<std::vector<Real>>& weights,
			 const std::pair <size_t, size_t> iterations);


	~Driver(){

	}


protected:

	MatrixXd initWeightsSUm(std::vector<std::vector<Real>>& weights,
						  const std::vector<Sample>& samples);
	std::vector<Sample> squeeze(const std::vector<Sample>& samples);//lines 305-312
	Fitting tri2full(Fitting fitting);
	Sample calcFsum(std::vector<Sample> f);//lines 329-347
	void ss2pr(Fitting fitting);


	std::pair<size_t, size_t> iterations_;



};

} /* namespaceVectorFitting */
