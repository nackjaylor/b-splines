/*
Part of B-Splines using Eigen - a class for fitting and sampling B-splines from data.
Copyright (C) 2021  Jack Naylor

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. 
*/

#pragma once

#include "types.h"
#include "Eigen/Core"
#include <vector>
using namespace Eigen;

class BSpline {

    public:

        PointArray control_p;
        VectorXd knot_vector;
        int num_cp = control_p.cols();

        BSpline();
        BSpline(const PointArray& control_p_, const VectorXd& knot_vector_);

        void fit(const PointArray& data, const int& degree = 3, const int& num_cp = 5);

        PointArray sample(const int& num_samples,const PointArray& control_p, const VectorXd& knot_vector, const int& degree = 3) const;
        PointArray sample(const int& num_samples, const int& degree = 3) const;


        VectorXd centripetal_method(const PointArray& data);

        VectorXd resample_spacing(const VectorXd& spacing, const int& number_points);

        VectorXd solve_knot_vector(const VectorXd& spacing, const int& degree, const int& num_knots);

        VectorXd coefficients(const double& u_val, const int& num_points, const VectorXd& knot_vec, const int& degree = 3) const;



};




