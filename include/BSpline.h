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
        /**
         * @brief 3xn array of control points
         * 
         */
        PointArray control_p;

        /**
         * @brief nx1 vector of knots
         * 
         */
        VectorXd knot_vector;

        /**
         * @brief Number of control points
         * 
         */
        int num_cp = control_p.cols();

        /**
         * @brief Construct a new BSpline object
         * 
         */
        BSpline();

        /**
         * @brief Construct a new BSpline object
         * 
         * @param control_p_ 3xn array of control points
         * @param knot_vector_ nx1 vector of knots
         */
        BSpline(const PointArray& control_p_, const VectorXd& knot_vector_);

        /**
         * @brief Fit a B-spline to input data of given degree
         * 
         * @param data 3xn array of *sorted* data (i.e. points are sequential)
         * @param degree Degree of B-Spline
         * @param num_cp Number of control points
         */
        void fit(const PointArray& data, const int& degree = 3, const int& num_cp = 5);

        /**
         * @brief Sample a B-Spline instance
         * 
         * @param num_samples Number of points on curve to sample
         * @param control_p Control point array
         * @param knot_vector nx1 knot vector
         * @param degree Degree of B-Spline curve
         * @return PointArray 3xn array of samples
         */
        static PointArray sample(const int& num_samples,const PointArray& control_p, const VectorXd& knot_vector, const int& degree = 3);

        /**
         * @brief Sample the current B-Spline
         * 
         * @param num_samples Number of points on curve to sample
         * @param degree Degree of B-Spline
         * @return PointArray 3xn array of samples
         */
        PointArray sample(const int& num_samples, const int& degree = 3) const;

        /**
         * @brief Invoke the centripetal method to calculate spacing of points
         * 
         * @param data 3xn data array of which to sample
         * @return VectorXd (n-1)x1 spacing vector
         */
        VectorXd centripetal_method(const PointArray& data);

        /**
         * @brief Resample spacing vector with new points using linear interpolation
         * 
         * @param spacing Current vector
         * @param number_points number of points to resample to
         * @return VectorXd (number_points)x1 vector
         */
        VectorXd resample_spacing(const VectorXd& spacing, const int& number_points);

        /**
         * @brief Solve knot vector given spacing
         * 
         * @param spacing Vector of data spacings
         * @param degree Degree of spline
         * @param num_knots Number of knots in vector
         * @return VectorXd (num_knots)x1 vector of knots
         */
        VectorXd solve_knot_vector(const VectorXd& spacing, const int& degree, const int& num_knots);

        /**
         * @brief Coefficients of B-Spline Curve
         * 
         * @param u_val Sample value
         * @param num_points Number of points
         * @param knot_vec Knot vector
         * @param degree Degree of B-Spline curve
         * @return VectorXd (num_points)x1 vector of coefficients
         */
        static VectorXd coefficients(const double& u_val, const int& num_points, const VectorXd& knot_vec, const int& degree = 3);



};




