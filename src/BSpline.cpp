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

#include "types.h"

#include "BSpline.h"
#include "closest_floor_index.h"
#include <cmath>


BSpline::BSpline () {

}

BSpline::BSpline(const PointArray& control_p_, const VectorXd& knot_vector_)
    : control_p(control_p_), knot_vector(knot_vector_) {}

void BSpline::fit(const PointArray& data, const int& degree, const int& num_cp) {
    int n = data.cols()-1;
    int h = num_cp -1;
    int m = h+1+degree;


    VectorXd spacing = centripetal_method(data);


    VectorXd new_spacing = resample_spacing(spacing, h+1);

    VectorXd knots = solve_knot_vector(new_spacing, degree, m+1);

    knot_vector = knots;


    MatrixXd N_matrix;
    N_matrix.resize(n-1,h-1);

    Matrix<double, 3, Dynamic> Q_points;
    Q_points.resize(NoChange,data.cols());

    for (int i = 0; i<n+1; i++) {
        VectorXd coeffs = coefficients(spacing(i), h+1, knots, degree);

        
        if (i>0 && i < n) {

            N_matrix.row(i-1) = coeffs(seq(1,h-1));
        }

        Q_points.col(i) << data.col(i) - coeffs(0)*data.col(0) - coeffs(h)*data.col(n);

    }


    Matrix<double,Dynamic,3> Q_matrix;
    Q_matrix.resize(h-1,NoChange);
    Q_matrix.setZero();

    for (int j = 0; j<h-1; j++) {

        for (int k = 1; k<n; k++) {

            Q_matrix.row(j) += N_matrix(k-1,j)*Q_points.col(k).transpose();
        }
    }

    MatrixXd M_matrix = N_matrix.transpose()*N_matrix;


    PointArray p_points = M_matrix.bdcSvd(ComputeThinU | ComputeThinV).solve(Q_matrix).transpose();

    PointArray output_points;
    output_points.resize(NoChange,num_cp);
    output_points.col(0) = data.col(0);
    output_points.col(h) = data.col(n);
    output_points.block(0,1,3,h-1) = p_points;

    control_p = output_points;


}

VectorXd BSpline::resample_spacing(const VectorXd& spacing, const int& number_points) {

    int prev_s_ind = 1;
    VectorXd curr_samples = VectorXd::LinSpaced(spacing.size(),0,1);
    VectorXd wanted_samples = VectorXd::LinSpaced(number_points,0,1);
    VectorXd resampled_spacing;
    resampled_spacing.resize(number_points);
    resampled_spacing(0) = 0;
    resampled_spacing(number_points-1) = 1;

    for (int i = 1; i< number_points-1; i++) {
        
        for (int j = prev_s_ind; j<spacing.size()-1; j++) {
            
            if (curr_samples(j) == wanted_samples(i)) {

                resampled_spacing(i) = spacing(j);

                prev_s_ind = j;

                break;

            } else if (curr_samples(j) < wanted_samples(i)) {
                
                resampled_spacing(i) = spacing(j) + (spacing(j+1)-spacing(j))*(wanted_samples(i)-curr_samples(j))/(curr_samples(j+1)-curr_samples(j));
                
                prev_s_ind = j;
                break;
            }
        }
    }

    return resampled_spacing;
}

VectorXd BSpline::centripetal_method(const PointArray& data) {

    VectorXd spacing_vec;
    spacing_vec.resize(data.cols());

    double pow_a = 1/2;

    VectorXd interm_vec;
    interm_vec.resize(data.cols()-1);

    for (int i = 1; i<data.cols()-1; i++) {
        interm_vec(i) = pow((data(all,i)-data(all,i-1)).norm(),pow_a);
    }

    double length_polygon = interm_vec.sum();

    spacing_vec(0) = 0;
    for (int i = 1; i<data.cols(); i++) {
        spacing_vec(i) = interm_vec(seq(0,i-1)).sum()/length_polygon;
    }

    return spacing_vec;

}

VectorXd BSpline::solve_knot_vector(const VectorXd& spacing, const int& degree, const int& num_knots) {


    VectorXd knots;

    knots.resize(num_knots);
    
    for (int i = 0; i<degree+1; i++) {
        knots(i) = 0;
    }



    for (int i = 1; i<num_knots-2*degree-1; i++) {

        knots(i+degree) = spacing(seq(i,i+degree-1)).sum()/degree;


    }

    for (int i = num_knots-degree-1; i < num_knots; i++) {
        knots(i) = 1;
    }

    return knots;
}

VectorXd BSpline::coefficients(const double& u_val, const int& num_points, const VectorXd& knot_vec, const int& degree) const {
    
    VectorXd solved_coeffs;
    solved_coeffs.resize(num_points);
    solved_coeffs.setZero();

    if (u_val == knot_vec(0)) {
 
        solved_coeffs(0) = 1;
        return solved_coeffs;
    } else if (u_val == knot_vec(knot_vec.size()-1)) {
        solved_coeffs(num_points-1) = 1;
        return solved_coeffs;
    }


    

    int knot_ind = closest_floor_index(knot_vec,u_val);

    solved_coeffs(knot_ind) = 1;
    for (int d = 1; d < degree+1; d++) {
        solved_coeffs(knot_ind-d) = (knot_vec(knot_ind+1)-u_val)/(knot_vec(knot_ind+1)-knot_vec(knot_ind-d+1))*solved_coeffs(knot_ind-d+1);

        for (int i = knot_ind-d+1; i < knot_ind; i++) {
            solved_coeffs(i) = (u_val-knot_vec(i))/(knot_vec(i+d)-knot_vec(i))*solved_coeffs(i)+ (knot_vec(i+d+1)-u_val)/(knot_vec(i+d+1)-knot_vec(i+1))*solved_coeffs(i+1);
        }

        solved_coeffs(knot_ind) = (u_val-knot_vec(knot_ind))/(knot_vec(knot_ind+d)-knot_vec(knot_ind))*solved_coeffs(knot_ind);
    }

    return solved_coeffs;
}

PointArray BSpline::sample(const int& num_samples, const PointArray& control_p_, const VectorXd& knot_vector_, const int& degree) const {
    Matrix<double, 3, Dynamic> spline_samples;
    spline_samples.resize(NoChange,num_samples);
    spline_samples.setZero();
    VectorXd sample_vec = VectorXd::LinSpaced(num_samples,0,1);

    for (int i = 0; i<num_samples; i++) {

        for (int j = 0; j < control_p_.cols(); j++) {

            spline_samples.col(i) += coefficients(sample_vec(i),control_p_.cols(),knot_vector_,degree)(j)*control_p_.col(j);
        }
    }

    return spline_samples;
}






PointArray BSpline::sample(const int& num_samples, const int& degree) const {

    Matrix<double, 3, Dynamic> spline_samples;
    spline_samples.resize(NoChange,num_samples);
    spline_samples.setZero();

    VectorXd sample_vec = VectorXd::LinSpaced(num_samples,0,1);

    for (int i = 0; i<num_samples; i++) {

        for (int j = 0; j < control_p.cols(); j++) {

            spline_samples.col(i) += coefficients(sample_vec(i),control_p.cols(),knot_vector,degree)(j)*control_p.col(j);
        }
    }

    return spline_samples;

}