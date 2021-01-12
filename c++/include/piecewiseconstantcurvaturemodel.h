/*
This code implements different approaches to model the kinematics/statics
of a two segment tendon driven continuum robot and is part of the following
publication:

How to model tendon-driven continuum robots and benchmark modelling performance
Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
frontiers in Robotics and AI 2021
DOI: 10.3389/frobt.2020.630245

Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga
*/

#ifndef PIECEWISECONSTANTCURVATUREMODEL_H
#define PIECEWISECONSTANTCURVATUREMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <array>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <numeric>

// This class implements a Piecewise Constant Curvature model to solve the forward kinematics of a TDCR.
// This code is written for a two segment TDCR with to three tendons per segment subject to an external load at the tip.
class PiecewiseConstantCurvatureModel 
{
public:

	// Constructor to set up the Piecewise Constant Curvature model with default parameters
    PiecewiseConstantCurvatureModel();
    
	// Simple destructor
	virtual ~PiecewiseConstantCurvatureModel();

    // This function runs the forward kinematics for the TDCR. It returns true, if the forward kinematics were solved successfully.
	//
	// Inputs:
	// q 					6x1 vector containing the actuation values for each tendon (either force in Newton or displacements in meter).
	// 	   					First three entries belong to the first segment, while the last three entries belong to the second segment.
	// f_ext				3x1 vector containing the external force acting at the tip/end-effector of the robot expressed in the base frame.	 
	// l_ext				3x1 vector containing the external moment acting at the tip/end-effector of the robot expressed in the base frame.	
	//		
	// Outputs:		
	// diskFrames			4x(4*n) matrix containing the current disk frames (4x4 transformation matrix for each discrete disk) of the robot as a stacked matrix.	
    bool forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext);
    
    // This function implements the boundary conditions that need to be solved for the Cosserat rod model (force and moment equilibirium at the tip of the robot).
	//
	// Inputs:
	// output				3nx1 vector containing the moment equilibirium at the end of each subsegment.
	//		
	// Outputs:		
	// input				3nx1 vector containing the current values for each subsegment (beta, gamma, epsilon).
    void get_res(Eigen::MatrixXd &output, Eigen::MatrixXd input);
    
    // This function enables a continuation mode.
	// That means, that every new run of the forward kinematics will use the final values for these variables obtained from the last forward kinematics solution as the new initial guess.
	// This makes sense in cases, where configurations of the robot only change marginally (thus the initial guesses would be similar), 
	// and can increase computation time a lot, since the algorithm will converge in just a couple of iterations.
	//
	// Inputs:
	// keep					Boolean value to indicate, if the continuation mode is used or not.
	//						If continuation mode is disabled, the default initial values are used for the initial guess
	//						(by default this is the straight, not bent state of the robot, but they can also be changed - see below). 
    void setKeepInits(bool keep);
    
    // This function allows to set and update the TDCR parameters.
	//
	// Inputs:
	// length				std::array that holds the length of each of the two segments of the TDCR.
	// youngs_modulus		Youngs modulus of the backbone of the TDCR.	 
	// routing				std::vector that holds the routing position of each tendon of the TDCR expressed as a 3x1 position vector in the local disk frame.
	// 	   					First three entries belong to the first segment, while the last three entries belong to the second segment.
	// number_disks			std::array that holds the number of disks for each of the two segments of the TDCR.
	// pradius_disks		std::array that holds the pitch radius of the disks (distance between tendon routing and backbone) for each of the two segments of the TDCR.
	// ro					Radius of the backbone of the TDCR.	
    void setRobotParameters(std::array<double,2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int,2> number_disks, std::array<double,2> pradius_disks, double ro);
    
    // This function returns the number of total disks throughout the whole robot structure.
    double getNumberOfTotalDisks();
    
    // This function allows to set the default initial values that are used as an initial guess for the implemented shooting method.
	//
	// Inputs:
	// inits				3nx1 vector containing the default values that are used as an initial guess for the values (beta, gamma, epsilon) for each subsegment, if continuation mode is disabled.
    void setDefaultInitValues(Eigen::MatrixXd inits);
    
	// This function returns the last values final values of the shooting method for the last solution of the forward kinematics.
	// It returns 3nx1 vector containing the last values (beta, gamma, epsilon) for each subsegment.
    Eigen::MatrixXd getFinalInitValues();
    
    // This function returns a double value containing the least squares resudial of the boundary conditions for the last solution of the forward kinematics.
    double getFinalResudial();

private:
	
	// Member function returning the transformation matrix between two disks indexed by idx_from and idx_to.
    Eigen::Matrix4d getDiskTransform(Eigen::MatrixXd var, Eigen::VectorXd length_ss, int idx_from, int idx_to);

	// Member variables
    std::array<double,2> m_length;
    double m_youngs_modulus;
    std::vector<Eigen::Vector3d> m_routing;
    std::array<int,2> m_number_disks;
    std::array<double,2> m_pradius_disks;
    double m_ro;

    bool m_keep_inits;

    Eigen::Vector3d m_f_ext;
    Eigen::Vector3d m_l_ext;

    Eigen::Matrix<double,6,1> m_tau;

    Eigen::MatrixXd m_last_inits;
    Eigen::MatrixXd m_default_inits;

    double m_residual_error;

};

#endif // PIECEWISECONSTANTCURVATUREMODEL_H
