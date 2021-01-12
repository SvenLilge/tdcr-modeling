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

#ifndef CONSTANTCURVATUREMODEL_H
#define CONSTANTCURVATUREMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <array>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <numeric>


// This class implements a simple Constant Curvature model to solve the forward kinematics of a TDCR.
// This code is explicitly written for a two segment TDCR with either three or two tendons per segment (two segments if two_tendons is true).
// The tendons are equally distributed around the backbone with the distance pradius_disks.
// The first tendon is located at the y-axis of the local disk frame, while the other tendons are offset by 120 degrees (three tendons) or 120 degrees (two tendons), respectively.
class ConstantCurvatureModel 
{
public:

	// Constructor to set up the Constant Curvature model with default parameters
	ConstantCurvatureModel();
	
	// Simple destructor
	virtual ~ConstantCurvatureModel();
	
	// This function runs the forward kinematics for the TDCR. It returns true, if the forward kinematics were solved successfully.
	//
	// Inputs:
	// q 					6x1 vector containing the actuation values for each tendon (displacements in meter).
	// 	   					First three entries belong to the first segment, while the last three entries belong to the second segment.
	//		
	// Outputs:		
	// diskFrames			4x(4*n) matrix containing the current disk frames (4x4 transformation matrix for each discrete disk) of the robot as a stacked matrix.
	bool forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q);
	
	// This function allows to set and update the TDCR parameters.
	//
	// Inputs:
	// length				std::array that holds the length of each of the two segments of the TDCR.
	// number_disks			std::array that holds the number of disks for each of the two segments of the TDCR.
	// pradius_disks		std::array that holds the pitch radius of the disks (distance between tendon routing and backbone) for each of the two segments of the TDCR.
	// two_tendons			Specifies, if only two tendons for each segment are employed and actuated.
	void setRobotParameters(std::array<double,2> length, std::array<int,2> number_disks, std::array<double,2> pradius_disks, bool two_tendons);

private:
	
	// Member variables
	std::array<double,2> m_length;
	std::vector<Eigen::Vector3d> m_routing;
	std::array<int,2> m_number_disks;
	std::array<double,2> m_pradius_disks;
	bool m_two_tendons;
};

#endif // CONSTANTCURVATUREMODEL_H
