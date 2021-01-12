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

#ifndef TENDONDRIVENROBOT_H
#define TENDONDRIVENROBOT_H

#include <Eigen/Core>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <array>


#include "cosseratrodmodel.h"
#include "constantcurvaturemodel.h"
#include "piecewiseconstantcurvaturemodel.h"
#include "pseudorigidbodymodel.h"
#include "subsegmentcosseratrodmodel.h"


// This class implements a simple wrapper for the different TDCR modeling approaches.
// Currently, the implementation is tailored to a two segment TDCR with three tendons per segment.
// Each modeling approach is implemented using its own class and can also be used
// individually.
class TendonDrivenRobot
{
public:
	// Constructor to set up a TDCR with default parameters
    TendonDrivenRobot();
    
    // Simple destructor
    ~TendonDrivenRobot();
	
	// Enumerated type to choose between different model implementations
    enum Model { CosseratRod, ConstantCurvature, PiecewiseConstantCurvature, PseudoRigidBody, SubsegmentCosseratRod };
	
	
	// This function runs the forward kinematics for the TDCR. It returns true, if the forward kinematics were solved successfully.
	//
	// Inputs:
	// q 					6x1 vector containing the actuation values for each tendon (either force in Newton or displacements in meter).
	// 	   					First three entries belong to the first segment, while the last three entries belong to the second segment.
	// f_ext				3x1 vector containing the external force acting at the tip/end-effector of the robot expressed in the base frame.	 
	// l_ext				3x1 vector containing the external moment acting at the tip/end-effector of the robot expressed in the base frame.	
	// model				The chosen modeling approach to use for solving the forward kinematics.
	//		
	// Outputs:		
	// ee_frame				4x4 matrix containing the end-effector pose as a 4x4 frame (R p; 0 0 0 1).	
    bool forwardKinematics(Eigen::Matrix4d &ee_frame, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext, Model model);
    
    
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
	// two_tendons			Specifies, if only two tendons for each segment are employed and actuated.
	//						Only affects the implementation of the Constant Curvature modeling approach (details can be found there).
    void setRobotParameters(std::array<double,2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int,2> number_disks, std::array<double,2> pradius_disks, double ro, bool two_tendons);
    
    // This function enables a continuation mode for the modeling approaches that are utilizing an optimization scheme that is based on an initial guess.
	// That means, that every new run of the forward kinematics will use the final values for these variables obtained from the last forward kinematics solution as the new initial guess.
	// This makes sense in cases, where configurations of the robot only change marginally (thus the initial guesses would be similar), 
	// and can increase computation time a lot, since the algorithm will converge in just a couple of iterations.
	//
	// Inputs:
	// keep					Boolean value to indicate, if the continuation mode is used or not.
	//						If continuation mode is disabled, the default initial state of the robot (straight, no bending) is assumed for the initial guess. 
    void keepInits(bool keep);
    
    // This function returns the current end-effector/tip frame (4x4 transformation matrix) of the robot (based on the last forward kinematics computation).
    Eigen::Matrix4d getEEFrame();
    
    // This function returns the current configuration/actuation values (6x1 vector) for each tendon of the TDCR (based on the last q that was handed to the forward kinematics computation).
    Eigen::MatrixXd getCurrentConfig();
    
    // This function returns the current disk frames (4x4 transformation matrix for each discrete disk) of the robot (based on the last forward kinematics computation) as a stacked matrix.
    // The returned matrix of size 4x(4*n), where n is the total number of disks and every nth four columns belong to the nth disk frame.
    Eigen::MatrixXd getDiskFrames();
    
    // This function returns the current tendon displacements (6x1 vector) for each tendon of the robot (based on the last forward kinematics computation).
    // This is done by evaluating the current shape of the robot and computing the length of each tendon based on the current disk frames.
    // This function can be used to obtain tendon displacements from the TDCR if it has currently been operated using a modeling approach that is based on tendon force acutation.
    Eigen::MatrixXd getTendonDisplacements();

private:

	// Pointers to the different modeling approaches
    CosseratRodModel* mp_cr_model;
    ConstantCurvatureModel* mp_cc_model;
    PiecewiseConstantCurvatureModel* mp_pcc_model;
    PseudoRigidBodyModel* mp_prb_model;
    SubsegmentCosseratRodModel* mp_sscr_model;
	
	// Member variables
    Eigen::MatrixXd m_actuation_values;
    Eigen::Vector3d m_ext_forces;
    Eigen::Vector3d m_ext_moments;
    Eigen::Matrix4d m_ee_frame;
    Eigen::MatrixXd m_disk_frames;
    bool m_keep_inits;

    //Robot Parameter
    std::array<double,2> m_length;
    double m_youngs_modulus;
    std::vector<Eigen::Vector3d> m_routing;
    std::array<int,2> m_number_disks;
    std::array<double,2> m_pradius_disks;
    double m_ro;
    bool m_two_tendons;

};

#endif // TENDONDRIVENROBOT_H
