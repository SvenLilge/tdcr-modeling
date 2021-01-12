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

#include "constantcurvaturemodel.h"

ConstantCurvatureModel::ConstantCurvatureModel()
{
	//Initialize member variables and set up default parameters for TDCR
    m_length[0]         = 0.1;
    m_length[1]         = 0.1;
    m_number_disks[0]   = 10;
    m_number_disks[1]   = 10;
    m_pradius_disks[0]  = 0.008;
    m_pradius_disks[1]  = 0.006;
    m_two_tendons       = false;

}

void ConstantCurvatureModel::setRobotParameters(std::array<double, 2> length, std::array<int, 2> number_disks, std::array<double, 2> pradius_disks, bool two_tendons)
{
	
	//Update parameters
    m_length[0]         = length[0];
    m_length[1]         = length[1];
    m_number_disks[0]   = number_disks[0];
    m_number_disks[1]   = number_disks[1];
    m_pradius_disks[0]  = pradius_disks[0];
    m_pradius_disks[1]  = pradius_disks[1];
    m_two_tendons       = two_tendons;

}

ConstantCurvatureModel::~ConstantCurvatureModel()
{
}


bool ConstantCurvatureModel::forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q)
{
    // Compute arc parameters for each segment
    double kappa1, kappa2;
    double phi1, phi2;
    double l1, l2;

    // Calculate tendon lengths
    Eigen::Vector3d l_t1;
    l_t1 <<     m_length[0] - q(0),
            m_length[0] - q(1),
            m_length[0] - q(2);

    Eigen::Vector3d l_t2;
    l_t2 <<     m_length[1] - q(3),
            m_length[1] - q(4),
            m_length[1] - q(5);


    // First segment
    
    // Check if tendon is straight
    double temp =std::sqrt(l_t1(0)*l_t1(0) + l_t1(1)*l_t1(1) + l_t1(2)*l_t1(2) - l_t1(0)*l_t1(1) - l_t1(0)*l_t1(2) - l_t1(1)*l_t1(2));

    if(m_two_tendons == false)
    {

        if(temp == 0 || std::isnan(temp))
        {
            kappa1 = 0;
            phi1 = 0;
            l1 = l_t1(0);
        }
        else
        {	
			// Compute arc parameters
            kappa1 = 2*temp/(m_pradius_disks[0]*(l_t1(0) + l_t1(1) + l_t1(2)));
            phi1 = std::atan2(sqrt(3)*(l_t1(1)+l_t1(2)-2*l_t1(0)),3*(l_t1(2)-l_t1(1)));
            l1 = ((m_number_disks[0]*m_pradius_disks[0]*(l_t1(0)+l_t1(1)+l_t1(2)))/temp)*std::asin(temp/(3*m_number_disks[0]*m_pradius_disks[0]));
        }

        // Second segment
        
		// Check if tendon is straight
        temp =std::sqrt(l_t2(0)*l_t2(0) + l_t2(1)*l_t2(1) + l_t2(2)*l_t2(2) - l_t2(0)*l_t2(1) - l_t2(0)*l_t2(2) - l_t2(1)*l_t2(2));

        if(temp == 0 || std::isnan(temp))
        {
            kappa2 = 0;
            phi2 = 0;
            l2 = l_t2(0);
        }
        else
        {
			// Compute arc parameters
            kappa2 = 2*temp/(m_pradius_disks[1]*(l_t2(0) + l_t2(1) + l_t2(2)));
            phi2 = std::atan2(sqrt(3)*(l_t2(1)+l_t2(2)-2*l_t2(0)),3*(l_t2(2)-l_t2(1)));
            l2 = ((m_number_disks[1]*m_pradius_disks[1]*(l_t2(0)+l_t2(1)+l_t2(2)))/temp)*std::asin(temp/(3*m_number_disks[1]*m_pradius_disks[1]));
        }
    }
    else
    {
        if(q(0) > 0)
        {
            phi1 = M_PI/2;
            kappa1 = q(0)/(m_pradius_disks[0]*m_length[0]);

        }
        else if(q(1) > 0)
        {
            phi1 = 3*M_PI/2;
            kappa1 = q(1)/(m_pradius_disks[0]*m_length[0]);
        }
        else
        {
            phi1 = 0;
            kappa1 = 0;
        }


        if(q(3) > 0)
        {
            phi2 = M_PI/2;
            kappa2 = q(3)/(m_pradius_disks[1]*m_length[1]);

        }
        else if(q(4) > 0)
        {
            phi2 = 3*M_PI/2;
            kappa2 = q(4)/(m_pradius_disks[1]*m_length[1]);
        }
        else
        {
            phi2 = 0;
            kappa2 = 0;
        }
    }

    //Calculate disk transformation matrices using robot independent mapping

    //First segment
    Eigen::MatrixXd disk_frames1;
    disk_frames1.resize(4,4*(m_number_disks[0]+1));

    //First disk is identity (base disk)
    disk_frames1.block(0,0,4,4) = Eigen::Matrix4d::Identity();

    Eigen::Matrix4d T;
    double c_p = std::cos(phi1);
    double s_p = std::sin(phi1);
    for(int i = 1; i <= m_number_disks[0]; i ++)
    {
        double s = i*m_length[0]/m_number_disks[0];
        double c_ks = std::cos(kappa1*s);
        double s_ks = std::sin(kappa1*s);
        if(kappa1 == 0)
        {
            T <<    c_p*c_ks, -s_p, c_p*s_ks, 0,
                    s_p*c_ks, c_p, s_p*s_ks, 0,
                    -s_ks, 0, c_ks, s,
                    0, 0, 0, 1;
        }
        else
        {
            T << c_p*c_ks, -s_p, c_p*s_ks, c_p*(1-c_ks)/kappa1,
                    s_p*c_ks, c_p, s_p*s_ks, s_p*(1-c_ks)/kappa1,
                    -s_ks, 0, c_ks, s_ks/kappa1,
                    0, 0, 0, 1;
        }
        //Apply Transformation to fix the z-axis (local z-orientation stays fixed for TDCR)
        Eigen::Matrix4d rot_phi;
        rot_phi << std::cos(phi1), std::sin(phi1), 0, 0,
                -std::sin(phi1), std::cos(phi1),0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
        T = T*rot_phi;
        disk_frames1.block(0,4*i,4,4) = T;

    }

    //Second segment
    Eigen::MatrixXd disk_frames2;
    disk_frames2.resize(4,4*m_number_disks[1]);

    c_p = std::cos(phi2);
    s_p = std::sin(phi2);
    for(int i = 1; i <= m_number_disks[1]; i ++)
    {
        double s = i*m_length[1]/m_number_disks[1];
        double c_ks = std::cos(kappa2*s);
        double s_ks = std::sin(kappa2*s);
        if(kappa2 == 0)
        {
            T <<    c_p*c_ks, -s_p, c_p*s_ks, 0,
                    s_p*c_ks, c_p, s_p*s_ks, 0,
                    -s_ks, 0, c_ks, s,
                    0, 0, 0, 1;
        }
        else
        {
            T << c_p*c_ks, -s_p, c_p*s_ks, c_p*(1-c_ks)/kappa2,
                    s_p*c_ks, c_p, s_p*s_ks, s_p*(1-c_ks)/kappa2,
                    -s_ks, 0, c_ks, s_ks/kappa2,
                    0, 0, 0, 1;
        }
        //Apply Transformation to fix the z-axis
        Eigen::Matrix4d rot_phi;
        rot_phi << std::cos(phi2), std::sin(phi2), 0, 0,
                -std::sin(phi2), std::cos(phi2),0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
        T = T*rot_phi;
        //Apply transformation of first segment
        T = disk_frames1.block(0,4*m_number_disks[0],4,4)*T;

        disk_frames2.block(0,4*(i-1),4,4) = T;

    }

    diskFrames.resize(4,disk_frames1.cols()+disk_frames2.cols());

    diskFrames << disk_frames1, disk_frames2;

    return true;

}




