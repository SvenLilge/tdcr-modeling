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

#include "tendondrivenrobot.h"



TendonDrivenRobot::TendonDrivenRobot()
{


    mp_cr_model = new CosseratRodModel();
    mp_cc_model = new ConstantCurvatureModel();
    mp_pcc_model = new PiecewiseConstantCurvatureModel();
    mp_prb_model = new PseudoRigidBodyModel();
    mp_sscr_model = new SubsegmentCosseratRodModel();

    m_keep_inits = false;

    m_length[0]         = 0.1;
    m_length[1]         = 0.1;
    m_youngs_modulus    = 60e9;
    m_number_disks[0]   = 10;
    m_number_disks[1]   = 10;
    m_pradius_disks[0]  = 0.008;
    m_pradius_disks[1]  = 0.006;
    m_ro                = 0.0005;

    std::vector<Eigen::Vector3d> routing;

    Eigen::Vector3d tendon1;
    tendon1 << 0,
            m_pradius_disks[0],
            0;

    Eigen::Vector3d tendon2;
    tendon2 <<  m_pradius_disks[0]*std::cos(-M_PI/6),
            m_pradius_disks[0]*std::sin(-M_PI/6),
            0;
    Eigen::Vector3d tendon3;
    tendon3 <<  m_pradius_disks[0]*std::cos(7*M_PI/6),
            m_pradius_disks[0]*std::sin(7*M_PI/6),
            0;

    routing.push_back(tendon1);
    routing.push_back(tendon2);
    routing.push_back(tendon3);

    tendon1 << 0,
            m_pradius_disks[1],
            0;

    tendon2 <<  m_pradius_disks[1]*std::cos(-M_PI/6),
            m_pradius_disks[1]*std::sin(-M_PI/6),
            0;

    tendon3 <<  m_pradius_disks[1]*std::cos(7*M_PI/6),
            m_pradius_disks[1]*std::sin(7*M_PI/6),
            0;

    routing.push_back(tendon1);
    routing.push_back(tendon2);
    routing.push_back(tendon3);

    m_routing           = routing;

	m_two_tendons = false;

    mp_cr_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);
    mp_cc_model->setRobotParameters(m_length,m_number_disks,m_pradius_disks,m_two_tendons);
    mp_pcc_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);
    mp_prb_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);
    mp_sscr_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);

    m_actuation_values.resize(6,1);
    m_actuation_values.setZero();
    m_ext_forces.setZero();
    m_ext_moments.setZero();


}


void TendonDrivenRobot::setRobotParameters(std::array<double, 2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int, 2> number_disks, std::array<double, 2> pradius_disks, double ro, bool two_tendons)
{
    m_length[0]         = length[0];
    m_length[1]         = length[1];
    m_youngs_modulus    = youngs_modulus;
    m_number_disks[0]   = number_disks[0];
    m_number_disks[1]   = number_disks[1];
    m_pradius_disks[0]  = pradius_disks[0];
    m_pradius_disks[1]  = pradius_disks[1];
    m_ro                = ro;
    m_routing           = routing;
	m_two_tendons		= two_tendons;


    //Update Parameters in model as well
    mp_cr_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);
    mp_cc_model->setRobotParameters(m_length,m_number_disks,m_pradius_disks,m_two_tendons);
    mp_pcc_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);
    mp_prb_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);
    mp_sscr_model->setRobotParameters(m_length,m_youngs_modulus,m_routing,m_number_disks,m_pradius_disks,m_ro);



}

TendonDrivenRobot::~TendonDrivenRobot()
{
    delete mp_cr_model;
    delete mp_cc_model;
    delete mp_pcc_model;
    delete mp_prb_model;
    delete mp_sscr_model;

}


// q is an 6x1 vector containing the forces for each tendon or displacement in m in the case of the constat curvature model. Tendons should be located with the specified pitch radius
bool TendonDrivenRobot::forwardKinematics(Eigen::Matrix4d &ee_frame, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext, Model model)
{

    //Set current actuation values
    m_actuation_values = q;
    m_ext_forces = f_ext;
    m_ext_moments = l_ext;

    Eigen::MatrixXd diskFrames;

    bool success;


    //Choose correct model
    switch(model) {
    case CosseratRod:                   success = mp_cr_model->forwardKinematics(diskFrames,q,f_ext,l_ext);
        break;

    case ConstantCurvature:             success = mp_cc_model->forwardKinematics(diskFrames,q);
        break;

    case PiecewiseConstantCurvature:    success = mp_pcc_model->forwardKinematics(diskFrames,q,f_ext,l_ext);
        break;

    case PseudoRigidBody:               success = mp_prb_model->forwardKinematics(diskFrames,q,f_ext,l_ext);
        break;

    case SubsegmentCosseratRod:         success = mp_sscr_model->forwardKinematics(diskFrames,q,f_ext,l_ext);
        break;

    default:                            success = mp_cr_model->forwardKinematics(diskFrames,q,f_ext,l_ext);
        break;
    }



    m_disk_frames = diskFrames;

    m_ee_frame = m_disk_frames.rightCols(4);

    ee_frame = m_ee_frame;


    return success;

}

void TendonDrivenRobot::keepInits(bool keep)
{
    m_keep_inits = keep;
    mp_cr_model->setKeepInits(keep);
    mp_pcc_model->setKeepInits(keep);
    mp_prb_model->setKeepInits(keep);
    mp_sscr_model->setKeepInits(keep);
}

Eigen::Matrix4d TendonDrivenRobot::getEEFrame()
{
    return m_ee_frame;
}



Eigen::MatrixXd TendonDrivenRobot::getCurrentConfig()
{
    return m_actuation_values;
}

Eigen::MatrixXd TendonDrivenRobot::getDiskFrames()
{
    return m_disk_frames;
}

Eigen::MatrixXd TendonDrivenRobot::getTendonDisplacements()
{
    Eigen::MatrixXd tendon_length;
    tendon_length.resize(6,1);
    tendon_length.setZero();

    //Iterate through all the disk frames
    for(int k = 1; k < m_disk_frames.cols()/4; k++)
    {
        Eigen::Matrix4d disk_frame;
        disk_frame = m_disk_frames.block(0,4*k,4,4);

        Eigen::Matrix4d disk_frame_prev;
        disk_frame_prev = m_disk_frames.block(0,4*(k-1),4,4);

        Eigen::Matrix4d routing1 = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d routing2 = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d routing3 = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d routing4 = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d routing5 = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d routing6 = Eigen::Matrix4d::Identity();

        routing1.block(0,3,3,1) = m_routing.at(0);
        routing2.block(0,3,3,1) = m_routing.at(1);
        routing3.block(0,3,3,1) = m_routing.at(2);
        routing4.block(0,3,3,1) = m_routing.at(3);
        routing5.block(0,3,3,1) = m_routing.at(4);
        routing6.block(0,3,3,1) = m_routing.at(5);

        if(k < m_number_disks[0] + 1)
        {
            tendon_length(0) += ((disk_frame*routing1).block(0,3,3,1) - (disk_frame_prev*routing1).block(0,3,3,1)).norm();
            tendon_length(1) += ((disk_frame*routing2).block(0,3,3,1) - (disk_frame_prev*routing2).block(0,3,3,1)).norm();
            tendon_length(2) += ((disk_frame*routing3).block(0,3,3,1) - (disk_frame_prev*routing3).block(0,3,3,1)).norm();
        }
        else
        {
            tendon_length(3) += ((disk_frame*routing4).block(0,3,3,1) - (disk_frame_prev*routing4).block(0,3,3,1)).norm();
            tendon_length(4) += ((disk_frame*routing5).block(0,3,3,1) - (disk_frame_prev*routing5).block(0,3,3,1)).norm();
            tendon_length(5) += ((disk_frame*routing6).block(0,3,3,1) - (disk_frame_prev*routing6).block(0,3,3,1)).norm();
        }
    }


    tendon_length(0) = m_length[0] - tendon_length(0);
    tendon_length(1) = m_length[0] - tendon_length(1);
    tendon_length(2) = m_length[0] - tendon_length(2);

    tendon_length(3) = m_length[1] - tendon_length(3);
    tendon_length(4) = m_length[1] - tendon_length(4);
    tendon_length(5) = m_length[1] - tendon_length(5);

    return tendon_length;
}








