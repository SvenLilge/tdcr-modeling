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
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>


#include "cosseratrodmodel.h"
#include "constantcurvaturemodel.h"
#include "piecewiseconstantcurvaturemodel.h"
#include "pseudorigidbodymodel.h"
#include "subsegmentcosseratrodmodel.h"


class TendonDrivenRobot
{
public:
    TendonDrivenRobot();
    ~TendonDrivenRobot();

    enum Model { CosseratRod, ConstantCurvature, PiecewiseConstantCurvature, PseudoRigidBody, SubsegmentCosseratRod };

    bool forwardKinematics(Eigen::Matrix4d &ee_frame, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext, Model model);
    void setRobotParameters(std::array<double,2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int,2> number_disks, std::array<double,2> pradius_disks, double ro, bool two_tendons);
    void keepInits(bool);
    Eigen::Matrix4d getEEFrame();
    Eigen::MatrixXd getCurrentConfig();
    Eigen::MatrixXd getDiskFrames();
    Eigen::MatrixXd getTendonDisplacements();

private:
    CosseratRodModel* mp_cr_model;
    ConstantCurvatureModel* mp_cc_model;
    PiecewiseConstantCurvatureModel* mp_pcc_model;
    PseudoRigidBodyModel* mp_prb_model;
    SubsegmentCosseratRodModel* mp_sscr_model;

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
