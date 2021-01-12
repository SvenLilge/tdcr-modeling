#ifndef SUBSEGMENTCOSSERATRODMODEL_H
#define SUBSEGMENTCOSSERATRODMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <array>
#include <gsl/gsl_odeiv2.h>
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

class SubsegmentCosseratRodModel 
{
public:
    SubsegmentCosseratRodModel();
    virtual ~SubsegmentCosseratRodModel();

    bool forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext);
    void cosserat_ode(Eigen::MatrixXd &dy, Eigen::MatrixXd y);
    void get_res(Eigen::MatrixXd &output, Eigen::MatrixXd input);
    void setKeepInits(bool keep);
    void setRobotParameters(std::array<double,2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int,2> number_disks, std::array<double,2> pradius_disks, double ro);
    double getNumberOfTotalDisks();
    void setDefaultInitValues(Eigen::MatrixXd inits);
    Eigen::MatrixXd getFinalInitValues();
    double getFinalResudial();


private:
    void run_IVP(Eigen::MatrixXd &states, Eigen::MatrixXd init_state);


    std::array<double,2> m_length;
    double m_youngs_modulus;
    std::vector<Eigen::Vector3d> m_routing;
    std::array<int,2> m_number_disks;
    std::array<double,2> m_pradius_disks;
    double m_ro;

    bool m_keep_inits;

    Eigen::Matrix3d m_Kbt;
    Eigen::Matrix3d m_Kse;

    Eigen::Vector3d m_f_ext;
    Eigen::Vector3d m_l_ext;

    Eigen::Matrix<double,6,1> m_tau;

    Eigen::MatrixXd m_last_inits;
    Eigen::MatrixXd m_default_inits;

    double m_residual_error;

    Eigen::MatrixXd m_states;

    Eigen::VectorXd m_length_ss;

    int m_current_segment;

    //Stuff for ODE solving
    gsl_odeiv2_system m_ode_sys;
    const gsl_odeiv2_step_type* m_ode_step_type;
    gsl_odeiv2_step* m_ode_step;
    gsl_odeiv2_control* m_ode_controller;
    gsl_odeiv2_evolve* m_ode_evolve;
};

#endif // SUBSEGMENTCOSSERATRODMODEL_H
