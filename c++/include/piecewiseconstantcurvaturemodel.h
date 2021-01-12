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

class PiecewiseConstantCurvatureModel 
{
public:
    PiecewiseConstantCurvatureModel();
    virtual ~PiecewiseConstantCurvatureModel();

    bool forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext);
    void get_res(Eigen::MatrixXd &output, Eigen::MatrixXd input);
    void setKeepInits(bool keep);
    void setRobotParameters(std::array<double,2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int,2> number_disks, std::array<double,2> pradius_disks, double ro);
    double getNumberOfTotalDisks();
    void setDefaultInitValues(Eigen::MatrixXd inits);
    Eigen::MatrixXd getFinalInitValues();
    double getFinalResudial();

private:

    Eigen::Matrix4d getDiskTransform(Eigen::MatrixXd var, Eigen::VectorXd length_ss, int idx_from, int idx_to);

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
