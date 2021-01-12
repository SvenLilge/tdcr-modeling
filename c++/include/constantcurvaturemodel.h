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

class ConstantCurvatureModel 
{
public:
    ConstantCurvatureModel();
    virtual ~ConstantCurvatureModel();

    bool forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q);
    void setRobotParameters(std::array<double,2> length, std::array<int,2> number_disks, std::array<double,2> pradius_disks, bool two_tendons);

private:

    std::array<double,2> m_length;
    std::vector<Eigen::Vector3d> m_routing;
    std::array<int,2> m_number_disks;
    std::array<double,2> m_pradius_disks;
    bool m_two_tendons;
};

#endif // CONSTANTCURVATUREMODEL_H
