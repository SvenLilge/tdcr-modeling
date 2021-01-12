//includes
#include "tendondrivenrobot.h"

//stl
#include <iostream>

//Eigen
#include <Eigen/Dense>


int main(int argc, char **argv)
{
	//Create tendon driven robot with default parameters
    TendonDrivenRobot robot;
    
    //Example on defining custom parameters for the robot
    
    //Length per segment
    std::array<double,2> length;
    length[0] = 0.2;
    length[1] = 0.2;

	//Youngs modulus
    double youngs_modulus = 54*1e9;

	//Number of disks per segment
    std::array<int,2> number_disks;
    number_disks[0] = 10;
    number_disks[1] = 10;
	
	//Pitch radius per segment (distance from tendons to backbone)
    std::array<double,2> pradius_disks;
    pradius_disks[0] = 10*1e-3;
    pradius_disks[1] = 10*1e-3;
	
	
	//Radius of the backbone
    double ro = 0.7*1e-3;
	
	
	
	//Tendon routing expressed as a three-dimensional position vector in the local (disk) frame of the robot for each tendon
	//First three tendons for first segment, last three tendons for second segment
	//Stored in a std::vector
    std::vector<Eigen::Vector3d> routing;

    Eigen::Vector3d tendon1;
    tendon1 << 0,
            pradius_disks[0],
            0;

    Eigen::Vector3d tendon2;
    tendon2 <<  pradius_disks[0]*std::cos(-M_PI/6),
            pradius_disks[0]*std::sin(-M_PI/6),
            0;
    Eigen::Vector3d tendon3;
    tendon3 <<  pradius_disks[0]*std::cos(7*M_PI/6),
            pradius_disks[0]*std::sin(7*M_PI/6),
            0;

    routing.push_back(tendon1);
    routing.push_back(tendon2);
    routing.push_back(tendon3);

    tendon1 << 0,
            pradius_disks[1],
            0;

    tendon2 <<  pradius_disks[1]*std::cos(-M_PI/6),
            pradius_disks[1]*std::sin(-M_PI/6),
            0;

    tendon3 <<  pradius_disks[1]*std::cos(7*M_PI/6),
            pradius_disks[1]*std::sin(7*M_PI/6),
            0;

    routing.push_back(tendon1);
    routing.push_back(tendon2);
    routing.push_back(tendon3);

    //Update the robot parameters
    robot.setRobotParameters(length,youngs_modulus,routing,number_disks,pradius_disks,ro,false);
	
	//Defining the actuation values
	//Every entry in q defined actuation of one tendon in the order used above
	//q is tendon force in Newton for static models tendon displacement in meter for Constant Curvature model
	Eigen::MatrixXd q;
    q.resize(6,1);
    q << 4,
		 0,
		 0,
		 0,
		 2,
		 0;
	
	//External force at tip set to zero
    Eigen::Vector3d f_ext;
    f_ext.setZero();

	//External monment at tip set to zero
    Eigen::Vector3d l_ext;
    l_ext.setZero();

    int success;
	
	//Variable to store the tip frame
    Eigen::Matrix4d ee_frame;

	//Running the forward kinematics using the different models
	success = robot.forwardKinematics(ee_frame,q,f_ext,l_ext,TendonDrivenRobot::Model::CosseratRod);
	if(success)
		std::cout << "Tip frame calculated using Cosserat Rod:" <<std::endl << ee_frame <<std::endl;
	
	success = robot.forwardKinematics(ee_frame,q,f_ext,l_ext,TendonDrivenRobot::Model::PiecewiseConstantCurvature);
	if(success)
		std::cout << "Tip frame calculated using Piecewise Constant Curvature:" <<std::endl << ee_frame <<std::endl;
			
	success = robot.forwardKinematics(ee_frame,q,f_ext,l_ext,TendonDrivenRobot::Model::PseudoRigidBody);
	if(success)
		std::cout << "Tip frame calculated using Pseudo Rigid Body:" <<std::endl << ee_frame <<std::endl;
	
	success = robot.forwardKinematics(ee_frame,q,f_ext,l_ext,TendonDrivenRobot::Model::SubsegmentCosseratRod);
	if(success)
		std::cout << "Tip frame calculated using Subsegment Cosserat Rod:" <<std::endl << ee_frame <<std::endl;
		
	//Use the current state of the robot (calculated with the Subsegment Cosserat Rod model, since this was last called) to calculate the current tendon displacements as an input for the Constant Curvature Model
	//Note that this is an approximation of the tendon displacements based on another model
	q = robot.getTendonDisplacements();
	
	success = robot.forwardKinematics(ee_frame,q,f_ext,l_ext,TendonDrivenRobot::Model::ConstantCurvature);	
	if(success)
		std::cout << "Tip frame calculated using Constant Curvature:" <<std::endl << ee_frame <<std::endl;
		
	return 1;
}
