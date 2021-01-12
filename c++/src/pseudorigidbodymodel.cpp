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

#include "pseudorigidbodymodel.h"

int evaluate_boundary_conditions_prb(const gsl_vector * x, void *param, gsl_vector * f)
{


    PseudoRigidBodyModel *model = (PseudoRigidBodyModel*)param;

    Eigen::MatrixXd inputs;

    inputs.resize((model->getNumberOfJoints()+1)*model->getNumberOfTotalDisks(),1);

    Eigen::MatrixXd outputs;

    for(int i = 0; i < (model->getNumberOfJoints()+1)*model->getNumberOfTotalDisks(); i++)
    {
        inputs(i,0) = gsl_vector_get (x, i);
    }


    //Function that needs to be minimized goes here

    model->get_res(outputs, inputs);



    //Functions that needs to be minimized ends here

    for(int i = 0; i < (model->getNumberOfJoints()+1)*model->getNumberOfTotalDisks(); i++)
    {

        gsl_vector_set (f, i, outputs(i,0));
    }


    return GSL_SUCCESS;
}

void
callback_prb(const size_t iter, void *param,
         const gsl_multifit_nlinear_workspace *w)
{

    PseudoRigidBodyModel *model = (PseudoRigidBodyModel*)param;

    //Pring current iteration
    std::cout << "Iteration " << iter << std::endl;

    gsl_vector * x = gsl_multifit_nlinear_position(w);
    gsl_vector * r = gsl_multifit_nlinear_residual(w);

    std::cout << "Current values:" << std::endl;
    for(int i = 0; i < (model->getNumberOfJoints()+1)*model->getNumberOfTotalDisks(); i++)
    {
        std::cout << " " << gsl_vector_get(x, i);
    }
    std::cout << std::endl << "Current residuals:" << std::endl;
    for(int i = 0; i < (model->getNumberOfJoints()+1)*model->getNumberOfTotalDisks(); i++)
    {
        std::cout << " " << gsl_vector_get(r, i);
    }

    double chisq;
    gsl_blas_ddot(r, r, &chisq);

    std::cout << std::endl << "Current cost: "<< chisq << std::endl << std::endl;

}

PseudoRigidBodyModel::PseudoRigidBodyModel()
{
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






    //4 joints per subsegment with relative positions/lengths
    m_joint_pos.resize(4);
    m_joint_pos <<  0.125,
                    0.35,
                    0.388,
                    0.136;

    m_joint_pos = m_joint_pos/m_joint_pos.sum();

    m_keep_inits = false;
    m_default_inits.resize(1,(m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]));
    for(int i = 0; i < m_default_inits.cols(); i++)
    {
        m_default_inits(0,i) = 0;
    }
    m_last_inits = m_default_inits;



    m_tau.setZero();

}

void PseudoRigidBodyModel::setRobotParameters(std::array<double, 2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int, 2> number_disks, std::array<double, 2> pradius_disks, double ro)
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


}

double PseudoRigidBodyModel::getNumberOfTotalDisks()
{
    return m_number_disks[0] + m_number_disks[1];
}

double PseudoRigidBodyModel::getNumberOfJoints()
{
    return m_joint_pos.size();
}

void PseudoRigidBodyModel::setDefaultInitValues(Eigen::MatrixXd inits)
{
    m_default_inits = inits;
}

Eigen::MatrixXd PseudoRigidBodyModel::getFinalInitValues()
{
    return m_last_inits;
}

double PseudoRigidBodyModel::getFinalResudial()
{
    return m_residual_error;
}

Eigen::Matrix4d PseudoRigidBodyModel::getDiskTransform(Eigen::MatrixXd &Trb, Eigen::MatrixXd var, Eigen::VectorXd length_ss, int idx_from, int idx_to)
{
    Eigen::Matrix4d diskTransform = Eigen::Matrix4d::Identity();

    int nrb = m_joint_pos.size();

    Trb.resize(4,4*nrb*(idx_from-idx_to));


    for(int i = idx_to + 1; i <= idx_from; i++)
    {
        Eigen::MatrixXd theta;
        theta = var.block((nrb+1)*i-(nrb+1),0,nrb-1,1);
        double phi = var((nrb+1)*i-(nrb+1) + 3,0);
        double epsi = var((nrb+1)*i-(nrb+1) + 4,0);

        Eigen::Matrix4d Rphi;
        Rphi << std::cos(phi), -1*std::sin(phi), 0, 0,
              std::sin(phi), std::cos(phi), 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;

        Eigen::Matrix4d Rphi_epsi;
        Rphi_epsi << std::cos(epsi-phi), -1*std::sin(epsi-phi), 0, 0,
              std::sin(epsi-phi), std::cos(epsi-phi), 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;


        Eigen::Matrix4d gamma_mat;
        gamma_mat <<    1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, m_joint_pos(0)*length_ss(i-1),
                        0, 0, 0, 1;

        Eigen::Matrix4d tempTransform;
        tempTransform = Rphi*gamma_mat;
        Trb.block(0,4*(nrb*(i-(idx_to+1))),4,4) = diskTransform*tempTransform;

        for(int k = 0; k < nrb-1; k++)
        {
            Eigen::Matrix4d theta_mat;
            theta_mat << std::cos(theta(k)), 0, std::sin(theta(k)), m_joint_pos(k+1)*length_ss(i-1)*std::sin(theta(k)),
                         0, 1, 0, 0,
                         -1*std::sin(theta(k)), 0, std::cos(theta(k)), m_joint_pos(k+1)*length_ss(i-1)*std::cos(theta(k)),
                         0, 0, 0, 1;
            tempTransform = tempTransform*theta_mat;
            Trb.block(0,4*(nrb*(i-(idx_to+1))+k+1),4,4) = diskTransform*tempTransform;
        }
        diskTransform = diskTransform*tempTransform*Rphi_epsi;

    }

    return diskTransform;

}

PseudoRigidBodyModel::~PseudoRigidBodyModel()
{

}


bool PseudoRigidBodyModel::forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext)
{
    m_tau = q;
    m_f_ext = f_ext;
    m_l_ext = l_ext;

    // Run BVP Shooting method
    Eigen::MatrixXd inits;

    if(m_keep_inits)
    {
        inits = m_last_inits;
    }
    else
    {
        inits = m_default_inits;
    }

    const size_t n = (m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]); //Number of equations (boundary conditions)
    const size_t p = (m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]); //Number of parameters (initial guesses)
    gsl_vector *x0 = gsl_vector_alloc(p);
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    fdf.f = evaluate_boundary_conditions_prb;
    fdf.df = NULL; //Provide your own Jacobian
    fdf.fvv = NULL; //Provide your own Hessian
    fdf.n = n;
    fdf.p = p;
    fdf.params = this;

    //Set initial values
    for(int i = 0; i < (m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]); i++)
    {
        gsl_vector_set(x0,i,inits(0,i));
    }

    //Choose solver (Levenberg-Marquardt)
    fdf_params.trs = gsl_multifit_nlinear_trs_lm;
    //Choose solver (Levenberg-Marquardt with acceleration)
    //fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
    //Choose solver (dogleg)
    //fdf_params.trs =  gsl_multifit_nlinear_trs_dogleg;

    //Set scaling parameter
    fdf_params.scale =   gsl_multifit_nlinear_scale_more;

    //Set solver for trust region
    fdf_params.solver =  gsl_multifit_nlinear_solver_qr;

    //Set trust region factor
    fdf_params.factor_up = 2.5; //default 3
    fdf_params.factor_down = 4; //default 2


    ///---Solve the system---
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    const size_t max_iter = 1000;
    const double xtol = 1.0e-10;
    const double gtol = 1.0e-10;
    const double ftol = 1.0e-10;
    gsl_multifit_nlinear_workspace *work = gsl_multifit_nlinear_alloc(T,&fdf_params, n, p);
    gsl_vector * f = gsl_multifit_nlinear_residual(work);
    gsl_vector * x = gsl_multifit_nlinear_position(work);
    int info;
    double chisq0, chisq, rcond;

    /* initialize solver */
    gsl_multifit_nlinear_init(x0, &fdf, work);

    /* store initial cost */
    gsl_blas_ddot(f, f, &chisq0);

    /* iterate until convergence */
    //gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,callback_prb, this, &info, work);
    //Without Callback
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,NULL, NULL, &info, work);

    /* store final cost */
    gsl_blas_ddot(f, f, &chisq);
    m_residual_error = chisq;

    /* store cond(J(x)) */
    gsl_multifit_nlinear_rcond(&rcond, work);

    ///---End Solve the system---

    Eigen::MatrixXd final_values;
    final_values.resize(1,(m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]));

    Eigen::MatrixXd final_boundcon;
    final_boundcon.resize(1,(m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]));


    for(int i = 0; i < (m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]); i ++)
    {
        final_values(0,i) = gsl_vector_get (x, i);
        final_boundcon(0,i) = gsl_vector_get (f, i);
    }

    m_last_inits = final_values;
    gsl_vector_free (x0);
    gsl_multifit_nlinear_free(work);




    Eigen::VectorXd length_ss;
    length_ss.resize(m_number_disks[0]+m_number_disks[1]);

    for(int i = 0; i < m_number_disks[0]; i++)
    {
        length_ss(i) = m_length[0]/m_number_disks[0];
    }

    for(int i = m_number_disks[0]; i < m_number_disks[0]+m_number_disks[1]; i++)
    {
        length_ss(i) = m_length[1]/m_number_disks[1];
    }

    diskFrames.resize(4,4*(m_number_disks[0]+m_number_disks[1]+1));

    diskFrames.block(0,0,4,4) = Eigen::Matrix4d::Identity();

    Eigen::MatrixXd Trb;

    for(int i = 1; i < m_number_disks[0]+m_number_disks[1] + 1; i++)
    {
        Eigen::Matrix4d frame = getDiskTransform(Trb, final_values.transpose(),length_ss,i,0);
        diskFrames.block(0,4*i,4,4) = frame;
    }


    if(m_residual_error > 1e-6)
    {
        return false;
    }


    return true;

}

void PseudoRigidBodyModel::get_res(Eigen::MatrixXd &output, Eigen::MatrixXd input)
{

    Eigen::VectorXd length_ss;
    length_ss.resize(m_number_disks[0]+m_number_disks[1]);

    for(int i = 0; i < m_number_disks[0]; i++)
    {
        length_ss(i) = m_length[0]/m_number_disks[0];
    }

    for(int i = m_number_disks[0]; i < m_number_disks[0]+m_number_disks[1]; i++)
    {
        length_ss(i) = m_length[1]/m_number_disks[1];
    }

    output.resize((m_joint_pos.size()+1)*(m_number_disks[0]+m_number_disks[1]),1);

    int nrb = m_joint_pos.size();

    Eigen::Vector4d F_prev;
    F_prev.setZero();
    Eigen::Vector3d M_prev;
    M_prev.setZero();

    for(int ss_i = m_number_disks[0]+m_number_disks[1]; ss_i > 0; ss_i--)
    {

        //Get transformation matrix between disk i and disk i-1
        Eigen::Matrix4d diskTransform;
        Eigen::MatrixXd Trb;
        diskTransform = getDiskTransform(Trb,input,length_ss,ss_i,ss_i-1);

        Eigen::MatrixXd theta;
        theta = input.block((nrb+1)*ss_i-(nrb+1),0,nrb-1,1);
        double phi = input((nrb+1)*ss_i-(nrb+1) + 3,0);
        double epsi = input((nrb+1)*ss_i-(nrb+1) + 4,0);
        Eigen::Vector3d ni(std::cos(phi+M_PI/2),std::sin(phi+M_PI/2),0);

        Eigen::Matrix<double,3,6> tendon_pos;
        for(int i = 0; i < 6; i++)
        {
            tendon_pos.col(i) = diskTransform.block(0,0,3,3)*m_routing.at(i) + diskTransform.block(0,3,3,1);
        }

        Eigen::MatrixXd pt_mat;
        pt_mat.resize(4,6);
        pt_mat << m_routing.at(0) - tendon_pos.col(0), m_routing.at(1) - tendon_pos.col(1), m_routing.at(2) - tendon_pos.col(2), m_routing.at(3) - tendon_pos.col(3), m_routing.at(4) - tendon_pos.col(4), m_routing.at(5) - tendon_pos.col(5),
                  1, 1, 1, 1, 1, 1;

        Eigen::VectorXd norm_ct1;
        norm_ct1.resize(6);
        for(int i = 0; i < 6; i++)
        {
            norm_ct1(i) = (m_routing.at(i) - tendon_pos.col(i)).norm();
        }

        Eigen::MatrixXd ct1_mat;
        ct1_mat.resize(4,6);
        ct1_mat << norm_ct1.transpose(),
                   norm_ct1.transpose(),
                   norm_ct1.transpose(),
                   norm_ct1.transpose();

        Eigen::MatrixXd Pi;
        Pi.resize(4,nrb);
        for(int k = 0; k < nrb; k++)
        {
            Pi.col(k) = Trb.block(0,4*k+3,4,1);
        }

        Eigen::MatrixXd F_disk;
        F_disk.resize(m_number_disks[0]+m_number_disks[1],6);
        F_disk.setZero();

        for(int i = 0; i < m_number_disks[0]; i++)
        {
            F_disk.block(i,0,1,6) = m_tau.transpose();
        }

        for(int i = m_number_disks[0]; i < m_number_disks[0] + m_number_disks[1]; i++)
        {
            F_disk.block(i,3,1,3) = m_tau.block(3,0,3,1).transpose();
        }

        Eigen::MatrixXd F_disk_mat;
        F_disk_mat.resize(4,6);
        F_disk_mat << F_disk.row(ss_i-1),
                      F_disk.row(ss_i-1),
                      F_disk.row(ss_i-1),
                      F_disk.row(ss_i-1);

        Eigen::MatrixXd Fi;
        Fi.resize(4,6);

        Eigen::Vector4d zi = diskTransform.col(2);


        Eigen::Matrix4d diskTransform1;
        Eigen::MatrixXd Trb1;

        if(ss_i < m_number_disks[0]+m_number_disks[1])
        {
            diskTransform1 = getDiskTransform(Trb1,input,length_ss,ss_i+1,ss_i-1);

            Eigen::Matrix<double,3,6> tendon_pos1;
            for(int i = 0; i < 6; i++)
            {
                tendon_pos1.col(i) = diskTransform1.block(0,0,3,3)*m_routing.at(i) + diskTransform1.block(0,3,3,1);
            }

            Eigen::MatrixXd pt_mat1;
            pt_mat1.resize(4,6);
            pt_mat1 << tendon_pos1.col(0) - tendon_pos.col(0), tendon_pos1.col(1) - tendon_pos.col(1), tendon_pos1.col(2) - tendon_pos.col(2), tendon_pos1.col(3) - tendon_pos.col(3), tendon_pos1.col(4) - tendon_pos.col(4), tendon_pos1.col(5) - tendon_pos.col(5),
                      1, 1, 1, 1, 1, 1;

            Eigen::VectorXd norm_ct2;
            norm_ct2.resize(6);
            for(int i = 0; i < 6; i++)
            {
                norm_ct2(i) = (tendon_pos1.col(i) - tendon_pos.col(i)).norm();
            }

            Eigen::MatrixXd ct2_mat;
            ct2_mat.resize(4,6);
            ct2_mat << norm_ct2.transpose(),
                       norm_ct2.transpose(),
                       norm_ct2.transpose(),
                       norm_ct2.transpose();

            Eigen::MatrixXd F_disk_mat1;
            F_disk_mat1.resize(4,6);
            F_disk_mat1 << F_disk.row(ss_i),
                          F_disk.row(ss_i),
                          F_disk.row(ss_i),
                          F_disk.row(ss_i);

            Fi = pt_mat.array()/ct1_mat.array()*F_disk_mat.array() + pt_mat1.array()/ct2_mat.array()*F_disk_mat1.array();
            if(ss_i == m_number_disks[0])
            {
                for(int i = 3; i < 6; i++)
                {
                    Fi.col(i) = Fi.col(i) - zi.dot(Fi.col(i))*zi/zi.norm()/zi.norm();
                }

            }
            else
            {
                for(int i = 0; i < 6; i++)
                {
                    Fi.col(i) = Fi.col(i) - zi.dot(Fi.col(i))*zi/zi.norm()/zi.norm();
                }
            }

        }
        else
        {

            Fi = pt_mat.array()/ct1_mat.array()*F_disk_mat.array();
            F_prev.setZero();
            M_prev.setZero();
        }


        Eigen::MatrixXd Mi;
        Mi.resize(3,6);

        for(int i = 0; i < 6; i++)
        {
            Eigen::Vector3d pos = Pi.block(0,nrb-1,3,1) - tendon_pos.col(i);
            Eigen::Vector3d F = Fi.block(0,i,3,1);
            Mi.col(i) = -1*pos.cross(F);
        }

        Eigen::Vector4d F_ext;
        F_ext.setZero();
        Eigen::Vector3d M_ext;
        M_ext.setZero();

        Eigen::MatrixXd temp;
        Eigen::Matrix4d Rt = getDiskTransform(temp,input,length_ss,ss_i-1,0);
        F_ext = Rt.inverse()*Eigen::Vector4d(m_f_ext(0),m_f_ext(1),m_f_ext(2),0);
        Eigen::Matrix4d R_ex = getDiskTransform(temp,input,length_ss,m_number_disks[0]+m_number_disks[1],ss_i-1);
        Eigen::Vector3d p_ex = R_ex.block(0,3,3,1);

        Eigen::Vector3d pos = Pi.block(0,nrb-1,3,1) - p_ex;
        Eigen::Vector3d F_ext_temp = F_ext.block(0,0,3,1);

        M_ext = Rt.block(0,0,3,3).transpose()*m_l_ext - pos.cross(F_ext_temp);

        F_prev = diskTransform*F_prev;
        M_prev = diskTransform.block(0,0,3,3)*M_prev;

        Eigen::Vector4d F_net = Fi.rowwise().sum() + F_prev;
        F_net(3) = 0;

        Eigen::Vector3d pos1 = diskTransform1.block(0,3,3,1) - Pi.block(0,3,3,1);
        Eigen::Vector3d F_temp1 = F_prev.block(0,0,3,1);
        Eigen::Vector3d M_net = Mi.rowwise().sum() + M_prev + pos1.cross(F_temp1);

        F_prev = F_net;
        M_prev = M_net;

        F_net += F_ext;
        F_net(3) = 0;
        M_net += M_ext;


        double I = 0.25*M_PI*(m_ro*m_ro*m_ro*m_ro);
        double stiffness = m_youngs_modulus*I/length_ss(ss_i-1);

        double G = m_youngs_modulus/(2*1.3);

        Eigen::Vector3d K(3.25*stiffness,2.84*stiffness,2.95*stiffness);

        Eigen::Vector3d F_temp = F_net.block(0,0,3,1);
        for(int k = 0; k < nrb - 1; k++)
        {
            Eigen::Vector3d pos_temp = Pi.block(0,nrb-1,3,1) - Pi.block(0,k,3,1);
            Eigen::Matrix3d Rb = Trb.block(0,4*(k+1),3,3);
            Eigen::Vector3d M_netb = Rb.transpose()*(pos_temp.cross(F_temp)+M_net);
            output((nrb+1)*ss_i-(nrb+1)+k) = K(k)*theta(k) - M_netb(1);
        }

        Eigen::Vector3d posi = diskTransform.block(0,3,3,1);
        // Net moment applied at the base of the subsegment
        Eigen::Vector3d M_net2 = posi.cross(F_temp)+M_net;
        // Moment used for the comptation of phi
        Eigen::Vector3d M_phi = M_net2;
        M_phi(2) = 0;

        output((nrb+1)*ss_i-(nrb+1) + 3) = ni.transpose()*M_phi - M_phi.norm(); //perpendicularity condition

        // Moment used for the computation of epsi
        Eigen::Matrix3d Ri = diskTransform.block(0,0,3,3);
        Eigen::Vector3d M_epsi = Ri.transpose()*M_net2;

        output((nrb+1)*ss_i-(nrb+1) + 4) = M_epsi(2)-2*I*G/length_ss(ss_i-1)*epsi;





    }



}

void PseudoRigidBodyModel::setKeepInits(bool keep)
{
    m_keep_inits = keep;
}


