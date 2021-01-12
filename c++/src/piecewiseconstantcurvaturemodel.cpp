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

#include "piecewiseconstantcurvaturemodel.h"

// Wrapper function for non-linear boundary conditions for GNU::GSL
int evaluate_boundary_conditions_pcc(const gsl_vector * x, void *param, gsl_vector * f)
{


    PiecewiseConstantCurvatureModel *model = (PiecewiseConstantCurvatureModel*)param;

    Eigen::MatrixXd inputs;

    inputs.resize(3*model->getNumberOfTotalDisks(),1);

    Eigen::MatrixXd outputs;

    for(int i = 0; i < 3*model->getNumberOfTotalDisks(); i++)
    {
        inputs(i,0) = gsl_vector_get (x, i);
    }



    model->get_res(outputs, inputs);




    for(int i = 0; i < 3*model->getNumberOfTotalDisks(); i++)
    {

        gsl_vector_set (f, i, outputs(i,0));
    }


    return GSL_SUCCESS;
}

// Callback function that can be used during non-linear squares solving
void
callback_pcc(const size_t iter, void *param,
         const gsl_multifit_nlinear_workspace *w)
{

    PiecewiseConstantCurvatureModel *model = (PiecewiseConstantCurvatureModel*)param;

    //Print current iteration
    std::cout << "Iteration " << iter << std::endl;

    gsl_vector * x = gsl_multifit_nlinear_position(w);
    gsl_vector * r = gsl_multifit_nlinear_residual(w);

    std::cout << "Current values:" << std::endl;
    for(int i = 0; i < 3*model->getNumberOfTotalDisks(); i++)
    {
        std::cout << " " << gsl_vector_get(x, i);
    }
    std::cout << std::endl << "Current residuals:" << std::endl;
    for(int i = 0; i < 3*model->getNumberOfTotalDisks(); i++)
    {
        std::cout << " " << gsl_vector_get(r, i);
    }

    double chisq;
    gsl_blas_ddot(r, r, &chisq);

    std::cout << std::endl << "Current cost: "<< chisq << std::endl << std::endl;

}

PiecewiseConstantCurvatureModel::PiecewiseConstantCurvatureModel()
{
	//Initialize member variables and set up default parameters for TDCR
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



    m_keep_inits = false;
    m_default_inits.resize(1,3*(m_number_disks[0]+m_number_disks[1]));
    for(int i = 0; i < m_default_inits.cols(); i++)
    {
        m_default_inits(0,i) = 0.01;
    }
    m_last_inits = m_default_inits;


    m_tau.setZero();

}

void PiecewiseConstantCurvatureModel::setRobotParameters(std::array<double, 2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int, 2> number_disks, std::array<double, 2> pradius_disks, double ro)
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

double PiecewiseConstantCurvatureModel::getNumberOfTotalDisks()
{
    return m_number_disks[0] + m_number_disks[1];
}

void PiecewiseConstantCurvatureModel::setDefaultInitValues(Eigen::MatrixXd inits)
{
    m_default_inits = inits;
}

Eigen::MatrixXd PiecewiseConstantCurvatureModel::getFinalInitValues()
{
    return m_last_inits;
}

double PiecewiseConstantCurvatureModel::getFinalResudial()
{
    return m_residual_error;
}

Eigen::Matrix4d PiecewiseConstantCurvatureModel::getDiskTransform(Eigen::MatrixXd var, Eigen::VectorXd length_ss, int idx_from, int idx_to)
{
    Eigen::Matrix4d diskTransform = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d tempTransform = Eigen::Matrix4d::Identity();


    for(int i = idx_to + 1; i <= idx_from; i++)
    {
        double beta = var(3*i-3);
        double gamma = var(3*i-2);
        double epsi = var(3*i-1);
        double k = std::sqrt(beta*beta+gamma*gamma);
        double phi = std::atan2(gamma,beta);
        double theta = k*length_ss(i-1);
        Eigen::Vector3d p;
        p <<    (1-std::cos(theta))*std::cos(phi)/k,
                (1-std::cos(theta))*std::sin(phi)/k,
                std::sin(theta)/k;
        Eigen::Matrix3d Rz;
        Rz << std::cos(phi), -1*std::sin(phi), 0,
              std::sin(phi), std::cos(phi), 0,
              0, 0, 1;
        Eigen::Matrix3d Ry;
        Ry << std::cos(theta), 0, std::sin(theta),
              0, 1, 0,
              -1*std::sin(theta), 0, std::cos(theta);
        Eigen::Matrix3d Rz2;
        Rz2 << std::cos(epsi-phi), -1*std::sin(epsi-phi), 0,
              std::sin(epsi-phi), std::cos(epsi-phi), 0,
              0, 0, 1;
        tempTransform << Rz*Ry*Rz2, p,
                         0, 0, 0, 1;

        diskTransform = diskTransform*tempTransform;

    }

    return diskTransform;

}

PiecewiseConstantCurvatureModel::~PiecewiseConstantCurvatureModel()
{

}


bool PiecewiseConstantCurvatureModel::forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext)
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

    const size_t n = 3*(m_number_disks[0]+m_number_disks[1]); //Number of equations (boundary conditions)
    const size_t p = 3*(m_number_disks[0]+m_number_disks[1]); //Number of parameters (initial guesses)
    gsl_vector *x0 = gsl_vector_alloc(p);
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    fdf.f = evaluate_boundary_conditions_pcc;
    fdf.df = NULL; //Provide your own Jacobian
    fdf.fvv = NULL; //Provide your own Hessian
    fdf.n = n;
    fdf.p = p;
    fdf.params = this;

    //Set initial values
    for(int i = 0; i < 3*(m_number_disks[0]+m_number_disks[1]); i++)
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
    //gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,callback_pcc, this, &info, work);
    //Without Callback
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,NULL, NULL, &info, work);

    /* store final cost */
    gsl_blas_ddot(f, f, &chisq);
    m_residual_error = chisq;

    /* store cond(J(x)) */
    gsl_multifit_nlinear_rcond(&rcond, work);

    ///---End Solve the system---

    Eigen::MatrixXd final_values;
    final_values.resize(1,3*(m_number_disks[0]+m_number_disks[1]));

    Eigen::MatrixXd final_boundcon;
    final_boundcon.resize(1,3*(m_number_disks[0]+m_number_disks[1]));


    for(int i = 0; i < 3*(m_number_disks[0]+m_number_disks[1]); i ++)
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

    for(int i = 1; i < m_number_disks[0]+m_number_disks[1] + 1; i++)
    {
        Eigen::Matrix4d frame = getDiskTransform(final_values.transpose(),length_ss,i,0);
        diskFrames.block(0,4*i,4,4) = frame;
    }



    if(m_residual_error > 1e-6)
    {
        return false;
    }

    return true;

}

void PiecewiseConstantCurvatureModel::get_res(Eigen::MatrixXd &output, Eigen::MatrixXd input)
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

    output.resize(3*(m_number_disks[0]+m_number_disks[1]),1);

    Eigen::Vector4d F_prev;
    F_prev.setZero();
    Eigen::Vector3d M_prev;
    M_prev.setZero();

    for(int ss_i = m_number_disks[0]+m_number_disks[1]; ss_i > 0; ss_i--)
    {

        //Get transformation matrix between disk i and disk i-1
        Eigen::Matrix4d diskTransform;
        double beta = input(3*ss_i-3);
        double gamma = input(3*ss_i-2);
        double epsi = input(3*ss_i-1);
        double k = std::sqrt(beta*beta+gamma*gamma);
        double phi = std::atan2(gamma,beta);
        double theta = k*length_ss(ss_i-1);
        Eigen::Vector3d p;
        p <<    (1-std::cos(theta))*std::cos(phi)/k,
                (1-std::cos(theta))*std::sin(phi)/k,
                std::sin(theta)/k;
        Eigen::Matrix3d Rz;
        Rz << std::cos(phi), -1*std::sin(phi), 0,
              std::sin(phi), std::cos(phi), 0,
              0, 0, 1;
        Eigen::Matrix3d Ry;
        Ry << std::cos(theta), 0, std::sin(theta),
              0, 1, 0,
              -1*std::sin(theta), 0, std::cos(theta);
        Eigen::Matrix3d Rz2;
        Rz2 << std::cos(epsi-phi), -1*std::sin(epsi-phi), 0,
              std::sin(epsi-phi), std::cos(epsi-phi), 0,
              0, 0, 1;
        diskTransform << Rz*Ry*Rz2, p,
                         0, 0, 0, 1;

        Eigen::Matrix<double,3,6> tendon_pos;
        for(int i = 0; i < 6; i++)
        {
            tendon_pos.col(i) = diskTransform.block(0,0,3,3)*m_routing.at(i) + diskTransform.block(0,3,3,1);
        }

        Eigen::MatrixXd pt_mat;
        pt_mat.resize(4,6);
        pt_mat << m_routing.at(0) - tendon_pos.col(0), m_routing.at(1) - tendon_pos.col(1), m_routing.at(2) - tendon_pos.col(2), m_routing.at(3) - tendon_pos.col(3), m_routing.at(4) - tendon_pos.col(4), m_routing.at(5) - tendon_pos.col(5),
                  0, 0, 0, 0, 0, 0;

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

        Eigen::MatrixXd F_rel;
        F_rel.resize(4,6);
        F_rel.setZero();

        Eigen::MatrixXd M_rel;
        M_rel.resize(3,6);
        M_rel.setZero();

        Eigen::Vector4d zi = diskTransform.col(2);

        if(ss_i < m_number_disks[0]+m_number_disks[1])
        {
            Eigen::Matrix4d tempTransform = getDiskTransform(input,length_ss,ss_i+1,ss_i-1);

            Eigen::Matrix<double,3,6> tendon_pos_temp;
            for(int i = 0; i < 6; i++)
            {
                tendon_pos_temp.col(i) = tempTransform.block(0,0,3,3)*m_routing.at(i) + tempTransform.block(0,3,3,1);
            }

            Eigen::MatrixXd pt1_mat;
            pt1_mat.resize(4,6);
            pt1_mat <<tendon_pos_temp.col(0) - tendon_pos.col(0), tendon_pos_temp.col(1) - tendon_pos.col(1), tendon_pos_temp.col(2) - tendon_pos.col(2), tendon_pos_temp.col(3) - tendon_pos.col(3), tendon_pos_temp.col(4) - tendon_pos.col(4),tendon_pos_temp.col(5) - tendon_pos.col(5),
                      0, 0, 0, 0, 0, 0;


            Eigen::VectorXd norm_ct2;
            norm_ct2.resize(6);
            for(int i = 0; i < 6; i++)
            {
                norm_ct2(i) = (tendon_pos_temp.col(i)- tendon_pos.col(i)).norm();
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

            F_rel = pt_mat.array()/ct1_mat.array()*F_disk_mat.array() + pt1_mat.array()/ct2_mat.array()*F_disk_mat1.array();
            if(ss_i == m_number_disks[0])
            {
                for(int i = 3; i < 6; i++)
                {
                    F_rel.col(i) = F_rel.col(i) - zi.dot(F_rel.col(i))*zi/zi.norm()/zi.norm();
                }

            }
            else
            {
                for(int i = 0; i < 6; i++)
                {
                    F_rel.col(i) = F_rel.col(i) - zi.dot(F_rel.col(i))*zi/zi.norm()/zi.norm();
                }
            }

        }
        else if(ss_i == m_number_disks[0]+m_number_disks[1])
        {
            F_rel = pt_mat.array()/ct1_mat.array()*F_disk_mat.array();
            F_prev.setZero();
            M_prev.setZero();

        }

        for(int i = 0; i < 6; i++)
        {
            Eigen::Vector3d pos = tendon_pos.col(i);
            Eigen::Vector3d F = F_rel.block(0,i,3,1);
            M_rel.col(i) = pos.cross(F);
        }

        Eigen::Vector4d F_ext;
        F_ext.setZero();
        Eigen::Vector3d M_ext;
        M_ext.setZero();
        if(ss_i == m_number_disks[0]+m_number_disks[1])
        {
            Eigen::Matrix4d Rt = getDiskTransform(input,length_ss,ss_i-1,0);

            F_ext = Rt.inverse()*Eigen::Vector4d(m_f_ext(0),m_f_ext(1),m_f_ext(2),0);
            Eigen::Vector3d F_ext_temp = F_ext.block(0,0,3,1);
            M_ext = Rt.block(0,0,3,3).transpose()*m_l_ext + p.cross(F_ext_temp);
        }

        F_prev = diskTransform*F_prev;
        M_prev = diskTransform.block(0,0,3,3)*M_prev;

        Eigen::Vector3d F_ext_temp = F_prev.block(0,0,3,1);

        Eigen::Vector4d F_net = F_rel.rowwise().sum() + F_prev + F_ext;
        F_net(3) = 0;
        Eigen::Vector3d M_net = M_rel.rowwise().sum() + M_prev + M_ext + p.cross(F_ext_temp);


        double I = 0.25*M_PI*(m_ro*m_ro*m_ro*m_ro);
        double G = m_youngs_modulus/(2*(1.3));

        Eigen::Vector3d M_bend = Rz*Eigen::Vector3d(0,k*m_youngs_modulus*I,0);
        Eigen::Vector3d M_tor = diskTransform.block(0,0,3,3)*Eigen::Vector3d(0, 0, 2*I*G*epsi/length_ss(ss_i-1));

        output.block(3*ss_i-3,0,3,1) = M_bend + M_tor - M_net;

        F_prev = F_net;
        M_prev = M_net;

    }



}

void PiecewiseConstantCurvatureModel::setKeepInits(bool keep)
{
    m_keep_inits = keep;
}


