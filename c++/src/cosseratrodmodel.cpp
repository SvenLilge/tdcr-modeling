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

#include "cosseratrodmodel.h"

// Wrapper function for ODEs for GNU::GSL
int differential_equations_cr(double t, const double input[2], double deriv[2], void *param)
{

    (void)(t); /* avoid unused parameter warning */

    Eigen::Matrix<double,1,19> y;
    Eigen::MatrixXd dy;

    //Map input to y
    for(int i = 0; i < 19; i++)
    {
        y(0,i) = input[i];
    }


    CosseratRodModel *robot = (CosseratRodModel*)param;
    robot->cosserat_ode(dy,y);

    //Map dy to deriv
    for(int i = 0; i < 19; i++)
    {
        deriv[i] = dy(0,i);
    }

    return GSL_SUCCESS;
}

// Wrapper function for non-linear boundary conditions for GNU::GSL
int evaluate_boundary_conditions_cr(const gsl_vector * x, void *param, gsl_vector * f)
{

    Eigen::Matrix<double,6,1> inputs;

    Eigen::MatrixXd outputs;

    for(int i = 0; i < 6; i++)
    {
        inputs(i,0) = gsl_vector_get (x, i);
    }

    CosseratRodModel *model = (CosseratRodModel*)param;
    model->get_res(outputs, inputs);

    for(int i = 0; i < 6; i++)
    {

        gsl_vector_set (f, i, outputs(i,0));
    }


    return GSL_SUCCESS;
}


// Callback function that can be used during non-linear squares solving
void
callback_cr(const size_t iter, void *params,
            const gsl_multifit_nlinear_workspace *w)
{

    //Pring current iteration
    std::cout << "Iteration " << iter << std::endl;

    gsl_vector * x = gsl_multifit_nlinear_position(w);
    gsl_vector * r = gsl_multifit_nlinear_residual(w);

    std::cout << "Current values:" << std::endl;
    for(int i = 0; i < 6; i++)
    {
        std::cout << " " << gsl_vector_get(x, i);
    }
    std::cout << std::endl << "Current residuals:" << std::endl;
    for(int i = 0; i < 6; i++)
    {
        std::cout << " " << gsl_vector_get(r, i);
    }

    double chisq;
    gsl_blas_ddot(r, r, &chisq);

    std::cout << std::endl << "Current cost: "<< chisq << std::endl << std::endl;

}

CosseratRodModel::CosseratRodModel()
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

    //Set stiffness matrices
    double G = m_youngs_modulus/(2*(1.3));
    double I = 0.25*M_PI*(m_ro*m_ro*m_ro*m_ro);
    double A = M_PI*(m_ro*m_ro);

    m_Kse << G*A, 0, 0,
            0, G*A, 0,
            0, 0, m_youngs_modulus*A;

    m_Kbt << m_youngs_modulus*I, 0, 0,
            0, m_youngs_modulus*I, 0,
            0, 0, 2*G*I;



    //Setup ODE stuff
    m_ode_sys.function = differential_equations_cr;
    m_ode_sys.jacobian = nullptr;
    m_ode_sys.dimension = 19;
    m_ode_sys.params = this;

    m_ode_step_type = gsl_odeiv2_step_rkf45;

    m_ode_step = gsl_odeiv2_step_alloc (m_ode_step_type, 19);
    m_ode_controller = gsl_odeiv2_control_y_new (1e-6, 1e-3);
    m_ode_evolve = gsl_odeiv2_evolve_alloc (19);

    m_keep_inits = false;
    m_default_inits.resize(1,6);
    m_default_inits << 0, 0, 1, 0, 0, 0;
    m_last_inits = m_default_inits;

    m_current_segment = 1;

    m_tau.setZero();

}

void CosseratRodModel::setRobotParameters(std::array<double, 2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int, 2> number_disks, std::array<double, 2> pradius_disks, double ro)
{
	//Update parameters
    m_length[0]         = length[0];
    m_length[1]         = length[1];
    m_youngs_modulus    = youngs_modulus;
    m_number_disks[0]   = number_disks[0];
    m_number_disks[1]   = number_disks[1];
    m_pradius_disks[0]  = pradius_disks[0];
    m_pradius_disks[1]  = pradius_disks[1];
    m_ro                = ro;
    m_routing           = routing;

    //Set stiffness matrices
    double G = m_youngs_modulus/(2*(1.3));
    double I = 0.25*M_PI*(m_ro*m_ro*m_ro*m_ro);
    double A = M_PI*(m_ro*m_ro);

    m_Kse << G*A, 0, 0,
            0, G*A, 0,
            0, 0, m_youngs_modulus*A;

    m_Kbt << m_youngs_modulus*I, 0, 0,
            0, m_youngs_modulus*I, 0,
            0, 0, 2*G*I;

}

void CosseratRodModel::setDefaultInitValues(Eigen::MatrixXd inits)
{
    m_default_inits = inits;
}

Eigen::MatrixXd CosseratRodModel::getFinalInitValues()
{
    return m_last_inits;
}

double CosseratRodModel::getFinalResudial()
{
    return m_residual_error;
}

CosseratRodModel::~CosseratRodModel()
{
    gsl_odeiv2_evolve_free (m_ode_evolve);
    gsl_odeiv2_control_free (m_ode_controller);
    gsl_odeiv2_step_free (m_ode_step);
}

void CosseratRodModel::run_IVP(Eigen::MatrixXd &states, Eigen::MatrixXd init_state)
{

    double t = 0;
    double t1 = m_length[m_current_segment-1];
    double y[19];

    //Variable step length
    //double h = 0.005;
    //Fixed step length
    double h = t1/m_number_disks[m_current_segment-1];
    double h_init = t1/m_number_disks[m_current_segment-1];





    for(int i = 0; i < 19; i++)
    {
        y[i] = init_state(0,i);
    }


    //Reset ODE stuff from previous run

    std::vector<Eigen::MatrixXd> y_stacked;
    for(int i = 0; i < m_number_disks[m_current_segment-1]; i++)
    {
        h = h_init;
        t1 = t+h_init;
        int status;
        while (t < t1)
        {
            //Variable step length
            status = gsl_odeiv2_evolve_apply (m_ode_evolve, m_ode_controller, m_ode_step, &m_ode_sys, &t, t1, &h, y);
            //Fixed step length
            //status = gsl_odeiv2_evolve_apply_fixed_step (m_ode_evolve, m_ode_controller, m_ode_step, &m_ode_sys, &t, h, y);

            if (status != GSL_SUCCESS)
                break;



        }
        if (!(status != GSL_SUCCESS))
        {
            Eigen::Matrix<double,19,1> temp;
            for(int k = 0; k < 19; k++)
            {
                temp(k) = y[k];
            }

            if(t <= t1) //Only add the last step if it's still within the given boundaries
                y_stacked.push_back(temp);
        }
    }


    if(m_current_segment == 1) // If first segment add the initial state as first disk (disk that is fixed at the base)
    {
        states.resize(y_stacked.size()+1,19);
        states.setZero();
        states.row(0) = init_state;
        for(int i = 0; i < y_stacked.size(); i++)
        {
            states.row(i+1) = y_stacked.at(i).transpose();
        }
    }
    else // Not needed for second segment, last disk of first segment is the same as first for second segment
    {
        states.resize(y_stacked.size(),19);
        states.setZero();
        for(int i = 0; i < y_stacked.size(); i++)
        {
            states.row(i) = y_stacked.at(i).transpose();
        }
    }

    //Reset stuff for next run
    gsl_odeiv2_step_reset(m_ode_step);
    gsl_odeiv2_evolve_reset(m_ode_evolve);

}

bool CosseratRodModel::forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext)
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

    const size_t n = 6; //Number of equations (boundary conditions)
    const size_t p = 6; //Number of parameters (initial guesses)
    gsl_vector *x0 = gsl_vector_alloc(p);
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    fdf.f = evaluate_boundary_conditions_cr;
    fdf.df = NULL; //Provide your own Jacobian
    fdf.fvv = NULL; //Provide your own Hessian
    fdf.n = n;
    fdf.p = p;
    fdf.params = this;

    //Set initial values
    for(int i = 0; i < 6; i++)
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
    //gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,callback_cr, NULL, &info, work);
    //Without Callback
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,NULL, NULL, &info, work);

    /* store final cost */
    gsl_blas_ddot(f, f, &chisq);
    m_residual_error = chisq;

    /* store cond(J(x)) */
    gsl_multifit_nlinear_rcond(&rcond, work);

    ///---End Solve the system---

    Eigen::Matrix<double,1,6> final_values;
    Eigen::Matrix<double,1,6> final_boundcon;


    for(int i = 0; i < 6; i ++)
    {
        final_values(0,i) = gsl_vector_get (x, i);
        final_boundcon(0,i) = gsl_vector_get (f, i);
    }

    m_last_inits = final_values;
    gsl_vector_free (x0);
    gsl_multifit_nlinear_free(work);

    //Get the last (initial) states
    Eigen::MatrixXd states1 = m_states1;
    Eigen::MatrixXd states2 = m_states2;

    Eigen::MatrixXd tempDisks;

    if(states1.rows() == 0 || states1.rows() == 1)
        return false;

    tempDisks.resize(states1.rows()+states2.rows(),12);
    tempDisks << states1.block(0,0,states1.rows(),12),
            states2.block(0,0,states2.rows(),12);


    diskFrames.resize(4,4*tempDisks.rows());

    for(int i = 0; i < tempDisks.rows(); i++)
    {
        Eigen::Matrix3d R;
        Eigen::Vector3d p;
        R << tempDisks.row(i).block(0,3,1,3),
                tempDisks.row(i).block(0,6,1,3),
                tempDisks.row(i).block(0,9,1,3);
        p = tempDisks.row(i).block(0,0,1,3).transpose();
        Eigen::Matrix4d frame;
        frame << R, p,
                0,0,0,1;
        diskFrames.block(0,4*i,4,4) = frame;
    }



    if(m_residual_error > 1e-6)
    {
        return false;
    }

    return true;

}


void CosseratRodModel::cosserat_ode(Eigen::MatrixXd &dy, Eigen::MatrixXd y)
{
    dy.resize(1,19);
    dy.setZero();


    Eigen::Vector3d u = y.block(0,15,1,3).transpose();
    Eigen::Vector3d v = y.block(0,12,1,3).transpose();


    Eigen::Matrix3d u_hat;
    u_hat << 0, -u(2), u(1),
            u(2), 0, -u(0),
            -u(1), u(0), 0;

    Eigen::Matrix3d v_hat;
    v_hat << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0), 0;

    Eigen::Matrix3d R;
    R << y.block(0,3,1,3),
            y.block(0,6,1,3),
            y.block(0,9,1,3);


    std::vector<Eigen::Matrix3d> r_hat;
    for(int i = 0; i < 6; i++) // Calculate stuff for each tendon
    {
        Eigen::Matrix3d temp;
        temp << 0, -m_routing.at(i)(2), m_routing.at(i)(1),
                m_routing.at(i)(2), 0, -m_routing.at(i)(0),
                -m_routing.at(i)(1), m_routing.at(i)(0), 0;

        r_hat.push_back(temp);
    }

    Eigen::MatrixXd pb_dot;
    pb_dot.resize(3,6);

    std::vector<Eigen::Matrix3d> A;

    Eigen::Matrix3d A_total;
    A_total.setZero();

    for(int i = 0; i < 6; i++)
    {
        pb_dot.block(0,i,3,1) = u_hat*m_routing.at(i) + v;

        Eigen::Matrix3d pb_dot_hat;
        pb_dot_hat <<0, -pb_dot(2,i), pb_dot(1,i),
                pb_dot(2,i), 0, -pb_dot(0,i),
                -pb_dot(1,i), pb_dot(0,i), 0;

        A.push_back(-(m_tau(i)/std::pow(pb_dot.block(0,i,3,1).norm(),3))*pb_dot_hat*pb_dot_hat);

        if(m_current_segment == 1 || i > 2) // Take into account all tendons for first segment, but only last four for second segment
            A_total += A.at(i);
    }

    std::vector<Eigen::Matrix3d> B;

    Eigen::Matrix3d B_total;
    B_total.setZero();

    for(int i = 0; i < 6; i++)
    {
        B.push_back(r_hat.at(i)*A.at(i));
        if(m_current_segment == 1 || i > 2) // Take into account all tendons for first segment, but only last four for second segment
            B_total += B.at(i);
    }

    Eigen::Matrix3d G;
    G.setZero();
    for(int i = 0; i < 6; i++)
    {
        if(m_current_segment == 1 || i > 2) // Take into account all tendons for first segment, but only last four for second segment
            G += -A.at(i)*r_hat.at(i);
    }

    Eigen::Matrix3d H;
    H.setZero();
    for(int i = 0; i < 6; i++)
    {
        if(m_current_segment == 1 || i > 2) // Take into account all tendons for first segment, but only last four for second segment
            H += -B.at(i)*r_hat.at(i);
    }

    std::vector<Eigen::Vector3d> a;
    Eigen::Vector3d a_total;
    a_total.setZero();
    for(int i = 0; i < 6; i++)
    {
        a.push_back(A.at(i)*u_hat*pb_dot.block(0,i,3,1));
        if(m_current_segment == 1 || i > 2) // Take into account all tendons for first segment, but only last four for second segment
            a_total += a.at(i);
    }

    std::vector<Eigen::Vector3d> b;
    Eigen::Vector3d b_total;
    b_total.setZero();
    for(int i = 0; i < 6; i++)
    {
        b.push_back(r_hat.at(i)*a.at(i));
        if(m_current_segment == 1 || i > 2) // Take into account all tendons for first segment, but only last four for second segment
            b_total += b.at(i);
    }

    Eigen::Vector3d le;
    le.setZero();

    Eigen::Vector3d fe;
    fe.setZero();

    Eigen::Vector3d c;
    c = -u_hat*m_Kbt*u-v_hat*m_Kse*(v - Eigen::Vector3d(0,0,1))-R.transpose()*le-b_total;

    Eigen::Vector3d d;
    d = -u_hat*m_Kse*(v - Eigen::Vector3d(0,0,1))-R.transpose()*fe-a_total;


    Eigen::VectorXd vu_dot;
    vu_dot.resize(6);
    Eigen::VectorXd dc;
    dc.resize(6);
    dc << d,
            c;

    Eigen::MatrixXd vu_mat;
    vu_mat.resize(6,6);

    vu_mat << m_Kse + A_total, G,
            B_total, m_Kbt+ H;

    vu_dot = vu_mat.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(dc);

    Eigen::Vector3d p_dot;
    p_dot = R*v;

    Eigen::Matrix3d R_dot;
    R_dot = R*u_hat;

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> R_row(R_dot);
    Eigen::Map<Eigen::RowVectorXd> R_flat(R_row.data(), R_row.size());


    dy << p_dot.transpose(),R_flat,vu_dot.transpose(),1;




}

void CosseratRodModel::get_res(Eigen::MatrixXd &output, Eigen::MatrixXd input)
{
    //Define variables
    Eigen::MatrixXd states1;
    Eigen::MatrixXd states2;

    //Run IVP for first segment
    Eigen::Matrix<double,1,19> init_state_segment1; //p, R, v, u, s
    init_state_segment1 << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, input.transpose(), 0;
    m_current_segment = 1; // Set current segment to first for integration along backbone
    run_IVP(states1,init_state_segment1);
    m_states1 = states1;

    //RUN IVP for second segment with initial values according to transition conditions between segments
    Eigen::Matrix<double,1,19> init_state_segment2; //p, R, v, u, s
    Eigen::Vector3d v_new;
    Eigen::Vector3d u_new;
    //calculate initial v and u for next segment by evaluating states at the end of the first segment and respecting transition/boundary conditions

    Eigen::Matrix3d R_end1;
    R_end1 << states1.bottomRows(1).block(0,3,1,3),
            states1.bottomRows(1).block(0,6,1,3),
            states1.bottomRows(1).block(0,9,1,3);

    Eigen::Vector3d u_end1 = states1.bottomRows(1).block(0,15,1,3).transpose();
    Eigen::Vector3d v_end1 = states1.bottomRows(1).block(0,12,1,3).transpose();

    Eigen::Vector3d p_dot1 = R_end1*(u_end1.cross(m_routing.at(0)) + v_end1);
    Eigen::Vector3d p_dot2 = R_end1*(u_end1.cross(m_routing.at(1)) + v_end1);
    Eigen::Vector3d p_dot3 = R_end1*(u_end1.cross(m_routing.at(2)) + v_end1);

    Eigen::Vector3d F_end1 = -m_tau(0)/p_dot1.norm()*p_dot1 -m_tau(1)/p_dot2.norm()*p_dot2 -m_tau(2)/p_dot3.norm()*p_dot3;
    Eigen::Vector3d L_end1 = -m_tau(0)/p_dot1.norm()*(R_end1*m_routing.at(0)).cross(p_dot1) -m_tau(1)/p_dot2.norm()*(R_end1*m_routing.at(1)).cross(p_dot2) -m_tau(2)/p_dot3.norm()*(R_end1*m_routing.at(2)).cross(p_dot3);

    v_new = v_end1 - m_Kse.inverse()*R_end1.transpose()*F_end1;
    u_new = u_end1 - m_Kbt.inverse()*R_end1.transpose()*L_end1;


    init_state_segment2 << states1.bottomRows(1).block(0,0,1,12), v_new.transpose(), u_new.transpose(), states1.bottomRows(1).block(0,18,1,1);
    m_current_segment = 2; // Set current segment to second for integration along backbone
    run_IVP(states2,init_state_segment2);
    m_states2 = states2;

    //Get states at the distal end of the robot
    Eigen::Vector3d n_end2;
    Eigen::Vector3d m_end2;
    Eigen::Vector3d u_end2;
    Eigen::Vector3d v_end2;
    Eigen::Vector3d F_end2;
    Eigen::Vector3d L_end2;
    Eigen::Matrix3d R_end2;

    v_end2 = states2.bottomRows(1).block(0,12,1,3).transpose();
    u_end2 = states2.bottomRows(1).block(0,15,1,3).transpose();
    R_end2 << states2.bottomRows(1).block(0,3,1,3),
            states2.bottomRows(1).block(0,6,1,3),
            states2.bottomRows(1).block(0,9,1,3);

    Eigen::Vector3d p_dot1_2 = R_end2*(u_end2.cross(m_routing.at(3)) + v_end2);
    Eigen::Vector3d p_dot2_2 = R_end2*(u_end2.cross(m_routing.at(4)) + v_end2);
    Eigen::Vector3d p_dot3_2 = R_end2*(u_end2.cross(m_routing.at(5)) + v_end2);


    F_end2 = -m_tau(3)/p_dot1_2.norm()*p_dot1_2 -m_tau(4)/p_dot2_2.norm()*p_dot2_2 -m_tau(5)/p_dot3_2.norm()*p_dot3_2;
    L_end2 = -m_tau(3)/p_dot1_2.norm()*(R_end2*m_routing.at(3)).cross(p_dot1_2) -m_tau(4)/p_dot2_2.norm()*(R_end2*m_routing.at(4)).cross(p_dot2_2) -m_tau(5)/p_dot3_2.norm()*(R_end2*m_routing.at(5)).cross(p_dot3_2);

    n_end2 = R_end2*m_Kse*(v_end2 - Eigen::Vector3d(0,0,1));
    m_end2 = R_end2*m_Kbt*u_end2;

    // 1) Force eqilibrium
    Eigen::Vector3d force_eq = n_end2 - F_end2 - m_f_ext;

    // 2) Moment eqilibrium
    Eigen::Vector3d moment_eq = m_end2-L_end2 - m_l_ext;

    // Set output to resudial error
    output.resize(6,1);
    output << force_eq,
            moment_eq;


}

void CosseratRodModel::setKeepInits(bool keep)
{
    m_keep_inits = keep;
}


