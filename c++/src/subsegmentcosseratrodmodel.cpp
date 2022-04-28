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

#include "subsegmentcosseratrodmodel.h"

// Wrapper function for ODEs for GNU::GSL
int differential_equations_sscr(double t, const double input[2], double deriv[2], void *param)
{

    (void)(t); /* avoid unused parameter warning */

    Eigen::Matrix<double,1,19> y;
    Eigen::MatrixXd dy;

    //Map input to y
    for(int i = 0; i < 19; i++)
    {
        y(0,i) = input[i];
    }


    SubsegmentCosseratRodModel *robot = (SubsegmentCosseratRodModel*)param;
    robot->cosserat_ode(dy,y);

    //Map dy to deriv
    for(int i = 0; i < 19; i++)
    {
        deriv[i] = dy(0,i);
    }

    return GSL_SUCCESS;
}

// Wrapper function for non-linear boundary conditions for GNU::GSL
int evaluate_boundary_conditions_sscr(const gsl_vector * x, void *param, gsl_vector * f)
{


    SubsegmentCosseratRodModel *model = (SubsegmentCosseratRodModel*)param;

    Eigen::MatrixXd inputs;

    inputs.resize(6*model->getNumberOfTotalDisks(),1);


    Eigen::MatrixXd outputs;

    for(int i = 0; i < 6*model->getNumberOfTotalDisks(); i++)
    {
        inputs(i,0) = gsl_vector_get (x, i);
    }


    //Function that needs to be minimized goes here

    model->get_res(outputs, inputs);



    //Functions that needs to be minimized ends here

    for(int i = 0; i < 6*model->getNumberOfTotalDisks(); i++)
    {

        gsl_vector_set (f, i, outputs(i,0));
    }


    return GSL_SUCCESS;
}

// Callback function that can be used during non-linear squares solving
void
callback_sscr(const size_t iter, void *param,
              const gsl_multifit_nlinear_workspace *w)
{
    SubsegmentCosseratRodModel *model = (SubsegmentCosseratRodModel*)param;


    //Pring current iteration
    std::cout << "Iteration " << iter << std::endl;

    gsl_vector * x = gsl_multifit_nlinear_position(w);
    gsl_vector * r = gsl_multifit_nlinear_residual(w);

    std::cout << "Current values:" << std::endl;
    for(int i = 0; i < 6*model->getNumberOfTotalDisks(); i++)
    {
        std::cout << " " << gsl_vector_get(x, i);
    }
    std::cout << std::endl << "Current residuals:" << std::endl;
    for(int i = 0; i < 6*model->getNumberOfTotalDisks(); i++)
    {
        std::cout << " " << gsl_vector_get(r, i);
    }

    double chisq;
    gsl_blas_ddot(r, r, &chisq);

    std::cout << std::endl << "Current cost: "<< chisq << std::endl << std::endl;

}

SubsegmentCosseratRodModel::SubsegmentCosseratRodModel()
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
    m_ode_sys.function = differential_equations_sscr;
    m_ode_sys.jacobian = nullptr;
    m_ode_sys.dimension = 19;
    m_ode_sys.params = this;

    m_ode_step_type = gsl_odeiv2_step_rkf45;

    m_ode_step = gsl_odeiv2_step_alloc (m_ode_step_type, 19);
    m_ode_controller = gsl_odeiv2_control_y_new (1e-6, 1e-3);
    m_ode_evolve = gsl_odeiv2_evolve_alloc (19);

    m_keep_inits = false;
    m_default_inits.resize(1,6*getNumberOfTotalDisks());
    for(int i = 0; i < getNumberOfTotalDisks();  i++)
    {
        m_default_inits.block(0,6*i,1,6) << 0, 0, 1, 0, 0, 0;

    }
    m_last_inits = m_default_inits;

    m_current_segment = 1;

    m_tau.setZero();

    m_length_ss.resize(m_number_disks[0]+m_number_disks[1]);

    for(int i = 0; i < m_number_disks[0]; i++)
    {
        m_length_ss(i) = m_length[0]/m_number_disks[0];
    }

    for(int i = m_number_disks[0]; i < m_number_disks[0]+m_number_disks[1]; i++)
    {
        m_length_ss(i) = m_length[1]/m_number_disks[1];
    }

}

double SubsegmentCosseratRodModel::getNumberOfTotalDisks()
{
    return m_number_disks[0] + m_number_disks[1];
}

void SubsegmentCosseratRodModel::setDefaultInitValues(Eigen::MatrixXd inits)
{
    m_default_inits = inits;
}

Eigen::MatrixXd SubsegmentCosseratRodModel::getFinalInitValues()
{
    return m_last_inits;
}

double SubsegmentCosseratRodModel::getFinalResudial()
{
    return m_residual_error;
}

void SubsegmentCosseratRodModel::setRobotParameters(std::array<double, 2> length, double youngs_modulus, std::vector<Eigen::Vector3d> routing, std::array<int, 2> number_disks, std::array<double, 2> pradius_disks, double ro)
{
    m_length[0]         = length[0];
    m_length[1]         = length[1];
    m_youngs_modulus    = youngs_modulus;
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


    m_length_ss.resize(m_number_disks[0]+m_number_disks[1]);

    for(int i = 0; i < m_number_disks[0]; i++)
    {
        m_length_ss(i) = m_length[0]/m_number_disks[0];
    }

    for(int i = m_number_disks[0]; i < m_number_disks[0]+m_number_disks[1]; i++)
    {
        m_length_ss(i) = m_length[1]/m_number_disks[1];
    }
	
    m_default_inits.resize(1,6*getNumberOfTotalDisks());
    for(int i = 0; i < getNumberOfTotalDisks();  i++)
    {
        m_default_inits.block(0,6*i,1,6) << 0, 0, 1, 0, 0, 0;

    }
    m_last_inits = m_default_inits;
}

SubsegmentCosseratRodModel::~SubsegmentCosseratRodModel()
{
    gsl_odeiv2_evolve_free (m_ode_evolve);
    gsl_odeiv2_control_free (m_ode_controller);
    gsl_odeiv2_step_free (m_ode_step);
}

void SubsegmentCosseratRodModel::run_IVP(Eigen::MatrixXd &states, Eigen::MatrixXd init_state)
{

    double t = 0;
    double t1 = m_length_ss[m_current_segment-1];
    double y[19];

    //Variable step length
    double h = m_length_ss[m_current_segment-1];
    //Fixed step length
    //double h = m_length_ss[m_current_segment-1];





    for(int i = 0; i < 19; i++)
    {
        y[i] = init_state(0,i);
    }


    //Reset ODE stuff from previous run

    std::vector<Eigen::MatrixXd> y_stacked;
    int status;
    while (t < t1)
    {
        //Variable step length
        status = gsl_odeiv2_evolve_apply (m_ode_evolve, m_ode_controller, m_ode_step, &m_ode_sys, &t, t1, &h, y);
        //Fixed step length
        //int status = gsl_odeiv2_evolve_apply_fixed_step (m_ode_evolve, m_ode_controller, m_ode_step, &m_ode_sys, &t, h, y);

        if (status != GSL_SUCCESS)
            break;
    }
    if(!(status != GSL_SUCCESS))
    {
        Eigen::Matrix<double,19,1> temp;
        for(int k = 0; k < 19; k++)
        {
            temp(k) = y[k];
        }

        if(t <= t1) //Only add the last step if it's still within the given boundaries
            y_stacked.push_back(temp);
    }

    states.resize(y_stacked.size(),19);
    states.setZero();
    for(int i = 0; i < y_stacked.size(); i++)
    {
        states.row(i) = y_stacked.at(i).transpose();
    }

    //Reset stuff for next run
    gsl_odeiv2_step_reset(m_ode_step);
    gsl_odeiv2_evolve_reset(m_ode_evolve);

}

bool SubsegmentCosseratRodModel::forwardKinematics(Eigen::MatrixXd &diskFrames, Eigen::MatrixXd q, Eigen::Vector3d f_ext, Eigen::Vector3d l_ext)
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

    const size_t n = 6*getNumberOfTotalDisks(); //Number of equations (boundary conditions)
    const size_t p = 6*getNumberOfTotalDisks(); //Number of parameters (initial guesses)
    gsl_vector *x0 = gsl_vector_alloc(p);
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    fdf.f = evaluate_boundary_conditions_sscr;
    fdf.df = NULL; //Provide your own Jacobian
    fdf.fvv = NULL; //Provide your own Hessian
    fdf.n = n;
    fdf.p = p;
    fdf.params = this;

    //Set initial values
    for(int i = 0; i < 6*getNumberOfTotalDisks(); i++)
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
    //gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,callback_sscr, this, &info, work);
    //Without Callback
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,NULL, NULL, &info, work);

    /* store final cost */
    gsl_blas_ddot(f, f, &chisq);
    m_residual_error = chisq;

    /* store cond(J(x)) */
    gsl_multifit_nlinear_rcond(&rcond, work);

    ///---End Solve the system---

    Eigen::MatrixXd final_values;
    final_values.resize(1,6*getNumberOfTotalDisks());

    Eigen::MatrixXd final_boundcon;
    final_boundcon.resize(1,6*getNumberOfTotalDisks());


    for(int i = 0; i < 6*getNumberOfTotalDisks(); i ++)
    {
        final_values(0,i) = gsl_vector_get (x, i);
        final_boundcon(0,i) = gsl_vector_get (f, i);
    }

    m_last_inits = final_values;
    gsl_vector_free (x0);
    gsl_multifit_nlinear_free(work);

    //Get the last (initial) states
    Eigen::MatrixXd states = m_states;

    Eigen::MatrixXd tempDisks;

    tempDisks.resize(states.rows(),12);
    tempDisks << states.block(0,0,states.rows(),12);


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


void SubsegmentCosseratRodModel::cosserat_ode(Eigen::MatrixXd &dy, Eigen::MatrixXd y)
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

    Eigen::Vector3d le;
    le.setZero();

    Eigen::Vector3d fe;
    fe.setZero();




    Eigen::Vector3d v_dot;
    Eigen::Vector3d u_dot;

    v_dot = -1*m_Kse.inverse()*(u_hat*m_Kse*(v - Eigen::Vector3d(0,0,1))+R.transpose()*fe);
    u_dot = -1*m_Kbt.inverse()*(u_hat*m_Kbt*u + v_hat*m_Kse*(v - Eigen::Vector3d(0,0,1))+R.transpose()*le);

    Eigen::Vector3d p_dot;
    p_dot = R*v;

    Eigen::Matrix3d R_dot;
    R_dot = R*u_hat;

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> R_row(R_dot);
    Eigen::Map<Eigen::RowVectorXd> R_flat(R_row.data(), R_row.size());


    dy << p_dot.transpose(),R_flat,v_dot.transpose(),u_dot.transpose(),1;




}

void SubsegmentCosseratRodModel::get_res(Eigen::MatrixXd &output, Eigen::MatrixXd input)
{

    output.resize(6*getNumberOfTotalDisks(),1);
    output.setZero();

    //Define variables
    Eigen::MatrixXd states_current;
    Eigen::MatrixXd states_next;
    Eigen::MatrixXd states_prev;
    m_states.resize(m_number_disks[0]+m_number_disks[1]+1,19);

    //Iterate through each segment and run IVP
    Eigen::Matrix<double,1,19> init_state_segment; //p, R, v, u, s
    init_state_segment << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, input.block(0,0,6,1).transpose(), 0;
    m_states.row(0) = init_state_segment;


    for(int i = 0; i < m_number_disks[0]+m_number_disks[1];i++)
    {
        //Run IVP for current segment
        m_current_segment = i+1;
        run_IVP(states_current,init_state_segment);

        //Save current state in state vector
        m_states.row(i+1) = states_current;


        //Evaluate force and moment equlibiurm at the end of current segment based on the interal forces and moments (n & m)
        Eigen::Matrix3d R_c; //current R and p
        R_c << states_current.block(0,3,1,3),
                states_current.block(0,6,1,3),
                states_current.block(0,9,1,3);

        Eigen::Vector3d u_c = states_current.block(0,15,1,3).transpose(); //current u and v
        Eigen::Vector3d v_c = states_current.block(0,12,1,3).transpose();
        Eigen::Vector3d u_n = Eigen::Vector3d::Zero(); //next u and v
        Eigen::Vector3d v_n;
        v_n << 0,
               0,
               1;
        if(i != m_number_disks[0]+m_number_disks[1] - 1) //If not at last segment
        {
            u_n = input.block(6*(i+1)+3,0,3,1); //next u and v
            v_n = input.block(6*(i+1),0,3,1);
        }

        Eigen::Vector3d n_c = R_c*m_Kse*(v_c - Eigen::Vector3d(0,0,1)); //current n and m
        Eigen::Vector3d m_c = R_c*m_Kbt*u_c;
        Eigen::Vector3d n_n = R_c*m_Kse*(v_n - Eigen::Vector3d(0,0,1)); //next n and m
        Eigen::Vector3d m_n = R_c*m_Kbt*u_n;

        output.block(6*i,0,3,1) = n_c - n_n; //Update force equilibiurm to include inner forces
        output.block(6*i+3,0,3,1) = m_c - m_n;  //Update moment equilibiurm to include inner moments



        //Update init state for next segment
        if(i != m_number_disks[0]+m_number_disks[1] - 1) //If not at last segment
        {
            init_state_segment = states_current; //Copy state at the end of current segment
            init_state_segment.block(0,12,1,6) = input.block(6*(i+1),0,6,1).transpose();//Set v and u to the guessed values
        }
    }


    //Iterate through all states again to include tendon and external forces in the BC (we need to do this afterwards, so that the shape and the tendon locations are known)
    for(int i = 1; i < m_number_disks[0]+m_number_disks[1]+1;i++)
    {
        //Get current, previous and next disk poses and tendon positions
        Eigen::Matrix4d T_cur;
        Eigen::Matrix<double,3,6> tendon_pos_cur;
        T_cur << m_states.row(i).block(0,3,1,3), m_states.row(i).block(0,0,1,1),
                m_states.row(i).block(0,6,1,3), m_states.row(i).block(0,1,1,1),
                m_states.row(i).block(0,9,1,3), m_states.row(i).block(0,2,1,1),
                0, 0, 0, 1;
        for(int k = 0; k < 6; k++)
        {
            tendon_pos_cur.col(k) = T_cur.block(0,0,3,3)*m_routing.at(k) + T_cur.block(0,3,3,1);
        }

        Eigen::Matrix4d T_prev;
        Eigen::Matrix<double,3,6> tendon_pos_prev;
        T_prev <<  m_states.row(i-1).block(0,3,1,3),  m_states.row(i-1).block(0,0,1,1),
                m_states.row(i-1).block(0,6,1,3),  m_states.row(i-1).block(0,1,1,1),
                m_states.row(i-1).block(0,9,1,3),  m_states.row(i-1).block(0,2,1,1),
                0, 0, 0, 1;
        for(int k = 0; k < 6; k++)
        {
            tendon_pos_prev.col(k) = T_prev.block(0,0,3,3)*m_routing.at(k) + T_prev.block(0,3,3,1);
        }

        Eigen::Matrix4d T_next;
        Eigen::Matrix<double,3,6> tendon_pos_next;

        if(i != m_number_disks[0]+m_number_disks[1]) //Next disk does not exist for last disk
        {
            T_next <<  m_states.row(i+1).block(0,3,1,3),  m_states.row(i+1).block(0,0,1,1),
                    m_states.row(i+1).block(0,6,1,3),  m_states.row(i+1).block(0,1,1,1),
                    m_states.row(i+1).block(0,9,1,3),  m_states.row(i+1).block(0,2,1,1),
                    0, 0, 0, 1;
            for(int k = 0; k < 6; k++)
            {
                tendon_pos_next.col(k) = T_next.block(0,0,3,3)*m_routing.at(k) + T_next.block(0,3,3,1);
            }
        }


        //Collect tendon forces acting on each disk in a matrix
        Eigen::MatrixXd F_disk;
        F_disk.resize(m_number_disks[0]+m_number_disks[1],6);
        F_disk.setZero();

        for(int k = 0; k < m_number_disks[0]; k++)
        {
            F_disk.block(k,0,1,6) = m_tau.transpose();
        }

        for(int k = m_number_disks[0]; k < m_number_disks[0] + m_number_disks[1]; k++)
        {
            F_disk.block(k,3,1,3) = m_tau.block(3,0,3,1).transpose();
        }

        //Define tendon forces and moments and normal vector of disk
        Eigen::Vector3d F_tendon(0,0,0);
        Eigen::Vector3d zi = T_cur.block(0,2,3,1);

        Eigen::Vector3d M_tendon(0,0,0);



        //Update force equilibrium of the disk w.r.t. tendon forces (and external forces if last segment)

        if(i == m_number_disks[0]+m_number_disks[1]) // if end of second segment
        {
            for(int k = 0; k < 6; k++) //There are no tendons routed to a next disk, cause there are none
            {

                Eigen::Vector3d force = (tendon_pos_prev.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i-1,k);
                Eigen::Vector3d pos = tendon_pos_cur.col(k) - T_cur.block(0,3,3,1);
                F_tendon += force;
                M_tendon += pos.cross(force);
            }

            //Update force and moment equilibrium with tendon and external forces (last disk)
            output.block(6*(i-1),0,3,1) += -F_tendon - m_f_ext; //Update force equilibiurm to include tendon and external forces
            output.block(6*(i-1)+3,0,3,1) += -M_tendon - m_l_ext;  //Update moment equilibiurm to include tendon and external moments

        }
        else if(i == m_number_disks[0]) // if end of first segment
        {
            for(int k = 0; k < 6; k++) // Tendon forces from previous and next disk
            {

                Eigen::Vector3d force1 = (tendon_pos_prev.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i-1,k);
                Eigen::Vector3d force2 = (tendon_pos_next.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i,k);
                Eigen::Vector3d pos = tendon_pos_cur.col(k) - T_cur.block(0,3,3,1);

                if(k > 2) // Eliminate the comoponent of the tendon force (only last three, cause first three are attached to this disk) that is normal to the disk (we do not include friction)
                {
                    force1 -= zi.dot((tendon_pos_prev.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i-1,k))*zi;
                    force2 -= zi.dot((tendon_pos_next.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i,k))*zi;
                }

                F_tendon += force1;
                F_tendon += force2;
                M_tendon += pos.cross(force1);
                M_tendon += pos.cross(force2);

            }
            output.block(6*(i-1),0,3,1) += -F_tendon; //Update force equilibiurm to include tendon and external forces
            output.block(6*(i-1)+3,0,3,1) += -M_tendon;  //Update moment equilibiurm to include tendon and external moments

        }
        else // all other segments
        {
            for(int k = 0; k < 6; k++) // Tendon forces from previous and next disk
            {
                Eigen::Vector3d force1 = (tendon_pos_prev.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i-1,k);
                Eigen::Vector3d force2 = (tendon_pos_next.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i,k);
                Eigen::Vector3d pos = tendon_pos_cur.col(k) - T_cur.block(0,3,3,1);

                // Eliminate the comoponent of the tendon force that is normal to the disk (we do not include friction)
                force1 -= zi.dot((tendon_pos_prev.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i-1,k))*zi;
                force2 -= zi.dot((tendon_pos_next.col(k)-tendon_pos_cur.col(k)).normalized()*F_disk(i,k))*zi;

                F_tendon += force1;
                F_tendon += force2;
                M_tendon += pos.cross(force1);
                M_tendon += pos.cross(force2);
            }
            output.block(6*(i-1),0,3,1) += -F_tendon; //Update force equilibiurm to include tendon and external forces
            output.block(6*(i-1)+3,0,3,1) += -M_tendon;  //Update moment equilibiurm to include tendon and external moments

        }




    }

}

void SubsegmentCosseratRodModel::setKeepInits(bool keep)
{
    m_keep_inits = keep;
}


