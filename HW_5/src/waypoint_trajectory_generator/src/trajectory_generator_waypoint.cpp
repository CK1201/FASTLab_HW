#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <math.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    MatrixXd M = MatrixXd::Zero(p_num1d * m, p_num1d * m);
    MatrixXd Mi;

    for (int k = 0; k < m; k++)
    {
        Mi = MatrixXd::Zero(p_num1d, p_num1d);
        //time=0
        Mi(0, 0) = 1;
        Mi(1, 1) = 1;
        Mi(2, 2) = 2;
        Mi(3, 3) = 6;
        //p
        for (int i = 0; i < p_num1d;i++){
            Mi(4, i) = pow(Time(k), i);
        }
        //v
        for (int i = 1; i < p_num1d;i++){
            Mi(5, i) = i * pow(Time(k), i - 1);
        }
        //a
        for (int i = 2; i < p_num1d;i++){
            Mi(6, i) = i * (i - 1) * pow(Time(k), i - 2);
        }
        //j
        for (int i = 3; i < p_num1d;i++){
            Mi(7, i) = i * (i - 1) * (i - 2) * pow(Time(k), i - 3);
        }
        // M.block<p_num1d, p_num1d>(k * p_num1d, k * p_num1d) = Mi;
        M.block(k * p_num1d, k * p_num1d, p_num1d, p_num1d) = Mi;
    }
    // cout << "M:" << endl
    //      << M << endl;

    /*   Produce the dereivatives in X, Y and Z axis directly.  */
    MatrixXd Ct = MatrixXd::Zero(2 * d_order * m, d_order * (m + 1));
    Ct(0, 0) = 1;
    Ct(1, 1) = 1;
    Ct(2, 2) = 1;
    Ct(3, 3) = 1;

    Ct(2 * d_order * m - 4, d_order + m - 1) = 1;
    Ct(2 * d_order * m - 3, d_order + m + 0) = 1;
    Ct(2 * d_order * m - 2, d_order + m + 1) = 1;
    Ct(2 * d_order * m - 1, d_order + m + 2) = 1;

    for (int i = 0; i < m-1; i++){
        
        Ct(d_order+2*i*d_order,d_order+i) = 1;
        Ct(d_order+2*i*d_order+d_order,d_order+i) = 1;

        Ct(d_order + 2 * i * d_order + 1, 2 * d_order + m + i * (d_order - 1) - 1) = 1; // v_end
        Ct(d_order+2*i*d_order+2,2*d_order+m+i*(d_order-1)+0)=1;// a_end
        Ct(d_order+2*i*d_order+3,2*d_order+m+i*(d_order-1)+1)=1;// j_end

        Ct(d_order+2*i*d_order+1+d_order,2*d_order+m+i*(d_order-1)-1)=1;// v_start
        Ct(d_order+2*i*d_order+2+d_order,2*d_order+m+i*(d_order-1)+0)=1;// a_start
        Ct(d_order+2*i*d_order+3+d_order,2*d_order+m+i*(d_order-1)+1)=1;// j_start
    }
    // cout << "Ct:" << endl
    //      << Ct << endl;
    auto C = Ct.transpose();

    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */
    MatrixXd Q = MatrixXd::Zero(p_num1d * m, p_num1d * m);
    MatrixXd Qi;
    for (int k = 0; k < m; k++)
    {
        Qi = MatrixXd::Zero(p_num1d, p_num1d);
        for (int i = 3; i < p_num1d; i++){
            for (int l = 3; l < p_num1d ; l++){
                Qi(i, l) = (i + 1) * i * (i - 1) * (i - 2) * (l + 1) * l * (l - 1) * (l - 2) / (i + l + 2 - 7) * pow(Time(k), i + l + 2 - 7);
            }
        }
        Q.block(k * p_num1d, k * p_num1d, p_num1d, p_num1d) = Qi;
    }
    // cout << "Q:" << endl
    //      << Q << endl;

    /*   Calculate PolyCoeff   */
    auto R = C * M.inverse().transpose() * Q * M.inverse() * Ct;

    int Fsize = Path.size() / 3 - 2 + 2 * 4;
    int Psize = R.rows() - Fsize;
    auto R_pp = R.block(Fsize, Fsize, Psize, Psize);
    auto R_fp = R.block(0, Fsize, Fsize, Psize);
    cout << "Path:" << endl
            << Path << endl;
    for (int i = 0; i < 3;i++){
        Vector4d start_cond(Path(0, i),Vel(0,i),Acc(1,i),0), end_cond(Path(Path.size() / 3-1, i),Vel(1,i),Acc(1,i),0);
        VectorXd dF(Fsize);
        dF << start_cond, Path.block(1,i,Path.size() / 3 - 2, 1), end_cond;
        auto dP = -R_pp.inverse() * R_fp.transpose() * dF;
        VectorXd d(Fsize + Psize);
        d << dF, dP;
        
        MatrixXd P1 = M.inverse() * Ct * d;
        MatrixXd P2 = P1.transpose();
        MatrixXd P(m, p_num1d);
        for (int j = 0; j < m; j++)
        {
            P.block(j,0,1,p_num1d) = P2.block(0, j * p_num1d, 1, p_num1d);
        }
        PolyCoeff.block(0, i * p_num1d, m, p_num1d) = P;
    }
    // cout << "PolyCoeff:" << endl
    //         << PolyCoeff << endl;
    cout << "PolyCoeff" << endl
         << PolyCoeff << endl;

    return PolyCoeff;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration_OOQP(
            const int d_order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time)
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    int nx = p_num1d * m;
    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */

    MatrixXd Q = MatrixXd::Zero(p_num1d * m, p_num1d * m);
    MatrixXd Qi;
    for (int k = 0; k < m; k++)
    {
        Qi = MatrixXd::Zero(p_num1d, p_num1d);
        for (int i = 3; i < p_num1d; i++){
            for (int l = 3; l < p_num1d ; l++){
                Qi(i, l) = (i + 1) * i * (i - 1) * (i - 2) * (l + 1) * l * (l - 1) * (l - 2) / (i + l + 2 - 7) * pow(Time(k), i + l + 2 - 7);
            }
        }
        Q.block(k * p_num1d, k * p_num1d, p_num1d, p_num1d) = Qi;
    }

    vector<pair<int, int>> index;
    vector<double> Qvalue;
    for (int i = 0; i < p_num1d * m; i++)
    {
        for (int j = 0; j <= i;j++){
            if(Q(i,j)!=0){
                index.push_back(make_pair(i,j));
                if(i==j)
                {
                    Qvalue.push_back(Q(i, j));
                }
                else
                {
                    Qvalue.push_back(Q(i, j) / 2);
                }
            }
        }
    }
    int nnzQ = Qvalue.size();
    int irowQ[nnzQ] = {0};
    int jcolQ[nnzQ] = {0};
    double dQ[nnzQ] = {0};
    for (int i = 0; i < nnzQ;i++){
        irowQ[i] = index[i].first;
        jcolQ[i] = index[i].second;
        dQ[i] = Qvalue[i];
    }
    double c[nx] = {0};
    // cout << "dQ[i]:" << dQ[0] << endl;
    for (int dim = 0; dim < 3;dim++){
        // int dim = 0;
        //Aeq b
        //start
        MatrixXd Aeq_start = MatrixXd::Zero(4, p_num1d * m);
        MatrixXd beq_start = MatrixXd::Zero(4, 1);
        Aeq_start(0, 0) = 1;//p
        Aeq_start(1, 1) = 1;//v
        Aeq_start(2, 2) = 2;//a
        Aeq_start(3, 3) = 6;//j
        beq_start(0) = Path(0, dim);
        beq_start(1) = Vel(0, dim);
        beq_start(2) = Acc(0, dim);
        beq_start(3) = 0;
        // cout << "Aeq_start:" << endl
        //      << Aeq_start << endl;

        //end
        MatrixXd Aeq_end = MatrixXd::Zero(4, p_num1d * m);
        MatrixXd beq_end = MatrixXd::Zero(4, 1);
        //p
        // cout << "***************"<< endl;
        for (int i = 0; i < p_num1d;i++){
            Aeq_end(0, p_num1d * m - p_num1d + i) = pow(Time(m-1), i);
            //cout << "***************"<< endl;
        }
        //cout << "***************"<< endl;
        //v
        for (int i = 1; i < p_num1d;i++){
            Aeq_end(1, p_num1d * m - p_num1d + i) = i * pow(Time(m-1), i - 1);
        }
        //a
        for (int i = 2; i < p_num1d;i++){
            Aeq_end(2, p_num1d * m - p_num1d + i) = i*(i-1)*pow(Time(m-1), i-2);
        }
        //j
        for (int i = 3; i < p_num1d;i++){
            Aeq_end(3, p_num1d * m - p_num1d + i) = i * (i - 1) * (i - 1) * pow(Time(m-1), i - 3);
        }
        beq_end(0) = Path(Path.rows() - 1, dim);
        beq_end(1) = Vel(1, dim);
        beq_end(2) = Acc(1, dim);
        beq_end(3) = 0;
        //cout << "----------------"<< endl;
        //position constrain
        MatrixXd Aeq_wp = MatrixXd::Zero(m-1, p_num1d * m);
        MatrixXd beq_wp = MatrixXd::Zero(m-1, 1);
        for (int i = 0; i < m - 1; i++)
        {
            for (int j = 0; j < p_num1d; j++)
            {
                Aeq_wp(i, i * p_num1d + j) = pow(Time(i), j);
            }
            beq_wp(i) = Path(i+1, dim);
        }
        //cout << "===================="<< endl;
        //position continuity constrain
        MatrixXd Aeq_con_p = MatrixXd::Zero(m-1, p_num1d * m);
        MatrixXd beq_con_p = MatrixXd::Zero(m-1, 1);
        for (int i = 0; i < m - 1; i++)
        {
            for (int j = 0; j < p_num1d; j++)
            {
                Aeq_con_p(i, i * p_num1d + j) = pow(Time(i), j);
            }
            Aeq_con_p(i, (i + 1) * p_num1d) = -1;
        }
        // cout << "Aeq_con_p:" << endl
        //      << Aeq_con_p << endl;
        MatrixXd Aeq_con_v = MatrixXd::Zero(m-1, p_num1d * m);
        MatrixXd beq_con_v = MatrixXd::Zero(m-1, 1);
        for (int i = 0; i < m - 1; i++)
        {
            for (int j = 1; j < p_num1d; j++)
            {
                Aeq_con_v(i, i * p_num1d + j) = j * pow(Time(i), j - 1);
            }
            //cout << "(i + 1) * p_num1d + 1:" << (i + 1) * p_num1d + 1 << endl;
            Aeq_con_v(i, (i + 1) * p_num1d + 1) = -1;
        }
        // cout << "Aeq_con_v:" << endl
        //      << Aeq_con_v << endl;
        //acceleration continuity constrain
        MatrixXd Aeq_con_a = MatrixXd::Zero(m-1, p_num1d * m);
        MatrixXd beq_con_a = MatrixXd::Zero(m-1, 1);
        for (int i = 0; i < m-1;i++)
        {
            for (int j = 2; j < p_num1d;j++){
                Aeq_con_a(i, i * p_num1d + j) =j*(j-1)*pow(Time(i), j-2);
            }
            Aeq_con_a(i, (i + 1) * p_num1d + 2) = -2;
        }
        // cout << "Aeq_con_a:" << endl
        //      << Aeq_con_a << endl;
        //jerk continuity constrain
        MatrixXd Aeq_con_j = MatrixXd::Zero(m-1, p_num1d * m);
        MatrixXd beq_con_j = MatrixXd::Zero(m-1, 1);
        for (int i = 0; i < m-1;i++)
        {
            for (int j = 3; j < p_num1d;j++){
                Aeq_con_j(i, i * p_num1d + j) =j*(j-1)*(j-2)*pow(Time(i), j-3);
            }
            Aeq_con_j(i, (i + 1) * p_num1d + 3) = -6;
        }
        // cout << "Aeq_con_j:" << endl
        //      << Aeq_con_j << endl;
        MatrixXd Aeq = MatrixXd::Zero(5 * (m - 1) + 2 * 4, p_num1d * m);
        MatrixXd beq = MatrixXd::Zero(5 * (m - 1) + 2 * 4, 1);
        Aeq.block(0, 0, 4, p_num1d * m) = Aeq_start;
        Aeq.block(4, 0, 4, p_num1d * m) = Aeq_end;
        Aeq.block(8, 0, m - 1, p_num1d * m) = Aeq_wp;
        Aeq.block(8 + 1 * (m - 1), 0, m - 1, p_num1d * m) = Aeq_con_p;
        Aeq.block(8 + 2 * (m - 1), 0, m - 1, p_num1d * m) = Aeq_con_v;
        Aeq.block(8 + 3 * (m - 1), 0, m - 1, p_num1d * m) = Aeq_con_a;
        Aeq.block(8 + 4 * (m - 1), 0, m - 1, p_num1d * m) = Aeq_con_j;

        beq.block(0, 0, 4, 1) = beq_start;
        beq.block(4, 0, 4, 1) = beq_end;
        beq.block(8, 0, m - 1, 1) = beq_wp;
        beq.block(8 + 1 * (m - 1), 0, m - 1, 1) = beq_con_p;
        beq.block(8 + 2 * (m - 1), 0, m - 1, 1) = beq_con_v;
        beq.block(8 + 3 * (m - 1), 0, m - 1, 1) = beq_con_a;
        beq.block(8 + 4 * (m - 1), 0, m - 1, 1) = beq_con_j;
        // cout << "Aeq:" << endl
        //      << Aeq << endl;
        // cout << "beq" << endl
        //      << beq << endl;

        index.clear();
        vector<double> Avalue;
        vector<double> bvalue;
        for (int i = 0; i < Aeq.rows(); i++)
        {

            for (int j = 0; j < Aeq.cols();j++){
                if(Aeq(i,j)!=0){
                    index.push_back(make_pair(i,j)); 
                    Avalue.push_back(Aeq(i, j));
                }
            }
            bvalue.push_back(beq(i, 0));
        }
        // cout << "-------------------" << endl;

        int my = bvalue.size();
        double b[my] = {0};
        for (int i = 0; i < my;i++){
            b[i] = bvalue[i];
        }

        int nnzA = Avalue.size();
        int irowA[nnzA] = {0};
        int jcolA[nnzA] = {0};
        double dA[nnzA] = {0};
        for (int i = 0; i < nnzA;i++){
            irowA[i] = index[i].first;
            jcolA[i] = index[i].second;
            dA[i] = Avalue[i];
        }
        // cout << "-------------------" << endl;
        double xupp[nx] = {0};  
        char ixupp[nx] = {0};

        double xlow[nx] = {0};
        char ixlow[nx] = {0};

        //neq constrain
        int mz = 0;
        double *clow = 0;
        char *iclow = 0;
        double *cupp = 0;
        char *icupp = 0;
        //neq 
        int nnzC = 0;
        int *irowC = 0;
        int *jcolC = 0;
        double *dC = 0;

        QpGenSparseMa27 *qp = new QpGenSparseMa27(nx, my, mz, nnzQ, nnzA, nnzC);
        // cout << "-------------------" << endl;
        QpGenData * prob = (QpGenData * ) qp->copyDataFromSparseTriple(
            c,      irowQ,  nnzQ,   jcolQ,  dQ,
            xlow,   ixlow,  xupp,   ixupp,
            irowA,  nnzA,   jcolA,  dA,     b,
            irowC,  nnzC,   jcolC,  dC,
            clow,   iclow,  cupp,   icupp );
        // cout << "-------------------" << endl;
        QpGenVars *vars = (QpGenVars *)qp->makeVariables(prob);
        // cout << "-------------------" << endl;
        QpGenResiduals *resid = (QpGenResiduals *)qp->makeResiduals(prob);
        // cout << "-------------------" << endl;
        GondzioSolver *s = new GondzioSolver(qp, prob);
        // cout << "-------------------" << endl;
        // if (!quiet)
        // s->monitorSelf();
        int ierr = s->solve(prob, vars, resid);

        // if (ierr == 0) {
        //     cout.precision(4);
        //     cout << "Solution: \n";
        //     vars->x->writefToStream(cout, "x[%{index}] = %{value}");
        // } else {
        //     cout << "Could not solve the problem.\n";
        // }
        double result[nx] = {0};
        vars->x->copyIntoArray(result);
        // double a = vars->x->copyIntoArray
        MatrixXd PolyCoeff_dim = MatrixXd::Zero(m, p_num1d);

        for (int i = 0; i < m;i++){
            for (int j = 0; j < p_num1d;j++){
                PolyCoeff_dim(i, j) = result[i * p_num1d + j];
                
            }
        }
        PolyCoeff.block(0, dim * p_num1d, m, p_num1d) = PolyCoeff_dim;
    }

    // MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);
    cout << "PolyCoeff_OOQP" << endl
         << PolyCoeff << endl;

    return PolyCoeff;
}