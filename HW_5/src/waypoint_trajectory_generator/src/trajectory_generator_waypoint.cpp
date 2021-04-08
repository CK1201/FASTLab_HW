#include "trajectory_generator_waypoint.h"
#include <stdio.h>
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
        P1 = P1.transpose();
        MatrixXd P(m, p_num1d);
        for (int j = 0; j < m; j++)
        {
            P.block(j,0,1,p_num1d) = P1.block(0, j * p_num1d, 1, p_num1d);
        }
        PolyCoeff.block(0, i * p_num1d, m, p_num1d) = P;
    }

    return PolyCoeff;
}

