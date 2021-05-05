#ifndef _TRAJECTORY_GENERATOR_WAYPOINT_H_
#define _TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <vector>
// #include <ooqp/QpGenData.h>
// #include <ooqp/QpGenVars.h>
// #include <ooqp/QpGenResiduals.h>
// #include <ooqp/GondzioSolver.h>
// #include <ooqp/QpGenSparseMa27.h>


#include "ooqp/QpGenData.h"
#include "ooqp/QpGenVars.h"
#include "ooqp/QpGenResiduals.h"
#include "ooqp/GondzioSolver.h"
#include "ooqp/QpGenSparseMa27.h"

class TrajectoryGeneratorWaypoint {
    private:
		double _qp_cost;
		Eigen::MatrixXd _Q;
		Eigen::VectorXd _Px, _Py, _Pz;
    public:
        TrajectoryGeneratorWaypoint();

        ~TrajectoryGeneratorWaypoint();

        Eigen::MatrixXd PolyQPGeneration(
            const int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time);

        Eigen::MatrixXd PolyQPGeneration_OOQP(
            const int d_order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time);
        
        int Factorial(int x);
};
        

#endif
