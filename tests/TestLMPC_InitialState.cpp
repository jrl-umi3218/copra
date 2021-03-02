/*
 * Author: Niels Dehio
 */

#include "doctest.h"
#include "systems.h"
#include "tools.h"
#include <Eigen/Core>
#include <algorithm>
#include "LMPC.h"
#include "InitialStateLMPC.h"
#include "PreviewSystem.h"
#include "constraints.h"
#include "costFunctions.h"
#include "QuadProgSolver.h"
#ifdef EIGEN_QLD_FOUND
#include "QLDSolver.h"
#endif
#ifdef EIGEN_LSSOL_FOUND
#include "LSSOLSolver.h"
#endif
#ifdef EIGEN_GUROBI_FOUND
#include "GUROBISolver.h"
#endif
#ifdef EIGEN_OSQP_FOUND
#include "OSQPSolver.h"
#endif

#include <memory>
#include <numeric>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions> // for matrix exponential

TEST_CASE_FIXTURE(BoundedSystem, "INITIAL_STATE_LMPC")
{
    //define costs and constraints
    const int numCost = 1;
    const int numCstr = 1;
    const int s_dim = 2;
    const int u_dim = 1;
    const double factor = 100;
    Eigen::MatrixXd tcMat = Eigen::MatrixXd::Random(numCost,s_dim);
    Eigen::VectorXd tcVec = factor*Eigen::VectorXd::Random(numCost);
    Eigen::MatrixXd tacMat = Eigen::MatrixXd::Random(numCost,s_dim);
    Eigen::VectorXd tacVec = factor*Eigen::VectorXd::Random(numCost);
    Eigen::MatrixXd ccMat = Eigen::MatrixXd::Random(numCost,u_dim);
    Eigen::VectorXd ccVec = factor*Eigen::VectorXd::Random(numCost);
    Eigen::MatrixXd mcStateMat = Eigen::MatrixXd::Random(numCost,s_dim);
    Eigen::MatrixXd mcControlMat = Eigen::MatrixXd::Random(numCost,u_dim);
    Eigen::VectorXd mcVec = factor*Eigen::VectorXd::Random(numCost);

    Eigen::MatrixXd tcstrMat = Eigen::MatrixXd::Random(numCstr,s_dim);
    Eigen::VectorXd tcstrVec = factor*Eigen::VectorXd::Random(numCstr);
    Eigen::MatrixXd ccstrMat = Eigen::MatrixXd::Random(numCstr,u_dim);
    Eigen::VectorXd ccstrVec = factor*Eigen::VectorXd::Random(numCstr);
    Eigen::MatrixXd mcstrStateMat = Eigen::MatrixXd::Random(numCstr,s_dim);
    Eigen::MatrixXd mcstrControlMat = Eigen::MatrixXd::Random(numCstr,u_dim);
    Eigen::VectorXd mcstrVec = factor*Eigen::VectorXd::Random(numCstr);
    Eigen::VectorXd tbcstrVec = factor*Eigen::VectorXd::Ones(s_dim);
    Eigen::VectorXd cbcstrVec = factor*Eigen::VectorXd::Ones(u_dim);

    std::shared_ptr<copra::TrajectoryCost> tc;
    std::shared_ptr<copra::TargetCost> tac;
    std::shared_ptr<copra::ControlCost> cc;
    std::shared_ptr<copra::MixedCost> mc;
    std::shared_ptr<copra::TrajectoryConstraint> tcstr;
    std::shared_ptr<copra::ControlConstraint> ccstr;
    std::shared_ptr<copra::MixedConstraint> mcstr;
    std::shared_ptr<copra::TrajectoryBoundConstraint> tbcstr;
    std::shared_ptr<copra::ControlBoundConstraint> cbcstr;

    tc = std::make_shared<copra::TrajectoryCost>(tcMat, tcVec);
    tac = std::make_shared<copra::TargetCost>(tacMat, tacVec);
    cc = std::make_shared<copra::ControlCost>(ccMat, ccVec);
    mc = std::make_shared<copra::MixedCost>(mcStateMat, mcControlMat, mcVec);
    tcstr = std::make_shared<copra::TrajectoryConstraint>(tcstrMat, tcstrVec, /* isInequalityConstraint = */ true);
    ccstr = std::make_shared<copra::ControlConstraint>(ccstrMat, ccstrVec, /* isInequalityConstraint = */ true);
    mcstr = std::make_shared<copra::MixedConstraint>(mcstrStateMat, mcstrControlMat, mcstrVec, /* isInequalityConstraint = */ true);
    tbcstr = std::make_shared<copra::TrajectoryBoundConstraint>(-tbcstrVec, tbcstrVec);
    cbcstr = std::make_shared<copra::ControlBoundConstraint>(-cbcstrVec, cbcstrVec);

    tc->weight(1);
    tac->weight(1);
    cc->weight(1);
    mc->weight(1);

    //initialize preview-system and lmpc
    std::shared_ptr<copra::PreviewSystem> previewSystem_;
    const copra::SolverFlag solverFlag = copra::SolverFlag::QLD; //or use another flag
    previewSystem_ = std::make_shared<copra::PreviewSystem>();
    copra::InitialStateLMPC lmpc = copra::InitialStateLMPC(previewSystem_, solverFlag);
    // copra::LMPC lmpc = copra::LMPC(previewSystem_, solverFlag); //for comparison
    const int nbSteps_ = 10;
    const double dt = 0.001;
    Eigen::MatrixXd combi_c = Eigen::MatrixXd::Zero(s_dim+u_dim,s_dim+u_dim);
    Eigen::MatrixXd combi_d = Eigen::MatrixXd::Zero(s_dim+u_dim,s_dim+u_dim);
    combi_c.topRightCorner(s_dim,s_dim).setIdentity();
    combi_d = (combi_c*dt).exp(); //matrix exponential
    Eigen::VectorXd biasVector2;
    biasVector2 = Eigen::VectorXd::Zero(s_dim);
    Eigen::VectorXd s_init;
    s_init = Eigen::VectorXd::Random(s_dim);
    previewSystem_->system( combi_d.topLeftCorner(s_dim, s_dim), 
                            combi_d.topRightCorner(s_dim, u_dim), 
                            biasVector2, s_init, nbSteps_ );

    //specify and solve lmpc
    lmpc.initializeController(previewSystem_);
    lmpc.resetInitialStateBounds(Eigen::VectorXd::Zero(s_dim), Eigen::VectorXd::Zero(s_dim));

    lmpc.addCost(tc);
    lmpc.addCost(tac);
    // lmpc.addCost(cc);
    lmpc.addCost(mc);
    lmpc.addConstraint(tcstr);
    lmpc.addConstraint(ccstr);
    lmpc.addConstraint(mcstr);
    // lmpc.addConstraint(tbcstr);
    // lmpc.addConstraint(cbcstr);

    bool solutionFound = lmpc.solve();

    lmpc.removeCost(tc);
    lmpc.removeCost(tac);
    // lmpc.removeCost(cc);
    lmpc.removeCost(mc);
    lmpc.removeConstraint(tcstr);
    lmpc.removeConstraint(ccstr);
    lmpc.removeConstraint(mcstr);
    // lmpc.removeConstraint(tbcstr);
    // lmpc.removeConstraint(cbcstr);

    if(solutionFound){
        std::cout << "commands: " << lmpc.control().transpose() << std::endl;
        std::cout << "states: " << lmpc.trajectory().transpose() << std::endl;
    }
    else{
        std::cout << "no solution found" << std::endl;
    }
}