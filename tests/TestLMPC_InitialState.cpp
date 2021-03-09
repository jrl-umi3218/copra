/*
 * Author: Niels Dehio
 */

#include "InitialStateLMPC.h"
#include "PreviewSystem.h"
#include "QuadProgSolver.h"
#include "constraints.h"
#include "costFunctions.h"
#include "doctest.h"
#include "systems.h"
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

void run_comparison_test(const bool fullSizeEntry)
{
    //define costs and constraints
    const int numCost = 1;
    const int numCstr = 1;
    const int xDim = 2;
    const int uDim = 1;
    const double factor = 10;
    const int nbSteps_ = 10;
    int U;
    int X;
    if (fullSizeEntry) {
        //this option results for all costs and constraints in fullSizeEntry==true
        U = nbSteps_; //NOTE: nbSteps_ == previewSystem_->fullUDim
        X = nbSteps_ + 1; //NOTE: nbSteps_+1 == previewSystem_->fullXDim
    } else {
        //this option results for all costs and constraints in fullSizeEntry==false
        U = 1;
        X = 1;
    }

    Eigen::MatrixXd tcMat = Eigen::MatrixXd::Ones(numCost, xDim);
    Eigen::VectorXd tcVec = factor * Eigen::VectorXd::Ones(numCost * X);
    Eigen::MatrixXd tacMat = Eigen::MatrixXd::Ones(numCost, xDim);
    Eigen::VectorXd tacVec = factor * Eigen::VectorXd::Ones(numCost);
    Eigen::MatrixXd ccMat = Eigen::MatrixXd::Ones(numCost, uDim);
    Eigen::VectorXd ccVec = factor * Eigen::VectorXd::Ones(numCost * U);
    Eigen::MatrixXd mcStateMat = Eigen::MatrixXd::Ones(numCost, xDim);
    Eigen::MatrixXd mcControlMat = Eigen::MatrixXd::Ones(numCost, uDim);
    Eigen::VectorXd mcVec = factor * Eigen::VectorXd::Ones(numCost * U);
    Eigen::MatrixXd tcstrMat = Eigen::MatrixXd::Ones(numCstr, xDim);
    Eigen::VectorXd tcstrVec = factor * Eigen::VectorXd::Ones(numCstr * X);
    Eigen::MatrixXd ccstrMat = Eigen::MatrixXd::Ones(numCstr, uDim);
    Eigen::VectorXd ccstrVec = factor * Eigen::VectorXd::Ones(numCstr * U);
    Eigen::MatrixXd mcstrStateMat = Eigen::MatrixXd::Ones(numCstr, xDim);
    Eigen::MatrixXd mcstrControlMat = Eigen::MatrixXd::Ones(numCstr, uDim);
    Eigen::VectorXd mcstrVec = factor * Eigen::VectorXd::Ones(numCstr * U);
    Eigen::VectorXd tbcstrVecL = (-std::numeric_limits<double>::infinity()) * Eigen::VectorXd::Ones(xDim * X);
    Eigen::VectorXd tbcstrVecU = (+std::numeric_limits<double>::infinity()) * Eigen::VectorXd::Ones(xDim * X);
    // Eigen::VectorXd tbcstrVecL = (-1000)*Eigen::VectorXd::Ones(xDim * X); //TODO why does this fail?
    // Eigen::VectorXd tbcstrVecU = (+1000)*Eigen::VectorXd::Ones(xDim * X); //TODO why does this fail?
    Eigen::VectorXd cbcstrVecL = (-3) * Eigen::VectorXd::Ones(uDim * U);
    Eigen::VectorXd cbcstrVecU = (+3) * Eigen::VectorXd::Ones(uDim * U);

    std::shared_ptr<copra::TrajectoryCost> tcA;
    std::shared_ptr<copra::TargetCost> tacA;
    std::shared_ptr<copra::ControlCost> ccA;
    std::shared_ptr<copra::MixedCost> mcA;
    std::shared_ptr<copra::TrajectoryConstraint> tcstrA;
    std::shared_ptr<copra::ControlConstraint> ccstrA;
    std::shared_ptr<copra::MixedConstraint> mcstrA;
    std::shared_ptr<copra::TrajectoryBoundConstraint> tbcstrA;
    std::shared_ptr<copra::ControlBoundConstraint> cbcstrA;
    tcA = std::make_shared<copra::TrajectoryCost>(tcMat, tcVec);
    tacA = std::make_shared<copra::TargetCost>(tacMat, tacVec);
    ccA = std::make_shared<copra::ControlCost>(ccMat, ccVec);
    mcA = std::make_shared<copra::MixedCost>(mcStateMat, mcControlMat, mcVec);
    tcstrA = std::make_shared<copra::TrajectoryConstraint>(tcstrMat, tcstrVec, /* isInequalityConstraint = */ true);
    ccstrA = std::make_shared<copra::ControlConstraint>(ccstrMat, ccstrVec, /* isInequalityConstraint = */ true);
    mcstrA = std::make_shared<copra::MixedConstraint>(mcstrStateMat, mcstrControlMat, mcstrVec, /* isInequalityConstraint = */ true);
    tbcstrA = std::make_shared<copra::TrajectoryBoundConstraint>(tbcstrVecL, tbcstrVecU);
    cbcstrA = std::make_shared<copra::ControlBoundConstraint>(cbcstrVecL, cbcstrVecU);
    tcA->autoSpan();
    tacA->autoSpan();
    ccA->autoSpan();
    mcA->autoSpan();
    tcstrA->autoSpan();
    ccstrA->autoSpan();
    mcstrA->autoSpan();
    tbcstrA->autoSpan();
    cbcstrA->autoSpan();
    tcA->weight(1);
    tacA->weight(1);
    ccA->weight(1);
    mcA->weight(1);

    std::shared_ptr<copra::TrajectoryCost> tcB;
    std::shared_ptr<copra::TargetCost> tacB;
    std::shared_ptr<copra::ControlCost> ccB;
    std::shared_ptr<copra::MixedCost> mcB;
    std::shared_ptr<copra::TrajectoryConstraint> tcstrB;
    std::shared_ptr<copra::ControlConstraint> ccstrB;
    std::shared_ptr<copra::MixedConstraint> mcstrB;
    std::shared_ptr<copra::TrajectoryBoundConstraint> tbcstrB;
    std::shared_ptr<copra::ControlBoundConstraint> cbcstrB;
    tcB = std::make_shared<copra::TrajectoryCost>(tcMat, tcVec);
    tacB = std::make_shared<copra::TargetCost>(tacMat, tacVec);
    ccB = std::make_shared<copra::ControlCost>(ccMat, ccVec);
    mcB = std::make_shared<copra::MixedCost>(mcStateMat, mcControlMat, mcVec);
    tcstrB = std::make_shared<copra::TrajectoryConstraint>(tcstrMat, tcstrVec, /* isInequalityConstraint = */ true);
    ccstrB = std::make_shared<copra::ControlConstraint>(ccstrMat, ccstrVec, /* isInequalityConstraint = */ true);
    mcstrB = std::make_shared<copra::MixedConstraint>(mcstrStateMat, mcstrControlMat, mcstrVec, /* isInequalityConstraint = */ true);
    tbcstrB = std::make_shared<copra::TrajectoryBoundConstraint>(tbcstrVecL, tbcstrVecU);
    cbcstrB = std::make_shared<copra::ControlBoundConstraint>(cbcstrVecL, cbcstrVecU);
    tcB->autoSpan();
    tacB->autoSpan();
    ccB->autoSpan();
    mcB->autoSpan();
    tcstrB->autoSpan();
    ccstrB->autoSpan();
    mcstrB->autoSpan();
    tbcstrB->autoSpan();
    cbcstrB->autoSpan();
    tcB->weight(1);
    tacB->weight(1);
    ccB->weight(1);
    mcB->weight(1);

    //initialize preview-system and lmpc
    const Eigen::MatrixXd combi = Eigen::MatrixXd::Ones(xDim + uDim, xDim + uDim);
    const Eigen::VectorXd biasVector = Eigen::VectorXd::Zero(xDim);
    const Eigen::VectorXd s_init = Eigen::VectorXd::Zero(xDim);
    std::shared_ptr<copra::PreviewSystem> previewSystem;
    previewSystem = std::make_shared<copra::PreviewSystem>();
    previewSystem->system(combi.topLeftCorner(xDim, xDim),
        combi.topRightCorner(xDim, uDim),
        biasVector, s_init, nbSteps_);
    const copra::SolverFlag solverFlag = copra::SolverFlag::QLD; //or use another flag

    //specify and solve lmpc
    copra::LMPC lmpcA = copra::LMPC(previewSystem, solverFlag);
    lmpcA.initializeController(previewSystem);
    lmpcA.addCost(tcA);
    lmpcA.addCost(tacA);
    lmpcA.addCost(ccA);
    lmpcA.addCost(mcA);
    lmpcA.addConstraint(tcstrA);
    lmpcA.addConstraint(ccstrA);
    lmpcA.addConstraint(mcstrA);
    lmpcA.addConstraint(tbcstrA);
    lmpcA.addConstraint(cbcstrA);
    REQUIRE(lmpcA.solve());
    lmpcA.removeCost(tcA);
    lmpcA.removeCost(tacA);
    lmpcA.removeCost(ccA);
    lmpcA.removeCost(mcA);
    lmpcA.removeConstraint(tcstrA);
    lmpcA.removeConstraint(ccstrA);
    lmpcA.removeConstraint(mcstrA);
    lmpcA.removeConstraint(tbcstrA);
    lmpcA.removeConstraint(cbcstrA);

    copra::InitialStateLMPC lmpcB = copra::InitialStateLMPC(previewSystem, solverFlag);
    lmpcB.initializeController(previewSystem);
    lmpcB.addCost(tcB);
    lmpcB.addCost(tacB);
    lmpcB.addCost(ccB);
    lmpcB.addCost(mcB);
    lmpcB.addConstraint(tcstrB);
    lmpcB.addConstraint(ccstrB);
    lmpcB.addConstraint(mcstrB);
    lmpcB.addConstraint(tbcstrB);
    lmpcB.addConstraint(cbcstrB);
    REQUIRE(lmpcB.solve());
    lmpcB.removeCost(tcB);
    lmpcB.removeCost(tacB);
    lmpcB.removeCost(ccB);
    lmpcB.removeCost(mcB);
    lmpcB.removeConstraint(tcstrB);
    lmpcB.removeConstraint(ccstrB);
    lmpcB.removeConstraint(mcstrB);
    lmpcB.removeConstraint(tbcstrB);
    lmpcB.removeConstraint(cbcstrB);

    Eigen::MatrixXd lmpcA_Q = lmpcA.Q();
    Eigen::MatrixXd lmpcA_Aineq = lmpcA.Aineq();
    Eigen::MatrixXd lmpcA_Aeq = lmpcA.Aeq();
    Eigen::VectorXd lmpcA_c = lmpcA.c();
    Eigen::VectorXd lmpcA_bineq = lmpcA.bineq();
    Eigen::VectorXd lmpcA_beq = lmpcA.beq();
    Eigen::VectorXd lmpcA_lb = lmpcA.lb();
    Eigen::VectorXd lmpcA_ub = lmpcA.ub();

    Eigen::MatrixXd lmpcB_Q = lmpcB.Q();
    Eigen::MatrixXd lmpcB_Aineq = lmpcB.Aineq();
    Eigen::MatrixXd lmpcB_Aeq = lmpcB.Aeq();
    Eigen::VectorXd lmpcB_c = lmpcB.c();
    Eigen::VectorXd lmpcB_bineq = lmpcB.bineq();
    Eigen::VectorXd lmpcB_beq = lmpcB.beq();
    Eigen::VectorXd lmpcB_lb = lmpcB.lb();
    Eigen::VectorXd lmpcB_ub = lmpcB.ub();

    REQUIRE_EQ(lmpcA.nrEqConstr(), lmpcB.nrEqConstr());
    REQUIRE_EQ(lmpcA.nrIneqConstr(), lmpcB.nrIneqConstr());

    REQUIRE_LE((lmpcA_Q - lmpcB_Q).cwiseAbs().maxCoeff(), 1e-6);
    REQUIRE_LE((lmpcA_c - lmpcB_c).cwiseAbs().maxCoeff(), 1e-6);
    REQUIRE_LE((lmpcA_lb - lmpcB_lb).cwiseAbs().maxCoeff(), 1e-6);
    REQUIRE_LE((lmpcA_ub - lmpcB_ub).cwiseAbs().maxCoeff(), 1e-6);
    if (lmpcA.nrEqConstr()) {
        REQUIRE_LE((lmpcA_Aeq - lmpcB_Aeq).cwiseAbs().maxCoeff(), 1e-6);
        REQUIRE_LE((lmpcA_beq - lmpcB_beq).cwiseAbs().maxCoeff(), 1e-6);
    }
    if (lmpcA.nrEqConstr()) {
        REQUIRE_LE((lmpcA_Aineq - lmpcB_Aineq).cwiseAbs().maxCoeff(), 1e-6);
        REQUIRE_LE((lmpcA_bineq - lmpcB_bineq).cwiseAbs().maxCoeff(), 1e-6);
    }

    const Eigen::VectorXd commandsA = lmpcA.control();
    const Eigen::VectorXd statesA = lmpcA.trajectory();
    const Eigen::VectorXd initialStateA = statesA.head(xDim);
    const Eigen::VectorXd commandsB = lmpcB.control();
    const Eigen::VectorXd statesB = lmpcB.trajectory();
    const Eigen::VectorXd initialStateB = statesB.head(xDim);
    std::cout << "commandsA: \n " << commandsA.transpose() << std::endl;
    std::cout << "commandsB: \n " << commandsB.transpose() << std::endl;
    std::cout << "statesA: \n " << statesA.transpose() << std::endl;
    std::cout << "statesB: \n " << statesB.transpose() << std::endl;
    std::cout << "initialStateA: \n " << initialStateA.transpose() << std::endl;
    std::cout << "initialStateB: \n " << initialStateB.transpose() << std::endl;

    //check if initial state is within its bounds
    for (auto i = 0; i < xDim; ++i) {
        CHECK_LE(initialStateA(i), s_init(i) + 1e-6);
        CHECK_LE(s_init(i), initialStateA(i) + 1e-6);
        REQUIRE_LE(initialStateA(i), s_init(i) + 1e-6);
        REQUIRE_LE(s_init(i), initialStateA(i) + 1e-6);

        CHECK_LE(initialStateB(i), s_init(i) + 1e-6);
        CHECK_LE(s_init(i), initialStateB(i) + 1e-6);
        REQUIRE_LE(initialStateB(i), s_init(i) + 1e-6);
        REQUIRE_LE(s_init(i), initialStateB(i) + 1e-6);
    }
}

TEST_CASE_FIXTURE(BoundedSystem, "LMPC_AND_INITIAL-STATE-LMPC_COMPARISON")
{
    //NOTE: there are no equality constraints in this test
    run_comparison_test(true);
    run_comparison_test(false);
}

// ########################################################################
// ########################################################################
// ########################################################################

void run_optimization_test(const bool fullSizeEntry)
{
    //define costs and constraints
    const int numCost = 1;
    const int numCstr = 1;
    const int xDim = 2;
    const int uDim = 1;
    const double factor = 10;
    const int nbSteps_ = 10;
    int U;
    int X;
    if (fullSizeEntry) {
        //this option results for all costs and constraints in fullSizeEntry==true
        U = nbSteps_; //NOTE: nbSteps_ == previewSystem_->fullUDim
        X = nbSteps_ + 1; //NOTE: nbSteps_+1 == previewSystem_->fullXDim
    } else {
        //this option results for all costs and constraints in fullSizeEntry==false
        U = 1;
        X = 1;
    }

    Eigen::MatrixXd tcMat = Eigen::MatrixXd::Ones(numCost, xDim);
    Eigen::VectorXd tcVec = factor * Eigen::VectorXd::Ones(numCost * X);
    Eigen::MatrixXd tacMat = Eigen::MatrixXd::Ones(numCost, xDim);
    Eigen::VectorXd tacVec = factor * Eigen::VectorXd::Ones(numCost);
    Eigen::MatrixXd ccMat = Eigen::MatrixXd::Ones(numCost, uDim);
    Eigen::VectorXd ccVec = factor * Eigen::VectorXd::Ones(numCost * U);
    Eigen::MatrixXd mcStateMat = Eigen::MatrixXd::Ones(numCost, xDim);
    Eigen::MatrixXd mcControlMat = Eigen::MatrixXd::Ones(numCost, uDim);
    Eigen::VectorXd mcVec = factor * Eigen::VectorXd::Ones(numCost * U);
    Eigen::MatrixXd tcstrMat = Eigen::MatrixXd::Ones(numCstr, xDim);
    Eigen::VectorXd tcstrVec = factor * Eigen::VectorXd::Ones(numCstr * X);
    Eigen::MatrixXd ccstrMat = Eigen::MatrixXd::Ones(numCstr, uDim);
    Eigen::VectorXd ccstrVec = factor * Eigen::VectorXd::Ones(numCstr * U);
    Eigen::MatrixXd mcstrStateMat = Eigen::MatrixXd::Ones(numCstr, xDim);
    Eigen::MatrixXd mcstrControlMat = Eigen::MatrixXd::Ones(numCstr, uDim);
    Eigen::VectorXd mcstrVec = factor * Eigen::VectorXd::Ones(numCstr * U);
    Eigen::VectorXd tbcstrVecL = (-std::numeric_limits<double>::infinity()) * Eigen::VectorXd::Ones(xDim * X);
    Eigen::VectorXd tbcstrVecU = (+std::numeric_limits<double>::infinity()) * Eigen::VectorXd::Ones(xDim * X);
    // Eigen::VectorXd tbcstrVecL = (-1000)*Eigen::VectorXd::Ones(xDim * X); //TODO why does this fail?
    // Eigen::VectorXd tbcstrVecU = (+1000)*Eigen::VectorXd::Ones(xDim * X); //TODO why does this fail?
    Eigen::VectorXd cbcstrVecL = (-3) * Eigen::VectorXd::Ones(uDim * U);
    Eigen::VectorXd cbcstrVecU = (+3) * Eigen::VectorXd::Ones(uDim * U);

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
    tbcstr = std::make_shared<copra::TrajectoryBoundConstraint>(tbcstrVecL, tbcstrVecU);
    cbcstr = std::make_shared<copra::ControlBoundConstraint>(cbcstrVecL, cbcstrVecU);
    tc->autoSpan();
    tac->autoSpan();
    cc->autoSpan();
    mc->autoSpan();
    tcstr->autoSpan();
    ccstr->autoSpan();
    mcstr->autoSpan();
    tbcstr->autoSpan();
    cbcstr->autoSpan();
    tc->weight(1);
    tac->weight(1);
    cc->weight(1);
    mc->weight(1);

    //initialize preview-system and lmpc
    const Eigen::MatrixXd combi = Eigen::MatrixXd::Ones(xDim + uDim, xDim + uDim);
    const Eigen::VectorXd biasVector = Eigen::VectorXd::Zero(xDim);
    const Eigen::VectorXd s_init = Eigen::VectorXd::Zero(xDim);
    std::shared_ptr<copra::PreviewSystem> previewSystem;
    previewSystem = std::make_shared<copra::PreviewSystem>();
    previewSystem->system(combi.topLeftCorner(xDim, xDim),
        combi.topRightCorner(xDim, uDim),
        biasVector, s_init, nbSteps_);
    const copra::SolverFlag solverFlag = copra::SolverFlag::QLD; //or use another flag

    //specify and solve lmpc
    copra::InitialStateLMPC lmpc = copra::InitialStateLMPC(previewSystem, solverFlag);
    lmpc.initializeController(previewSystem);
    const Eigen::VectorXd initialState_lowerBound = (-1) * Eigen::VectorXd::Ones(xDim); //this allows optimization of the initial state
    const Eigen::VectorXd initialState_upperBound = (+1) * Eigen::VectorXd::Ones(xDim); //this allows prevents optimization of the initial state
    lmpc.resetInitialStateBounds(initialState_lowerBound, initialState_upperBound);
    Eigen::MatrixXd R = (1e-6) * Eigen::MatrixXd::Identity(xDim, xDim); // ensure that R is positive definite
    Eigen::VectorXd r = Eigen::VectorXd::Zero(xDim);
    lmpc.resetInitialStateCost(R, r);
    lmpc.addCost(tc);
    lmpc.addCost(tac);
    lmpc.addCost(cc);
    lmpc.addCost(mc);
    lmpc.addConstraint(tcstr);
    lmpc.addConstraint(ccstr);
    lmpc.addConstraint(mcstr);
    lmpc.addConstraint(tbcstr);
    lmpc.addConstraint(cbcstr);
    REQUIRE(lmpc.solve());
    lmpc.removeCost(tc);
    lmpc.removeCost(tac);
    lmpc.removeCost(cc);
    lmpc.removeCost(mc);
    lmpc.removeConstraint(tcstr);
    lmpc.removeConstraint(ccstr);
    lmpc.removeConstraint(mcstr);
    lmpc.removeConstraint(tbcstr);
    lmpc.removeConstraint(cbcstr);

    const Eigen::VectorXd commands = lmpc.control();
    const Eigen::VectorXd states = lmpc.trajectory();
    const Eigen::VectorXd initialState = states.head(xDim);
    std::cout << "commands: \n " << commands.transpose() << std::endl;
    std::cout << "states: \n " << states.transpose() << std::endl;
    std::cout << "initialState: \n " << initialState.transpose() << std::endl;

    //check if initial state is within its bounds
    for (auto i = 0; i < xDim; ++i) {
        CHECK_LE(initialState(i), initialState_upperBound(i) + 1e-6);
        CHECK_LE(initialState_lowerBound(i), initialState(i) + 1e-6);
        REQUIRE_LE(initialState(i), initialState_upperBound(i) + 1e-6);
        REQUIRE_LE(initialState_lowerBound(i), initialState(i) + 1e-6);
    }
}

TEST_CASE_FIXTURE(BoundedSystem, "INITIAL-STATE-OPTIMIZATION")
{
    //NOTE: there are no equality constraints in this test
    run_optimization_test(true);
    run_optimization_test(false);
}
