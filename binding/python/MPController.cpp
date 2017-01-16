// This file is part of ModelPreviewController.

// ModelPreviewController is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ModelPreviewController is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ModelPreviewController.  If not, see
// <http://www.gnu.org/licenses/>.

#include "Constraints.h"
#include "PreviewController.h"
#include "PreviewSystem.h"
#include "solverUtils.h"
#include <boost/python.hpp>
#include <boost/timer/timer.hpp>
#include <memory>
#include <string>

//TODO: Add python docstring

namespace boost {

// Make boost::python understand std::shared_ptr
template <class T>
T* get_pointer(std::shared_ptr<T> p)
{
    return p.get();
}
}

namespace mpc {

std::shared_ptr<PreviewSystem> NewPreviewSystem()
{
    return std::make_shared<PreviewSystem>();
}

std::shared_ptr<TrajectoryConstraint> NewTrajectoryConstraint1(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    return std::make_shared<TrajectoryConstraint>(E, f);
}

std::shared_ptr<TrajectoryConstraint> NewTrajectoryConstraint2(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint)
{
    return std::make_shared<TrajectoryConstraint>(E, f, isInequalityConstraint);
}

std::shared_ptr<ControlConstraint> NewControlConstraint1(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    return std::make_shared<ControlConstraint>(E, f);
}

std::shared_ptr<ControlConstraint> NewControlConstraint2(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint = true)
{
    return std::make_shared<ControlConstraint>(E, f, isInequalityConstraint);
}

std::shared_ptr<TrajectoryBoundConstraint> NewTrajectoryBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
{
    return std::make_shared<TrajectoryBoundConstraint>(lower, upper);
}

std::shared_ptr<ControlBoundConstraint> NewControlBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
{
    return std::make_shared<ControlBoundConstraint>(lower, upper);
}
} // namespace mpc

BOOST_PYTHON_MODULE(_mpcontroller)
{
    using namespace mpc;
    using namespace boost::python;

    def("NewPreviewSystem", &NewPreviewSystem, "Create a new instance of a PrevieSystem shared_ptr");
    def("NewTrajectoryConstraint", &NewTrajectoryConstraint1, "Create a new instance of a TrajectoryConstraint shared_ptr");
    def("NewTrajectoryConstraint", &NewTrajectoryConstraint2, "Create a new instance of a TrajectoryConstraint shared_ptr");
    def("NewControlConstraint", &NewControlConstraint1, "Create a new instance of a ControlConstraint shared_ptr");
    def("NewControlConstraint", &NewControlConstraint2, "Create a new instance of a ControlConstraint shared_ptr");
    def("NewTrajectoryBoundConstraint", &NewTrajectoryBoundConstraint, "Create a new instance of a TrajectoryBoundConstraint shared_ptr");
    def("NewControlBoundConstraint", &NewControlBoundConstraint, "Create a new instance of a ControlBoundConstraint shared_ptr");

    enum_<SolverFlag>("SolverFlag", "Flags to qp solver")
        .value("DEFAULT", SolverFlag::DEFAULT, "Default flag (QuadProgDense)")
#ifdef LSSOL_SOLVER_FOUND
        .value("LSSOL", SolverFlag::LSSOL, "LSSOL flag")
#endif
#ifdef GUROBI_SOLVER_FOUND
        .value("GUROBIDense", SolverFlag::GUROBIDense, "Dense version of Gurobi flag")
#endif
#ifdef QLD_SOLVER_FOUND
        .value("QLD", SolverFlag::QLD, "QLD flag")
#endif
        .value("QuadProgDense", SolverFlag::QuadProgDense, "Dense version of QuadProg flag");

    // Access Solvers from python
    def("pythonSolverFactory", &pythonSolverFactory, return_value_policy<manage_new_object>(), "Return a solver corresponding to the giden flag");

    // Preview System
    void (PreviewSystem::*sys1)(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
        const Eigen::VectorXd&, const Eigen::VectorXd&, int)
        = &PreviewSystem::system;
    void (PreviewSystem::*sys2)(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
        const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, int)
        = &PreviewSystem::system;

    class_<PreviewSystem>("PreviewSystem", no_init, "The PreviewSystem is a read-write structure that holds all the information of the system.")
        .def("system", sys1)
        .def("system", sys2)
        .def_readwrite("isUpdated", &PreviewSystem::isUpdated)
        .def_readwrite("nrStep", &PreviewSystem::nrStep)
        .def_readwrite("xDim", &PreviewSystem::xDim)
        .def_readwrite("uDim", &PreviewSystem::uDim)
        .def_readwrite("fullXDim", &PreviewSystem::fullXDim)
        .def_readwrite("fullUDim", &PreviewSystem::fullUDim)
        .def_readwrite("x0", &PreviewSystem::x0)
        .def_readwrite("xd", &PreviewSystem::xd)
        .def_readwrite("A", &PreviewSystem::A)
        .def_readwrite("B", &PreviewSystem::B)
        .def_readwrite("d", &PreviewSystem::d)
        .def_readwrite("Phi", &PreviewSystem::Phi)
        .def_readwrite("Psi", &PreviewSystem::Psi)
        .def_readwrite("xi", &PreviewSystem::xi);

    //Constraint
    enum_<ConstraintFlag>("ConstraintFlag", "Flag to constraint type")
        .value("Constraint", ConstraintFlag::Constraint, "Global constraint flag")
        .value("EqualityConstraint", ConstraintFlag::EqualityConstraint, "Flag to equality constraint")
        .value("InequalityConstraint", ConstraintFlag::InequalityConstraint, "Flag to inequality constraint")
        .value("BoundConstraint", ConstraintFlag::BoundConstraint, "flag to bound constraint");

    struct ConstraintWrap : Constraint, wrapper<Constraint> {
        using Constraint::Constraint;

        void initializeConstraint(const PreviewSystem& ps)
        {
            this->get_override("initializeConstraint")(ps);
        }

        void update(const PreviewSystem& ps)
        {
            this->get_override("update")(ps);
        }

        ConstraintFlag constraintType()
        {
            return this->get_override("constraintType")();
        }
    };

    class_<ConstraintWrap, boost::noncopyable>("Constraint", init<const std::string&>())
        .def("initializeConstraint", pure_virtual(&Constraint::initializeConstraint))
        .def("update", pure_virtual(&Constraint::update))
        .def("constraintType", pure_virtual(&Constraint::constraintType))
        .def("name", &Constraint::name, return_internal_reference<>())
        .def("nrConstr", &Constraint::nrConstr);

    struct EqIneqConstraintWrap : EqIneqConstraint, wrapper<EqIneqConstraint> {
        using EqIneqConstraint::EqIneqConstraint;

        void initializeConstraint(const PreviewSystem& ps)
        {
            this->get_override("initializeConstraint")(ps);
        }

        void update(const PreviewSystem& ps)
        {
            this->get_override("update")(ps);
        }

        ConstraintFlag constraintType()
        {
            return this->get_override("constraintType")();
        }
    };

    class_<EqIneqConstraintWrap, boost::noncopyable, bases<Constraint> >("EqIneqConstraint", init<const std::string&, bool>())
        .def("A", &EqIneqConstraint::A, return_internal_reference<>())
        .def("b", &EqIneqConstraint::b, return_internal_reference<>());

    // Delete constructor to enforce call of New<Name_of_constraint> function that return a shared_ptr.
    class_<TrajectoryConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryConstraint", no_init, "Trajectory constraint. The object if instansiable through a NewTrajectoryConstraint function")
        .def("trajectory", &TrajectoryConstraint::trajectory);
    class_<ControlConstraint, boost::noncopyable, bases<EqIneqConstraint> >("ControlConstraint", no_init, "Control constraint. The object if instansiable through a NewControlConstraint function")
        .def("control", &ControlConstraint::control);
    class_<TrajectoryBoundConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryBoundConstraint", no_init, "Trajectory Bound constraint. The object if instansiable through a NewTrajectoryBoundConstraint function")
        .def("trajectoryBound", &TrajectoryBoundConstraint::trajectoryBound);
    class_<ControlBoundConstraint, boost::noncopyable, bases<Constraint> >("ControlBoundConstraint", no_init, "Control Bound constraint. The object if instansiable through a NewControlBoundConstraint function")
        .def("controlBound", &ControlBoundConstraint::controlBound)
        .def("lower", &ControlBoundConstraint::lower, return_internal_reference<>())
        .def("upper", &ControlBoundConstraint::upper, return_internal_reference<>());

    //MPCTypeFull
    struct MPCTypeFullWrap : MPCTypeFull, wrapper<MPCTypeFull> {
        using MPCTypeFull::MPCTypeFull;

        void initializeController(const std::shared_ptr<PreviewSystem>& ps)
        {
            if (override initializeController = this->get_override("initializeController"))
                initializeController(ps);
            else
                MPCTypeFull::initializeController(ps);
        }

        void default_initializeController(const std::shared_ptr<PreviewSystem>& ps)
        {
            this->MPCTypeFull::initializeController(ps);
        }

        void eigenWeights(const Eigen::VectorXd& wx, const Eigen::VectorXd& wu)
        {
            if (override weights = this->get_override("weights"))
                weights(wx, wu);
            else
                MPCTypeFull::weights(wx, wu);
        }

        void doubleWeight(double wx, double wu)
        {
            this->MPCTypeFull::weights(wx, wu);
        }

        void default_eigenWeights(const Eigen::VectorXd& wx, const Eigen::VectorXd& wu)
        {
            this->MPCTypeFull::weights(wx, wu);
        }
    };

    // The default copy-ctor is implicitely deleted due to ctor overloading
    class_<MPCTypeFullWrap, boost::noncopyable>("MPCTypeFull",
        init<optional<SolverFlag> >(), "MPC. This class runs the mpc with the desired QP and fills the PreviewSystem it is attach to")
        .def(init<const std::shared_ptr<PreviewSystem>&, optional<SolverFlag> >())
        .def("selectQPSolver", &MPCTypeFull::selectQPSolver)
        .def("initializeController", &MPCTypeFull::initializeController, &MPCTypeFullWrap::default_initializeController)
        .def("solve", &MPCTypeFull::solve)
        .def("solveTime", &MPCTypeFull::solveTime)
        .def("solveAndBuildTime", &MPCTypeFull::solveAndBuildTime)
        .def("weights", &MPCTypeFullWrap::eigenWeights, &MPCTypeFullWrap::default_eigenWeights)
        .def("weights", &MPCTypeFullWrap::doubleWeight)
        .def("control", &MPCTypeFull::control, return_internal_reference<>())
        .def("trajectory", &MPCTypeFull::trajectory)
        .def("addConstraint", &MPCTypeFull::addConstraint)
        .def("resetConstraints", &MPCTypeFull::resetConstraints);

    //MPCTypeLast
    struct MPCTypeLastWrap : MPCTypeLast, wrapper<MPCTypeLast> {
        using MPCTypeLast::MPCTypeLast;

        void eigenWeights(const Eigen::VectorXd& wx, const Eigen::VectorXd& wu)
        {
            if (override weights = this->get_override("weights"))
                weights(wx, wu);
            else
                MPCTypeLast::weights(wx, wu);
        }

        void doubleWeight(double wx, double wu)
        {
            this->MPCTypeLast::weights(wx, wu);
        }

        void default_eigenWeights(const Eigen::VectorXd& wx, const Eigen::VectorXd& wu)
        {
            this->MPCTypeLast::weights(wx, wu);
        }
    };

    class_<MPCTypeLastWrap, boost::noncopyable, bases<MPCTypeFull> >("MPCTypeLast",
        init<optional<SolverFlag> >(), "Faster version of the FullType but neglect parts before final time")
        .def(init<const std::shared_ptr<PreviewSystem>&, optional<SolverFlag> >())
        .def("weights", &MPCTypeLastWrap::eigenWeights, &MPCTypeLastWrap::default_eigenWeights)
        .def("weights", &MPCTypeLastWrap::doubleWeight);

    //cpu_times
    using namespace boost::timer;
    class_<cpu_times>("cpu_times", "Allow the use of boost measuring time")
        .def_readwrite("wall", &cpu_times::wall)
        .def_readwrite("user", &cpu_times::user)
        .def_readwrite("system", &cpu_times::system)
        .def("clear", &cpu_times::clear);

    register_ptr_to_python<std::shared_ptr<PreviewSystem> >();
    register_ptr_to_python<std::shared_ptr<TrajectoryConstraint> >();
    register_ptr_to_python<std::shared_ptr<ControlConstraint> >();
    register_ptr_to_python<std::shared_ptr<TrajectoryBoundConstraint> >();
    register_ptr_to_python<std::shared_ptr<ControlBoundConstraint> >();
    implicitly_convertible<std::shared_ptr<TrajectoryConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<ControlConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<TrajectoryBoundConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<ControlBoundConstraint>, std::shared_ptr<Constraint> >();
}