// This file is part of mpc.

// mpc is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with mpc.  If not, see
// <http://www.gnu.org/licenses/>.

#include "constraints.h"
#include "PreviewSystem.h"
#include "previewController.h"
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

} // namespace boost

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

BOOST_PYTHON_MODULE(_mpc)
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
        .value("DEFAULT", SolverFlag::DEFAULT)
#ifdef EIGEN_LSSOL_FOUND
        .value("LSSOL", SolverFlag::LSSOL)
#endif
#ifdef EIGEN_GUROBI_FOUND
        .value("GUROBIDense", SolverFlag::GUROBIDense)
#endif
#ifdef EIGEN_QLD_FOUND
        .value("QLD", SolverFlag::QLD)
#endif
        .value("QuadProgDense", SolverFlag::QuadProgDense);

    // Access Solvers from python
    def("pythonSolverFactory", &pythonSolverFactory, return_value_policy<manage_new_object>(), "Return a solver corresponding to the giden flag");

    // Preview System
    void (PreviewSystem::*sys1)(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
        const Eigen::VectorXd&, const Eigen::VectorXd&, int)
        = &PreviewSystem::system;
    void (PreviewSystem::*sys2)(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
        const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, int)
        = &PreviewSystem::system;

    class_<PreviewSystem>("PreviewSystem", "The PreviewSystem is a read-write structure that holds all the information of the system.", no_init)
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
        .value("Constraint", ConstraintFlag::Constraint)
        .value("EqualityConstraint", ConstraintFlag::EqualityConstraint)
        .value("InequalityConstraint", ConstraintFlag::InequalityConstraint)
        .value("BoundConstraint", ConstraintFlag::BoundConstraint);

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
    class_<TrajectoryConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryConstraint", "Trajectory constraint. The object if instansiable through a NewTrajectoryConstraint function", no_init)
        .def("trajectory", &TrajectoryConstraint::reset);
    class_<ControlConstraint, boost::noncopyable, bases<EqIneqConstraint> >("ControlConstraint", "Control constraint. The object if instansiable through a NewControlConstraint function", no_init)
        .def("control", &ControlConstraint::reset);
    class_<TrajectoryBoundConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryBoundConstraint", "Trajectory Bound constraint. The object if instansiable through a NewTrajectoryBoundConstraint function", no_init)
        .def("trajectoryBound", &TrajectoryBoundConstraint::reset);
    class_<ControlBoundConstraint, boost::noncopyable, bases<Constraint> >("ControlBoundConstraint", "Control Bound constraint. The object if instansiable through a NewControlBoundConstraint function", no_init)
        .def("controlBound", &ControlBoundConstraint::reset)
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
        "MPC. This class runs the mpc with the desired QP and fills the PreviewSystem it is attach to", init<optional<SolverFlag> >())
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
        "Faster version of the FullType but neglect parts before final time", init<optional<SolverFlag> >())
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