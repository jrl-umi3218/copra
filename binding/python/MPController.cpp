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

std::shared_ptr<TrajectoryConstraint> NewTrajectoryConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint = true)
{
    return std::make_shared<TrajectoryConstraint>(E, f, isInequalityConstraint);
}

std::shared_ptr<ControlConstraint> NewControlConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint = true)
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

    def("NewPreviewSystem", &NewPreviewSystem);
    def("NewTrajectoryConstraint", &NewTrajectoryConstraint);
    def("NewControlConstraint", &NewControlConstraint);
    def("NewTrajectoryBoundConstraint", &NewTrajectoryBoundConstraint);
    def("NewControlBoundConstraint", &NewControlBoundConstraint);

    enum_<SolverFlag>("SolverFlag")
        .value("DEFAULT", SolverFlag::DEFAULT)
#ifdef LSSOL_SOLVER_FOUND
        .value("LSSOL", SolverFlag::LSSOL)
#endif
#ifdef GUROBI_SOLVER_FOUND
        .value("GUROBIDense", SolverFlag::GUROBIDense)
#endif
#ifdef QLD_SOLVER_FOUND
        .value("QLD", SolverFlag::QLD)
#endif
        .value("QuadProgDense", SolverFlag::QuadProgDense);

    // Solvers
    def("pythonSolverFactory", &pythonSolverFactory, return_value_policy<manage_new_object>());

    // Preview System
    void (PreviewSystem::*sys1)(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
        const Eigen::VectorXd&, const Eigen::VectorXd&, int)
        = &PreviewSystem::system;
    void (PreviewSystem::*sys2)(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
        const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, int)
        = &PreviewSystem::system;

    class_<PreviewSystem>("PreviewSystem", no_init)
        .def("system", sys1)
        .def("system", sys2)
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
    enum_<ConstraintFlag>("ConstraintFlag")
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
    class_<TrajectoryConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryConstraint", no_init);
    class_<ControlConstraint, boost::noncopyable, bases<EqIneqConstraint> >("ControlConstraint", no_init);
    class_<TrajectoryBoundConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryBoundConstraint", no_init);
    class_<ControlBoundConstraint, boost::noncopyable, bases<Constraint> >("ControlBoundConstraint", no_init)
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

        void weights(const Eigen::VectorXd& wx, const Eigen::VectorXd& wu)
        {
            if (override weights = this->get_override("weights"))
                weights(wx, wu);
            else
                MPCTypeFull::weights(wx, wu);
        }

        void default_weights(const Eigen::VectorXd& wx, const Eigen::VectorXd& wu)
        {
            this->MPCTypeFull::weights(wx, wu);
        }
    };

    // The default copy-ctor is implicitely deleted due to ctor overloading
    class_<MPCTypeFullWrap, boost::noncopyable>("MPCTypeFull",
        init<optional<SolverFlag> >())
        .def(init<const std::shared_ptr<PreviewSystem>&, optional<SolverFlag> >())
        .def("selectQPSolver", &MPCTypeFull::selectQPSolver)
        .def("initializeController", &MPCTypeFull::initializeController, &MPCTypeFullWrap::default_initializeController)
        .def("solve", &MPCTypeFull::solve)
        .def("solveTime", &MPCTypeFull::solveTime)
        .def("solveAndBuildTime", &MPCTypeFull::solveAndBuildTime)
        .def("weights", &MPCTypeFull::weights, &MPCTypeFullWrap::default_weights)
        .def("control", &MPCTypeFull::control, return_internal_reference<>())
        .def("trajectory", &MPCTypeFull::trajectory)
        .def("addConstraint", &MPCTypeFull::addConstraint)
        .def("resetConstraints", &MPCTypeFull::resetConstraints);

    //MPCTypeLast
    class_<MPCTypeLast, boost::noncopyable, bases<MPCTypeFull> >("MPCTypeLast",
        init<optional<SolverFlag> >())
        .def(init<const std::shared_ptr<PreviewSystem>&, optional<SolverFlag> >());

    //cpu_times
    using namespace boost::timer;
    class_<cpu_times>("cpu_times")
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