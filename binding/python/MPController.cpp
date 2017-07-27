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

#include "LMPC.h"
#include "AutoSpan.h"
#include "PreviewSystem.h"
#include "constraints.h"
#include "costFunctions.h"
#include "solverUtils.h"
#include <boost/python.hpp>
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

template <typename T, typename... Args>
std::shared_ptr<T> createSharedPointer(Args... args)
{
    return std::make_shared<T>(args...);
}

auto NewPreviewSystem1 = &createSharedPointer<PreviewSystem>;
auto NewPreviewSystem2 = &createSharedPointer<PreviewSystem, const Eigen::MatrixXd&, const Eigen::MatrixXd&,
    const Eigen::VectorXd&, const Eigen::VectorXd&, int>;
auto NewTrajectoryConstraint1 = &createSharedPointer<TrajectoryConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto NewTrajectoryConstraint2 = &createSharedPointer<TrajectoryConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&, bool>;
auto NewControlConstraint1 = &createSharedPointer<ControlConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto NewControlConstraint2 = &createSharedPointer<ControlConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&, bool>;
auto NewMixedConstraint1 = &createSharedPointer<MixedConstraint, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto NewMixedConstraint2 = &createSharedPointer<MixedConstraint, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&, bool>;
auto NewTrajectoryBoundConstraint = &createSharedPointer<TrajectoryBoundConstraint, const Eigen::VectorXd&, const Eigen::VectorXd&>;
auto NewControlBoundConstraint = &createSharedPointer<ControlBoundConstraint, const Eigen::VectorXd&, const Eigen::VectorXd&>;
auto NewTrajectoryCost = &createSharedPointer<TrajectoryCost, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto NewTargetCost = &createSharedPointer<TargetCost, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto NewControlCost = &createSharedPointer<ControlCost, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto NewMixedCost = &createSharedPointer<MixedCost, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
} // namespace mpc

BOOST_PYTHON_MODULE(_mpc)
{
    using namespace mpc;
    using namespace boost::python;

    def("NewPreviewSystem", NewPreviewSystem1, "Create a new instance of a PrevieSystem shared_ptr with default constructor");
    def("NewPreviewSystem", NewPreviewSystem2, "Create a new instance of a PrevieSystem shared_ptr");
    def("NewTrajectoryConstraint", NewTrajectoryConstraint1, "Create a new instance of a TrajectoryConstraint shared_ptr");
    def("NewTrajectoryConstraint", NewTrajectoryConstraint2, "Create a new instance of a TrajectoryConstraint shared_ptr");
    def("NewControlConstraint", NewControlConstraint1, "Create a new instance of a ControlConstraint shared_ptr");
    def("NewControlConstraint", NewControlConstraint2, "Create a new instance of a ControlConstraint shared_ptr");
    def("NewMixedConstraint", NewMixedConstraint1, "Create a new instance of a MixedConstraint shared_ptr");
    def("NewMixedConstraint", NewMixedConstraint2, "Create a new instance of a MixedConstraint shared_ptr");
    def("NewTrajectoryBoundConstraint", NewTrajectoryBoundConstraint, "Create a new instance of a TrajectoryBoundConstraint shared_ptr");
    def("NewControlBoundConstraint", NewControlBoundConstraint, "Create a new instance of a ControlBoundConstraint shared_ptr");
    def("NewTrajectoryCost", NewTrajectoryCost, "Create a new instance of a TrajectoryCost shared_ptr");
    def("NewTargetCost", NewTargetCost, "Create a new instance of a TargetCost shared_ptr");
    def("NewControlCost", NewControlCost, "Create a new instance of a ControlCost shared_ptr");
    def("NewMixedCost", NewMixedCost, "Create a new instance of a NewMixedCost shared_ptr");

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
    class_<PreviewSystem>("PreviewSystem", "The PreviewSystem is a read-write structure that holds all the information of the system.", no_init)
        .def("system", &PreviewSystem::system)
        .def("updateSystem", &PreviewSystem::updateSystem)
        .def_readwrite("isUpdated", &PreviewSystem::isUpdated)
        .def_readwrite("nrUStep", &PreviewSystem::nrUStep)
        .def_readwrite("nrXStep", &PreviewSystem::nrXStep)
        .def_readwrite("xDim", &PreviewSystem::xDim)
        .def_readwrite("uDim", &PreviewSystem::uDim)
        .def_readwrite("fullXDim", &PreviewSystem::fullXDim)
        .def_readwrite("fullUDim", &PreviewSystem::fullUDim)
        .def_readwrite("x0", &PreviewSystem::x0)
        .def_readwrite("A", &PreviewSystem::A)
        .def_readwrite("B", &PreviewSystem::B)
        .def_readwrite("d", &PreviewSystem::d)
        .def_readwrite("Phi", &PreviewSystem::Phi)
        .def_readwrite("Psi", &PreviewSystem::Psi)
        .def_readwrite("xi", &PreviewSystem::xi);

    //AutoSpan
    class_<AutoSpan, boost::noncopyable>("AutoSpan", "Helper functions to automatically extend a matrix to a desired dimension.", no_init)
        .def("spanMatrix", &AutoSpan::spanMatrix).staticmethod("spanMatrix")
        .def("spanVector", &AutoSpan::spanVector).staticmethod("spanVector");

    //Constraint
    enum_<ConstraintFlag>("ConstraintFlag", "Flag to constraint type")
        .value("Constraint", ConstraintFlag::Constraint)
        .value("EqualityConstraint", ConstraintFlag::EqualityConstraint)
        .value("InequalityConstraint", ConstraintFlag::InequalityConstraint)
        .value("BoundConstraint", ConstraintFlag::BoundConstraint);

    struct ConstraintWrap : Constraint, wrapper<Constraint> {
        using Constraint::Constraint;

        void autoSpan()
        {
            this->get_override("autoSpan")();
        }

        void initializeConstraint(const PreviewSystem& ps)
        {
            this->get_override("initializeConstraint")(ps);
        }

        void update(const PreviewSystem& ps)
        {
            this->get_override("update")(ps);
        }

        ConstraintFlag constraintType() const noexcept
        {
            return this->get_override("constraintType")();
        }
    };

    class_<ConstraintWrap, boost::noncopyable>("Constraint", no_init) // Disable the constructor because of move semantics. No need anyway.
        .def("autoSpan", pure_virtual(&Constraint::autoSpan))
        .def("initializeConstraint", pure_virtual(&Constraint::initializeConstraint))
        .def("update", pure_virtual(&Constraint::update))
        .def("constraintType", pure_virtual(&Constraint::constraintType))
        .def("name", &Constraint::name, return_internal_reference<>())
        .def("nrConstr", &Constraint::nrConstr);

    struct EqIneqConstraintWrap : EqIneqConstraint, wrapper<EqIneqConstraint> {
        using EqIneqConstraint::EqIneqConstraint;

        void autoSpan()
        {
            this->get_override("autoSpan")();
        }

        void initializeConstraint(const PreviewSystem& ps)
        {
            this->get_override("initializeConstraint")(ps);
        }

        void update(const PreviewSystem& ps)
        {
            this->get_override("update")(ps);
        }

        ConstraintFlag constraintType() const noexcept
        {
            return this->get_override("constraintType")();
        }
    };

    class_<EqIneqConstraintWrap, boost::noncopyable, bases<Constraint> >("EqIneqConstraint", init<const std::string&, bool>())
        .def("A", &EqIneqConstraint::A, return_internal_reference<>())
        .def("b", &EqIneqConstraint::b, return_internal_reference<>());

    // Delete constructor to enforce call of New<Name_of_constraint> function that return a shared_ptr.
    class_<TrajectoryConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryConstraint", "Trajectory constraint. The object is instansiable through a NewTrajectoryConstraint function", no_init);
    class_<ControlConstraint, boost::noncopyable, bases<EqIneqConstraint> >("ControlConstraint", "Control constraint. The object is instansiable through a NewControlConstraint function", no_init);
    class_<MixedConstraint, boost::noncopyable, bases<EqIneqConstraint> >("MixedConstraint", "Mixed constraint. The object is instansiable through a NewMixedConstraint function", no_init);
    class_<TrajectoryBoundConstraint, boost::noncopyable, bases<EqIneqConstraint> >("TrajectoryBoundConstraint", "Trajectory Bound constraint. The object is instansiable through a NewTrajectoryBoundConstraint function", no_init);
    class_<ControlBoundConstraint, boost::noncopyable, bases<Constraint> >("ControlBoundConstraint", "Control Bound constraint. The object is instansiable through a NewControlBoundConstraint function", no_init)
        .def("lower", &ControlBoundConstraint::lower, return_internal_reference<>())
        .def("upper", &ControlBoundConstraint::upper, return_internal_reference<>());

    //Cost Functions
    struct CostFunctionWrap : CostFunction, wrapper<CostFunction> {
        using CostFunction::CostFunction;

        void autoSpan()
        {
            if (override autoSpan = this->get_override("autoSpan"))
                autoSpan();
            else
                CostFunction::autoSpan();
        }

        void default_autoSpan()
        {
            this->get_override("autoSpan");
        }

        void initializeCost(const PreviewSystem& ps)
        {
            if (override initializeConstraint = this->get_override("initializeCost"))
                initializeCost(ps);
            else
                CostFunction::initializeCost(ps);
        }

        void default_initializeCost(const PreviewSystem& ps)
        {
            this->get_override("initializeCost")(ps);
        }

        ConstraintFlag constraintType() const noexcept
        {
            return this->get_override("constraintType")();
        }

        void update(const PreviewSystem& ps)
        {
            this->get_override("update")(ps);
        }
    };

    class_<CostFunctionWrap, boost::noncopyable>("CostFunction", no_init) // Disable the constructor because of move semantics. No need anyway.
        .def("initializeConstraint", &CostFunction::initializeCost, &CostFunctionWrap::default_initializeCost)
        .def("update", pure_virtual(&CostFunction::update))
        .def("name", &CostFunction::name, return_internal_reference<>())
        .def("Q", &CostFunction::Q, return_internal_reference<>())
        .def("c", &CostFunction::c, return_internal_reference<>())
        .def("weights", &CostFunction::weights<const Eigen::VectorXd&>)
        .def("weights", &CostFunction::weights<double>);

    // Delete constructor to enforce call of New<Name_of_cost> function that return a shared_ptr.
    class_<TrajectoryCost, boost::noncopyable, bases<CostFunction> >("TrajectoryCost", "Trajectory cost. The object is instansiable through a NewTrajectoryCost function", no_init);
    class_<TargetCost, boost::noncopyable, bases<CostFunction> >("TargetCost", "Target cost. The object is instansiable through a NewTargetCost function", no_init);
    class_<ControlCost, boost::noncopyable, bases<CostFunction> >("ControlCost", "Control cost. The object is instansiable through a NewControlCost function", no_init);
    class_<MixedCost, boost::noncopyable, bases<CostFunction> >("MixedCost", "Mixed cost. The object is instansiable through a NewMixedCost function", no_init);

    //LMPC
    class_<LMPC, boost::noncopyable>("LMPC",
        "LMPC. This class runs the mpc with the desired QP and fills the PreviewSystem it is attach to", init<optional<SolverFlag> >())
        .def(init<const std::shared_ptr<PreviewSystem>&, optional<SolverFlag> >())
        .def("selectQPSolver", &LMPC::selectQPSolver)
        .def("initializeController", &LMPC::initializeController)
        .def("solve", &LMPC::solve)
        .def("solveTime", &LMPC::solveTime)
        .def("solveAndBuildTime", &LMPC::solveAndBuildTime)
        .def("control", &LMPC::control, return_internal_reference<>())
        .def("trajectory", &LMPC::trajectory)
        .def("addCost", &LMPC::addCost)
        .def("addConstraint", &LMPC::addConstraint)
        .def("resetConstraints", &LMPC::resetConstraints);

    // Register pointer
    register_ptr_to_python<std::shared_ptr<ControlBoundConstraint> >();
    register_ptr_to_python<std::shared_ptr<ControlConstraint> >();
    register_ptr_to_python<std::shared_ptr<ControlCost> >();
    register_ptr_to_python<std::shared_ptr<MixedConstraint> >();
    register_ptr_to_python<std::shared_ptr<MixedCost> >();
    register_ptr_to_python<std::shared_ptr<PreviewSystem> >();
    register_ptr_to_python<std::shared_ptr<TargetCost> >();
    register_ptr_to_python<std::shared_ptr<TrajectoryBoundConstraint> >();
    register_ptr_to_python<std::shared_ptr<TrajectoryConstraint> >();
    register_ptr_to_python<std::shared_ptr<TrajectoryCost> >();

    // Implicit conversion of pointers
    implicitly_convertible<std::shared_ptr<ControlBoundConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<ControlConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<ControlCost>, std::shared_ptr<CostFunction> >();
    implicitly_convertible<std::shared_ptr<MixedConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<MixedCost>, std::shared_ptr<CostFunction> >();
    implicitly_convertible<std::shared_ptr<TargetCost>, std::shared_ptr<CostFunction> >();
    implicitly_convertible<std::shared_ptr<TrajectoryBoundConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<TrajectoryConstraint>, std::shared_ptr<Constraint> >();
    implicitly_convertible<std::shared_ptr<TrajectoryCost>, std::shared_ptr<CostFunction> >();
}