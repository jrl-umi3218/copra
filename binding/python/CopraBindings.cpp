/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "AutoSpan.h"
#include "LMPC.h"
#include "PreviewSystem.h"
#include "constraints.h"
#include "costFunctions.h"
#include "solverUtils.h"
#include <boost/python.hpp>
#include <memory>
#include <pygen/converters.h>
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

namespace copra {

template <typename T, typename... Args>
std::shared_ptr<T> createSharedPointer(Args... args)
{
    return std::make_shared<T>(args...);
}

auto PreviewSystem1 = &createSharedPointer<PreviewSystem>;
auto PreviewSystem2 = &createSharedPointer<PreviewSystem, const Eigen::MatrixXd&, const Eigen::MatrixXd&,
    const Eigen::VectorXd&, const Eigen::VectorXd&, int>;
auto TrajectoryConstraint1 = &createSharedPointer<TrajectoryConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto TrajectoryConstraint2 = &createSharedPointer<TrajectoryConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&, bool>;
auto ControlConstraint1 = &createSharedPointer<ControlConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto ControlConstraint2 = &createSharedPointer<ControlConstraint, const Eigen::MatrixXd&, const Eigen::VectorXd&, bool>;
auto MixedConstraint1 = &createSharedPointer<MixedConstraint, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto MixedConstraint2 = &createSharedPointer<MixedConstraint, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&, bool>;
auto TrajectoryBoundConstraint1 = &createSharedPointer<TrajectoryBoundConstraint, const Eigen::VectorXd&, const Eigen::VectorXd&>;
auto ControlBoundConstraint1 = &createSharedPointer<ControlBoundConstraint, const Eigen::VectorXd&, const Eigen::VectorXd&>;
auto TrajectoryCost1 = &createSharedPointer<TrajectoryCost, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto TargetCost1 = &createSharedPointer<TargetCost, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto ControlCost1 = &createSharedPointer<ControlCost, const Eigen::MatrixXd&, const Eigen::VectorXd&>;
auto MixedCost1 = &createSharedPointer<MixedCost, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&>;

// AutoSpan wrapper functions
Eigen::MatrixXd autoSpanMatrix0(Eigen::MatrixXd mat, Eigen::Index new_dim)
{
    AutoSpan::spanMatrix(mat, new_dim);
    return mat;
}

Eigen::MatrixXd autoSpanMatrix1(Eigen::MatrixXd mat, Eigen::Index new_dim, int addCols)
{
    AutoSpan::spanMatrix(mat, new_dim, addCols);
    return mat;
}

Eigen::MatrixXd autoSpanVector(Eigen::VectorXd vec, Eigen::Index new_dim)
{
    AutoSpan::spanVector(vec, new_dim);
    return vec;
}

} // namespace copra

BOOST_PYTHON_MODULE(pyCopra)
{
    using namespace copra;
    namespace py = boost::python;
    namespace np = boost::python::numpy;

    Py_Initialize();
    np::initialize();
    pygen::convertMatrix<Eigen::MatrixXd>();
    pygen::convertVector<Eigen::VectorXd>();

    py::enum_<SolverFlag>("SolverFlag", "Flags to qp solver")
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
    py::def("pythonSolverFactory", &pythonSolverFactory, py::return_value_policy<py::manage_new_object>(), "Return a solver corresponding to the giden flag");

    // Preview System
    py::class_<PreviewSystem, std::shared_ptr<PreviewSystem>>("PreviewSystem", "The PreviewSystem is a read-write structure that holds all the information of the system.", py::no_init)
        .def("__init__", py::make_constructor(PreviewSystem1))
        .def("__init__", py::make_constructor(PreviewSystem2))
        .def("system", &PreviewSystem::system)
        .def("update_system", &PreviewSystem::updateSystem)
        .def_readwrite("is_updated", &PreviewSystem::isUpdated)
        .def_readwrite("nr_u_Step", &PreviewSystem::nrUStep)
        .def_readwrite("nr_x_Step", &PreviewSystem::nrXStep)
        .def_readwrite("x_dim", &PreviewSystem::xDim)
        .def_readwrite("u_dim", &PreviewSystem::uDim)
        .def_readwrite("full_x_dim", &PreviewSystem::fullXDim)
        .def_readwrite("full_u_dim", &PreviewSystem::fullUDim)
        .add_property("x0", py::make_getter(&PreviewSystem::x0, py::return_value_policy<py::copy_non_const_reference>()), py::make_setter(&PreviewSystem::x0))
        .add_property("A", py::make_getter(&PreviewSystem::A, py::return_value_policy<py::copy_non_const_reference>()), py::make_setter(&PreviewSystem::A))
        .add_property("B", py::make_getter(&PreviewSystem::B, py::return_value_policy<py::copy_non_const_reference>()), py::make_setter(&PreviewSystem::B))
        .add_property("d", py::make_getter(&PreviewSystem::d, py::return_value_policy<py::copy_non_const_reference>()), py::make_setter(&PreviewSystem::d))
        .add_property("Phi", py::make_getter(&PreviewSystem::Phi, py::return_value_policy<py::copy_non_const_reference>()), py::make_setter(&PreviewSystem::Phi))
        .add_property("Psi", py::make_getter(&PreviewSystem::Psi, py::return_value_policy<py::copy_non_const_reference>()), py::make_setter(&PreviewSystem::Psi))
        .add_property("xi", py::make_getter(&PreviewSystem::xi, py::return_value_policy<py::copy_non_const_reference>()), py::make_setter(&PreviewSystem::xi));

    //AutoSpan
    py::class_<AutoSpan, boost::noncopyable>("AutoSpan", "Helper functions to automatically extend a matrix to a desired dimension.", py::no_init)
        .def("span_matrix", py::make_function(autoSpanMatrix0))
        .def("span_matrix", py::make_function(autoSpanMatrix1))
        .staticmethod("span_matrix")
        .def("span_vector", py::make_function(autoSpanVector))
        .staticmethod("span_vector");

    //Constraint
    py::enum_<ConstraintFlag>("ConstraintFlag", "Flag to constraint type")
        .value("Constraint", ConstraintFlag::Constraint)
        .value("EqualityConstraint", ConstraintFlag::EqualityConstraint)
        .value("InequalityConstraint", ConstraintFlag::InequalityConstraint)
        .value("BoundConstraint", ConstraintFlag::BoundConstraint);

    struct ConstraintWrap : Constraint, py::wrapper<Constraint> {
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

    py::class_<ConstraintWrap, boost::noncopyable>("Constraint", py::no_init) // Disable the constructor because of move semantics. No need anyway.
        .def("auto_span", py::pure_virtual(&Constraint::autoSpan))
        .def("initialize_constraint", py::pure_virtual(&Constraint::initializeConstraint))
        .def("update", py::pure_virtual(&Constraint::update))
        .def("constraint_type", py::pure_virtual(&Constraint::constraintType))
        .def("name", &Constraint::name, py::return_value_policy<py::copy_const_reference>())
        .def("nr_constr", &Constraint::nrConstr);

    struct EqIneqConstraintWrap : EqIneqConstraint, py::wrapper<EqIneqConstraint> {
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

    py::class_<EqIneqConstraintWrap, boost::noncopyable, py::bases<Constraint>>("EqIneqConstraint", py::init<const std::string&, bool>())
        .def("A", &EqIneqConstraint::A, py::return_value_policy<py::copy_const_reference>())
        .def("b", &EqIneqConstraint::b, py::return_value_policy<py::copy_const_reference>());

    // Delete constructor to enforce call of New<Name_of_constraint> function that return a shared_ptr.
    py::class_<TrajectoryConstraint, std::shared_ptr<TrajectoryConstraint>, boost::noncopyable, py::bases<EqIneqConstraint>>("TrajectoryConstraint", "Trajectory constraint. The object is instansiable through a NewTrajectoryConstraint function", py::no_init)
        .def("__init__", py::make_constructor(TrajectoryConstraint1))
        .def("__init__", py::make_constructor(TrajectoryConstraint2));
    py::class_<ControlConstraint, std::shared_ptr<ControlConstraint>, boost::noncopyable, py::bases<EqIneqConstraint>>("ControlConstraint", "Control constraint. The object is instansiable through a NewControlConstraint function", py::no_init)
        .def("__init__", py::make_constructor(ControlConstraint1))
        .def("__init__", py::make_constructor(ControlConstraint2));
    py::class_<MixedConstraint, std::shared_ptr<MixedConstraint>, boost::noncopyable, py::bases<EqIneqConstraint>>("MixedConstraint", "Mixed constraint. The object is instansiable through a NewMixedConstraint function", py::no_init)
        .def("__init__", py::make_constructor(MixedConstraint1))
        .def("__init__", py::make_constructor(MixedConstraint2));
    py::class_<TrajectoryBoundConstraint, std::shared_ptr<TrajectoryBoundConstraint>, boost::noncopyable, py::bases<EqIneqConstraint>>("TrajectoryBoundConstraint", "Trajectory Bound constraint. The object is instansiable through a NewTrajectoryBoundConstraint function", py::no_init)
        .def("__init__", py::make_constructor(TrajectoryBoundConstraint1));
    py::class_<ControlBoundConstraint, std::shared_ptr<ControlBoundConstraint>, boost::noncopyable, py::bases<Constraint>>("ControlBoundConstraint", "Control Bound constraint. The object is instansiable through a NewControlBoundConstraint function", py::no_init)
        .def("__init__", py::make_constructor(ControlBoundConstraint1))
        .def("lower", &ControlBoundConstraint::lower, py::return_value_policy<py::copy_const_reference>())
        .def("upper", &ControlBoundConstraint::upper, py::return_value_policy<py::copy_const_reference>());

    //Cost Functions
    struct CostFunctionWrap : CostFunction, py::wrapper<CostFunction> {
        using CostFunction::CostFunction;

        void autoSpan()
        {
            if (py::override autoSpan = this->get_override("autoSpan"))
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
            if (py::override initializeConstraint = this->get_override("initializeCost"))
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

    py::class_<CostFunctionWrap, boost::noncopyable>("CostFunction", py::no_init) // Disable the constructor because of move semantics. No need anyway.
        .def("initialize_constraint", &CostFunction::initializeCost, &CostFunctionWrap::default_initializeCost)
        .def("update", py::pure_virtual(&CostFunction::update))
        .def("name", &CostFunction::name, py::return_value_policy<py::copy_const_reference>())
        .def("Q", &CostFunction::Q, py::return_value_policy<py::copy_const_reference>())
        .def("c", &CostFunction::c, py::return_value_policy<py::copy_const_reference>())
        .def("weights", &CostFunction::weights<const Eigen::VectorXd&>)
        .def("weight", &CostFunction::weight<double>);

    // Costs
    py::class_<TrajectoryCost, std::shared_ptr<TrajectoryCost>, boost::noncopyable, py::bases<CostFunction>>("TrajectoryCost", "Trajectory cost. The object is instansiable through a NewTrajectoryCost function", py::no_init)
        .def("__init__", py::make_constructor(TrajectoryCost1));
    py::class_<TargetCost, std::shared_ptr<TargetCost>, boost::noncopyable, py::bases<CostFunction>>("TargetCost", "Target cost. The object is instansiable through a NewTargetCost function", py::no_init)
        .def("__init__", py::make_constructor(TargetCost1));
    py::class_<ControlCost, std::shared_ptr<ControlCost>, boost::noncopyable, py::bases<CostFunction>>("ControlCost", "Control cost. The object is instansiable through a NewControlCost function", py::no_init)
        .def("__init__", py::make_constructor(ControlCost1));
    py::class_<MixedCost, std::shared_ptr<MixedCost>, boost::noncopyable, py::bases<CostFunction>>("MixedCost", "Mixed cost. The object is instansiable through a NewMixedCost function", py::no_init)
        .def("__init__", py::make_constructor(MixedCost1));

    // LMPC
    py::class_<LMPC, boost::noncopyable>("LMPC", "LMPC. This class runs the lmpc with the desired QP and fills the PreviewSystem it is attach to", py::init<py::optional<SolverFlag>>())
        .def(py::init<const std::shared_ptr<PreviewSystem>&, py::optional<SolverFlag>>())
        .def("select_qp_Solver", &LMPC::selectQPSolver)
        .def("initialize_controller", &LMPC::initializeController)
        .def("solve", &LMPC::solve)
        .def("solve_time", &LMPC::solveTime)
        .def("solve_and_build_time", &LMPC::solveAndBuildTime)
        .def("control", &LMPC::control, py::return_value_policy<py::copy_const_reference>())
        .def("trajectory", &LMPC::trajectory)
        .def("add_cost", &LMPC::addCost)
        .def("add_constraint", &LMPC::addConstraint)
        .def("reset_constraints", &LMPC::resetConstraints);

    // Implicit conversion of pointers
    py::implicitly_convertible<std::shared_ptr<ControlBoundConstraint>, std::shared_ptr<Constraint>>();
    py::implicitly_convertible<std::shared_ptr<ControlConstraint>, std::shared_ptr<Constraint>>();
    py::implicitly_convertible<std::shared_ptr<ControlCost>, std::shared_ptr<CostFunction>>();
    py::implicitly_convertible<std::shared_ptr<MixedConstraint>, std::shared_ptr<Constraint>>();
    py::implicitly_convertible<std::shared_ptr<MixedCost>, std::shared_ptr<CostFunction>>();
    py::implicitly_convertible<std::shared_ptr<TargetCost>, std::shared_ptr<CostFunction>>();
    py::implicitly_convertible<std::shared_ptr<TrajectoryBoundConstraint>, std::shared_ptr<Constraint>>();
    py::implicitly_convertible<std::shared_ptr<TrajectoryConstraint>, std::shared_ptr<Constraint>>();
    py::implicitly_convertible<std::shared_ptr<TrajectoryCost>, std::shared_ptr<CostFunction>>();
}
