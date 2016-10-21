#include "solverUtils.h"
#include "Constrains.h"
#include "PreviewController.h"
#include "PreviewSystem.h"
#include <boost/python.hpp>

BOOST_PYTHON_MODULE(_mpcontroller)
{
    using namespace mpc;
    using namespace boost::python;

    enum_<SolverFlag>("SolverFlag")
        .value("DEFAULT", SolverFlag::DEFAULT)
#ifdef LSSOL_SOLVER_FOUND
        .value("LSSOL", SolverFlag::LSSOL)
#endif
        .value("QLD", SolverFlag::QLD)
        .value("QuadProgDense", SolverFlag::QuadProgDense);

    // Preview System
    void (PreviewSystem::*sys1)(const Eigen::MatrixXd &, const Eigen::MatrixXd &,
                                const Eigen::VectorXd &, const Eigen::VectorXd &, int) = &PreviewSystem::system;
    void (PreviewSystem::*sys2)(const Eigen::MatrixXd &, const Eigen::MatrixXd &,
                                const Eigen::VectorXd &, const Eigen::VectorXd &, const Eigen::VectorXd &, int) = &PreviewSystem::system;

    class_<PreviewSystem>("PreviewSystem")
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

    //Constrain
    struct ConstrainWrap : Constrain, wrapper<Constrain>
    {
        using Constrain::Constrain;

        void initializeConstrain(const PreviewSystem &ps)
        {
            this->get_override("initializeConstrain")(ps);
        }

        void update(const PreviewSystem &ps)
        {
            this->get_override("update")(ps);
        }
    };

    class_<ConstrainWrap, boost::noncopyable>("Constrain", init<const Eigen::MatrixXd &, const Eigen::VectorXd &>())
        .def("initializeConstrain", pure_virtual(&Constrain::initializeConstrain))
        .def("update", pure_virtual(&Constrain::update))
        .def("nrConstr", &Constrain::nrConstr)
        .def("A", &Constrain::A, return_internal_reference<>())
        .def("b", &Constrain::b, return_internal_reference<>());

    class_<TrajectoryConstrain, boost::noncopyable, bases<ConstrainWrap>>("TrajectoryConstrain",
                                                                          init<const Eigen::MatrixXd &, const Eigen::VectorXd &>());

    class_<TrajectoryConstrain, boost::noncopyable, bases<ConstrainWrap>>("ControlConstrain",
                                                                          init<const Eigen::MatrixXd &, const Eigen::VectorXd &>());

    //MPCTypeFull
    struct MPCTypeFullWrap : MPCTypeFull, wrapper<MPCTypeFull>
    {
        using MPCTypeFull::MPCTypeFull;

        void weights(const PreviewSystem &ps, const Eigen::VectorXd &wx, const Eigen::VectorXd &wu)
        {
            if (override weights = this->get_override("weights"))
                weights(ps, wx, wu);
            MPCTypeFull::weights(ps, wx, wu);
        }

        void default_weights(const PreviewSystem &ps, const Eigen::VectorXd &wx, const Eigen::VectorXd &wu)
        {
            this->MPCTypeFull::weights(ps, wx, wu);
        }
    };

    // The default copy-ctor is implicitely deleted due to ctor overloading
    class_<MPCTypeFullWrap, boost::noncopyable>("MPCTypeFull",
                                                init<const PreviewSystem &, SolverFlag>())
        .def("selectQPSolver", &MPCTypeFull::selectQPSolver)
        .def("updateSystem", &MPCTypeFull::updateSystem) //Need an internal ref
        .def("solve", &MPCTypeFull::solve)
        .def("weights", &MPCTypeFull::weights, &MPCTypeFullWrap::default_weights)
        .def("control", &MPCTypeFull::control, return_internal_reference<>())
        .def("trajectory", &MPCTypeFull::trajectory)
        .def("addConstrain", &MPCTypeFull::addConstrain)
        .def("resetConstrains", &MPCTypeFull::resetConstrains);

    //MPCTypeLast
    class_<MPCTypeLast, boost::noncopyable, bases<MPCTypeFullWrap>>("MPCTypeFull",
                                                                    init<const PreviewSystem &, SolverFlag>());
}