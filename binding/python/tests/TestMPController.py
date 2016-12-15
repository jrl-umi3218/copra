import unittest
import mpcontroller as mpc
from minieigen import *
from sys import float_info


class TestMPC(unittest.TestCase):

    def setUp(self):
        self.timestep = 0.005
        self.mass = 5
        self.nbStep = 300
        self.A = MatrixXd.Identity(2, 2)
        self.A[0, 1] = self.timestep
        self.B = Vector2d(0.5 * self.timestep * self.timestep /
                          self.mass, self.timestep / self.mass)
        self.c = Vector2d(0, -9.81 * self.timestep)
        self.x0 = Vector2d(0, -5)
        self.xd = Vector2d.Zero
        self.wu = VectorXd.Ones(1) * 1e-4
        self.wx = Vector2d(10, 10000)

        # Inequality constraints
        self.Gineq = MatrixXd.Ones(1, 1)
        self.hineq = VectorXd.Ones(1) * 200.
        self.Eineq = MatrixXd.Zero(1, 2)
        self.Eineq[0, 1] = 1
        self.fineq = VectorXd.Zero(1)

        # Bound constraints
        self.uLower = VectorXd.Zero(1)
        self.uUpper = VectorXd.Zero(1)
        self.xLower = VectorXd.Zero(2)
        self.xUpper = VectorXd.Zero(2)
        self.uLower[0] = -float('Inf')
        self.uUpper[0] = 200
        self.xLower[0] = -float('Inf')
        self.xLower[1] = -float('Inf')
        self.xUpper[0] = float('Inf')
        self.xUpper[1] = 0

        # Equality constraints
        self.x0Eq = Vector2d(0, 0)
        self.xdEq = Vector2d(0, 0)
        self.Geq = MatrixXd.Ones(1, 1)
        self.heq = VectorXd.Ones(1) * 200.
        self.Eeq = MatrixXd.Zero(2, 2)
        self.Eeq[0, 0] = 1
        self.feq = self.x0Eq

    def test_mpcLast_ineq(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        trajConstr = mpc.NewTrajectoryConstraint(self.Eineq, self.fineq, True)
        contConstr = mpc.NewControlConstraint(self.Gineq, self.hineq, True)

        controller.addConstraint(trajConstr)
        controller.addConstraint(contConstr)

        controller.weights(self.wx, self.wu)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = len(fullTraj) / 2
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in xrange(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        self.assertLessEqual(control.maxCoeff(), self.hineq[0])

        print "Test mpc last with inequalities"
        print controller.solveTime().wall * 1e-6
        print controller.solveAndBuildTime().wall * 1e-6
        print

    def test_mpcLast_bound(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        trajConstr = mpc.NewTrajectoryBoundConstraint(self.xLower, self.xUpper)
        contConstr = mpc.NewControlBoundConstraint(self.uLower, self.uUpper)

        controller.addConstraint(trajConstr)
        controller.addConstraint(contConstr)

        controller.weights(self.wx, self.wu)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = len(fullTraj) / 2
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in xrange(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        self.assertLessEqual(max(velTraj), self.xUpper[1] + 1e-6)
        self.assertLessEqual(max(control), self.uUpper[0] + 1e-6)

        print "Test mpc last with bounds"
        print controller.solveTime().wall * 1e-6
        print controller.solveAndBuildTime().wall * 1e-6
        print

    def test_mpcLast_eq(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0Eq, self.xdEq, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        trajConstr = mpc.NewTrajectoryConstraint(self.Eeq, self.feq, False)

        controller.addConstraint(trajConstr)

        controller.weights(self.wx, self.wu)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = len(fullTraj) / 2
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in xrange(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0] + 1e-6)
        self.assertLessEqual(max(velTraj), self.feq[0] + 1e-6)

        print "Test mpc last with equalities"
        print controller.solveTime().wall * 1e-6
        print controller.solveAndBuildTime().wall * 1e-6
        print

    def test_constructors_initialisations(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        c1 = mpc.MPCTypeFull(ps)
        mpc.MPCTypeFull(mpc.SolverFlag.QuadProgDense)
        mpc.MPCTypeFull(ps)
        mpc.MPCTypeFull(ps, mpc.SolverFlag.QuadProgDense)
        c2 = mpc.MPCTypeLast(ps)
        mpc.MPCTypeLast(mpc.SolverFlag.QuadProgDense)
        mpc.MPCTypeLast(ps)
        mpc.MPCTypeLast(ps, mpc.SolverFlag.QuadProgDense)

        c1.initializeController(ps)
        c2.initializeController(ps)

    @unittest.expectedFailure
    def test_fail_construct_previewsystem(self):
        ps = mpc.PreviewSystem()

    @unittest.expectedFailure
    def test_fail_construct_trajectory(self):
        ps = mpc.TrajectoryConstraint()

    @unittest.expectedFailure
    def test_fail_construct_control(self):
        ps = mpc.ControlConstraint()

    @unittest.expectedFailure
    def test_fail_construct_control(self):
        ps = mpc.TrajectoryBoundConstraint()

    @unittest.expectedFailure
    def test_fail_construct_control(self):
        ps = mpc.ControlBoundConstraint()

    def test_constraint_dangling_pointer(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        trajConstr = mpc.NewTrajectoryConstraint(self.Eineq, self.fineq, True)
        contConstr = mpc.NewControlConstraint(self.Gineq, self.hineq, True)
        trajEqConstr = mpc.NewTrajectoryConstraint(self.Eeq, self.feq, False)
        contEqConstr = mpc.NewControlConstraint(self.Geq, self.heq, False)
        trajBdConstr = mpc.NewTrajectoryBoundConstraint(self.xLower, self.xUpper)
        contBdConstr = mpc.NewControlBoundConstraint(self.uLower, self.uUpper)
        
        controller.addConstraint(trajConstr)
        controller.addConstraint(contConstr)
        controller.addConstraint(trajEqConstr)
        controller.addConstraint(contEqConstr)
        controller.addConstraint(trajBdConstr)
        controller.addConstraint(contBdConstr)

        del trajConstr
        del contConstr

        controller.weights(self.wx, self.wu)

        del trajEqConstr
        del contEqConstr
        del trajBdConstr
        del contBdConstr

        self.assertTrue(controller.solve())

    def test_preview_systeme_still_exist(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        del ps
        trajConstr = mpc.NewTrajectoryConstraint(self.Eineq, self.fineq, True)
        contConstr = mpc.NewControlConstraint(self.Gineq, self.hineq, True)

        controller.addConstraint(trajConstr)
        controller.addConstraint(contConstr)

        controller.weights(self.wx, self.wu)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = len(fullTraj) / 2
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in xrange(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        self.assertLessEqual(control.maxCoeff(), self.hineq[0])


if __name__ == '__main__':
    unittest.main()
