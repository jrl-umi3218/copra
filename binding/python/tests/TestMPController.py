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
        trajConstr = mpc.NewTrajectoryConstraint(self.Eineq, self.fineq)
        contConstr = mpc.NewControlConstraint(self.Gineq, self.hineq)

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
        trajConstr = mpc.NewTrajectoryConstraint(self.Eineq, self.fineq)
        contConstr = mpc.NewControlConstraint(self.Gineq, self.hineq)
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

        controller.weights(self.wx, self.wu)

        del trajEqConstr
        del contEqConstr
        del trajBdConstr
        del contBdConstr

        self.assertFalse(controller.solve())
        self.assertTrue(controller.solve()) # Has kept the contConstr

    def test_preview_systeme_still_exist(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        del ps
        trajConstr = mpc.NewTrajectoryConstraint(self.Eineq, self.fineq)
        contConstr = mpc.NewControlConstraint(self.Gineq, self.hineq)

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

    def test_throw_handler(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        # Test trajectory constraint throws
        with self.assertRaises(RuntimeError):
            mpc.NewTrajectoryConstraint(MatrixXd.Identity(5, 5), VectorXd.Ones(2))
        
        trajConstr = mpc.NewTrajectoryConstraint(self.Eineq, self.fineq)
        with self.assertRaises(RuntimeError):
            trajConstr.trajectory(MatrixXd.Identity(5, 5), self.fineq)
        with self.assertRaises(RuntimeError):
            trajConstr.trajectory(self.Eineq, VectorXd.Ones(2))
        
        trajConstr = mpc.NewTrajectoryConstraint(MatrixXd.Identity(5, 5), VectorXd.Ones(5))
        with self.assertRaises(RuntimeError):
            controller.addConstraint(trajConstr)

        # Test control constraint throws
        with self.assertRaises(RuntimeError):
            mpc.NewControlConstraint(MatrixXd.Identity(5, 5), VectorXd.Ones(2))
        
        contConstr = mpc.NewControlConstraint(self.Gineq, self.hineq)
        with self.assertRaises(RuntimeError):
            contConstr.control(MatrixXd.Identity(5, 5), self.hineq)
        with self.assertRaises(RuntimeError):
            contConstr.control(self.Gineq, VectorXd.Ones(2))
        
        contConstr = mpc.NewControlConstraint(MatrixXd.Identity(5, 5), VectorXd.Ones(5))
        with self.assertRaises(RuntimeError):
            controller.addConstraint(contConstr)

        # Test trajectory bound constraint throws
        with self.assertRaises(RuntimeError):
            mpc.NewTrajectoryBoundConstraint(VectorXd.Ones(3), VectorXd.Ones(2))
        
        tbConstr = mpc.NewTrajectoryBoundConstraint(self.xLower, self.xUpper)
        with self.assertRaises(RuntimeError):
            tbConstr.trajectoryBound(VectorXd.Ones(3), self.xUpper)
        with self.assertRaises(RuntimeError):
            tbConstr.trajectoryBound(self.xLower, VectorXd.Ones(3))
        
        tbConstr = mpc.NewTrajectoryBoundConstraint(VectorXd.Ones(3), VectorXd.Ones(3))
        with self.assertRaises(RuntimeError):
            controller.addConstraint(tbConstr)

        # Test control bound constraint throws
        with self.assertRaises(RuntimeError):
            mpc.NewControlBoundConstraint(VectorXd.Ones(3), VectorXd.Ones(2))
        
        cbConstr = mpc.NewControlBoundConstraint(self.uLower, self.uUpper)
        with self.assertRaises(RuntimeError):
            cbConstr.controlBound(VectorXd.Ones(3), self.uUpper)
        with self.assertRaises(RuntimeError):
            cbConstr.controlBound(self.uLower, VectorXd.Ones(3))
        
        cbConstr = mpc.NewControlBoundConstraint(VectorXd.Ones(3), VectorXd.Ones(3))
        with self.assertRaises(RuntimeError):
            controller.addConstraint(cbConstr)

    def test_weights(self):
        ps = mpc.NewPreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, 2)

        controllerF = mpc.MPCTypeFull(ps)
        controllerL = mpc.MPCTypeLast(ps)

        controllerF.weights(2.3, 2.6)
        controllerL.weights(2.3, 2.6)
        controllerF.weights(self.wx, self.wu)
        controllerL.weights(self.wx, self.wu)
        wx = VectorXd([1, 2, 3, 4])
        wu = VectorXd([3, 4])
        controllerF.weights(wx, wu)
        wx = VectorXd([1, 2])
        controllerL.weights(wx, wu)
        wx = VectorXd([3, 4, 5])
        with self.assertRaises(RuntimeError):
            controllerF.weights(wx, wu)
        with self.assertRaises(RuntimeError):
            controllerL.weights(wx, wu)
        wx = VectorXd([1, 2])
        wu = VectorXd([3, 4, 7])
        with self.assertRaises(RuntimeError):
            controllerF.weights(wx, wu)
        with self.assertRaises(RuntimeError):
            controllerL.weights(wx, wu)
        # Just to check that a message is printed
        try:
            controllerL.weights(wx, wu)
            self.assertTrue(False)
        except RuntimeError as e:
            print 'Test error message handler'
            print e
            self.assertTrue(True)
            

    def test_dynamic_walk(self):
        A = MatrixXd([
            [      1,      0,      0,0.11699999999999999,      0,      0],
            [      0,      1,      0,      0,0.11699999999999999,      0],
            [      0,      0,      1,      0,      0,0.11699999999999999],
            [      0,      0,      0,      1,      0,      0],
            [      0,      0,      0,      0,      1,      0],
            [      0,      0,      0,      0,      0,      1]])
        B = MatrixXd([
	        [0.006844499999999999,      0,      0],
            [      0,0.006844499999999999,      0],
            [      0,      0,0.006844499999999999],
            [0.11699999999999999,      0,      0],
            [      0,0.11699999999999999,      0],
            [      0,      0,0.11699999999999999]])

        c = VectorXd([0,0,0, 0,0,0])
        x_init = VectorXd([1.5842778860957882,0.3422260214935311,2.289067474385933, 0,0,0])
        x_goal = VectorXd([1.627772868473883,0.4156386515475985,2.3984423755527136, 0.06745225960685897,0.3882830795737303,0.06845759848745198])
        nb_steps = 10
        G = MatrixXd([
            [     -1,9.946646523934742,-4.870790074510924,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [-18.826459196882055,3.4468275392859393,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [-1.374181960437557,-8.028252906078723,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [-9.936224732113594,5.000580301294253,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [8.750597187695343,-4.7538557382857105,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [18.65430148414319,      1,-5.084871935334947,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [2.1775137248880574e-15,-5.443784312220143e-16,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,     -1,9.946646523934742,-4.870790074510924,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,-18.826459196882055,3.4468275392859393,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,-1.374181960437557,-8.028252906078723,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,-9.936224732113594,5.000580301294253,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,8.750597187695343,-4.7538557382857105,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,18.65430148414319,      1,-5.084871935334947,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,2.1775137248880574e-15,-5.443784312220143e-16,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,     -1,9.946646523934742,-4.870790074510924,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,-18.826459196882055,3.4468275392859393,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,-1.374181960437557,-8.028252906078723,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,-9.936224732113594,5.000580301294253,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,8.750597187695343,-4.7538557382857105,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,18.65430148414319,      1,-5.084871935334947,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,2.1775137248880574e-15,-5.443784312220143e-16,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,     -1,9.946646523934742,-4.870790074510924,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,-18.826459196882055,3.4468275392859393,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,-1.374181960437557,-8.028252906078723,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,-9.936224732113594,5.000580301294253,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,8.750597187695343,-4.7538557382857105,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,18.65430148414319,      1,-5.084871935334947,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,2.1775137248880574e-15,-5.443784312220143e-16,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,     -1,9.946646523934742,-4.870790074510924,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-18.826459196882055,3.4468275392859393,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-1.374181960437557,-8.028252906078723,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-9.936224732113594,5.000580301294253,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,8.750597187695343,-4.7538557382857105,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,18.65430148414319,      1,-5.084871935334947,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,2.1775137248880574e-15,-5.443784312220143e-16,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,     -1,9.946646523934742,-4.870790074510924,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-18.826459196882055,3.4468275392859393,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-1.374181960437557,-8.028252906078723,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-9.936224732113594,5.000580301294253,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,8.750597187695343,-4.7538557382857105,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,18.65430148414319,      1,-5.084871935334947,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,2.1775137248880574e-15,-5.443784312220143e-16,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,8.750597218241072,-4.753855754313641,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-1.3741819771739483,-8.028252929943818,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-18.82645925631406,3.4468275193254927,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,1.1824397341134247,7.4638136184143935,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,14.006645137157978,-2.2569159229140494,     -1,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,2.1775137248880578e-15,-5.443784312220144e-16,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,8.750597218241072,-4.753855754313641,     -1,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-1.3741819771739483,-8.028252929943818,     -1,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-18.82645925631406,3.4468275193254927,     -1,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,1.1824397341134247,7.4638136184143935,     -1,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,14.006645137157978,-2.2569159229140494,     -1,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,2.1775137248880578e-15,-5.443784312220144e-16,      1,      0,      0,      0,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,8.750597218241072,-4.753855754313641,     -1,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-1.3741819771739483,-8.028252929943818,     -1,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-18.82645925631406,3.4468275193254927,     -1,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,1.1824397341134247,7.4638136184143935,     -1,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,14.006645137157978,-2.2569159229140494,     -1,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,2.1775137248880578e-15,-5.443784312220144e-16,      1,      0,      0,      0],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,8.750597218241072,-4.753855754313641,     -1],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-1.3741819771739483,-8.028252929943818,     -1],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,-18.82645925631406,3.4468275193254927,     -1],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,1.1824397341134247,7.4638136184143935,     -1],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,14.006645137157978,-2.2569159229140494,     -1],
            [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,2.1775137248880578e-15,-5.443784312220144e-16,      1]])

        h = VectorXd([47.76613348420254,9.80665,9.80665, 9.80665,9.806649999999998,49.86555936465245, 9.806649999999994,47.76613348420254,9.80665, 9.80665,9.80665,9.806649999999998, 49.86555936465245,9.806649999999994,47.76613348420254, 9.80665,9.80665,9.80665, 9.806649999999998,49.86555936465245,9.806649999999994, 47.76613348420254,9.80665,9.80665, 9.80665,9.806649999999998,49.86555936465245, 9.806649999999994,47.76613348420254,9.80665, 9.80665,9.80665,9.806649999999998, 49.86555936465245,9.806649999999994,47.76613348420254, 9.80665,9.80665,9.80665, 9.806649999999998,49.86555936465245,9.806649999999994, 9.806650000000007,9.806650000000008,9.80665, 9.806650000000001,9.806650000000007,9.806649999999996, 9.806650000000007,9.806650000000008,9.80665, 9.806650000000001,9.806650000000007,9.806649999999996, 9.806650000000007,9.806650000000008,9.80665, 9.806650000000001,9.806650000000007,9.806649999999996, 9.806650000000007,9.806650000000008,9.80665, 9.806650000000001,9.806650000000007,9.806649999999996])

        ps = mpc.NewPreviewSystem()
        ps.system(A, B, c, x_init, x_goal, nb_steps)

        controller = mpc.MPCTypeLast(ps)
        contConstr = mpc.NewControlConstraint(G, h)

        controller.addConstraint(contConstr)

        controller.solve()
        print controller.solveTime().wall*1e-6

nb_steps = 10

if __name__ == '__main__':
    unittest.main()
