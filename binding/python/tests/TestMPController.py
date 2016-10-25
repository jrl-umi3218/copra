import unittest
import mpcontroller as mpc
from minieigen import *

class TestMPC(unittest.TestCase):
    def setUp(self):
        self.timestep = 0.005
        self.mass = 5
        self.nbStep = 300
        self.A = MatrixXd.Identity(2, 2)
        self.A[0,1] = self.timestep
        self.B = Vector2d(0.5 * self.timestep * self.timestep / self.mass, self.timestep / self.mass)
        self.c = Vector2d(0, -9.81 * self.timestep)
        self.G = MatrixXd.Ones(1, 1)
        self.h = VectorXd.Ones(1) * 200.
        self.E = MatrixXd.Zero(1, 2)
        self.E[0,1] = 1
        self.f = VectorXd.Zero(1)
        self.x0 = Vector2d(0, -5)
        self.xd = Vector2d.Zero
        self.wx = Vector2d(10, 10000)
        self.wu = VectorXd.Ones(1) * 1e-4

    def test_mpcLast(self):
        ps = mpc.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        controller = mpc.MPCTypeLast(ps)
        trajConstr = mpc.TrajectoryConstrain(self.E, self.f)
        contConstr = mpc.ControlConstrain(self.G, self.h)

        controller.addConstrain(ps, trajConstr)
        controller.addConstrain(ps, contConstr)

        controller.weights(ps, self.wx, self.wu)

        controller.updateSystem(ps)
        
        self.assertTrue(controller.solve(ps))
        control = controller.control()
        fullTraj = controller.trajectory(ps)
        fTLen = len(fullTraj)/2
        posTraj = [0.]*fTLen
        velTraj = [0.]*fTLen
        for i in xrange(fTLen):
            posTraj[i] = fullTraj[2*i]
            velTraj[i] = fullTraj[2*i+1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        self.assertLessEqual(control.maxCoeff(), self.h[0])

        print
        print controller.solveTime().wall*1e-6
        print controller.solveAndBuildTime().wall*1e-6
        print

    def test_constructors_initialisations(self):
        ps = mpc.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.xd, self.nbStep)

        c1 = mpc.MPCTypeFull()
        mpc.MPCTypeFull(mpc.SolverFlag.LSSOL)
        mpc.MPCTypeFull(ps)
        mpc.MPCTypeFull(ps, mpc.SolverFlag.LSSOL)
        c2 = mpc.MPCTypeLast()
        mpc.MPCTypeLast(mpc.SolverFlag.LSSOL)
        mpc.MPCTypeLast(ps)
        mpc.MPCTypeLast(ps, mpc.SolverFlag.LSSOL)

        c1.initializeController(ps)
        c2.initializeController(ps)

if __name__ == '__main__':
    unittest.main()