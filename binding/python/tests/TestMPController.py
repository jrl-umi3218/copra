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

        controller.selectQPSolver(mpc.SolverFlag.LSSOL)
        controller.updateSystem(ps)
        self.assertTrue(controller.solve(ps))
        fullTraj = controller.trajectory(ps)
        print fullTraj[599]

if __name__ == '__main__':
    unittest.main()