import unittest
import MPController as mpc
from minieigen import *

class TestMPC(unittest.TestCase):
    def setUp(self):
        self.timestep = 0.005
        self.mass = 5
        self.nbStep = 300
        self.A = MatrixXd.Identity(2, 2)
        self.A[0,1] = self.timestep
        self.B = Vector2d(0.5 * self.timestep * self.timestep / self.mass, self.timestep / self.mass)
        self.c = Vector2d(0; -9.81 * self.timestep)
        self.G = MatrixXd.Identity(1, 1)
        self.h = VectorXd.ones(1) * 200.
        self.E = MatrixXd.Zero(1, 2)
        self.E[0,1] = 1
        self.f = VectorXd.zero(1)
        self.x0 = Vector2d(0, -5)
        self.xd = Vector2d.Zero()
        self.wx = Vector2d(10, 10000)
        self.wu = VectorXd.Ones(1) * 1e-4

    def test_mpcLast(self):
        


if __name__ == '__main__':
    unittest.main()