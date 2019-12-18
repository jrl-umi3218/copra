#
# Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

import unittest2
import pyCopra as copra
import numpy as np
from sys import float_info


class TestCopra(unittest2.TestCase):
    def setUp(self):
        self.timestep = 0.005
        self.mass = 5
        self.nbStep = 300
        self.A = np.identity(2)
        self.A[0, 1] = self.timestep
        self.B = np.array([[0.5 * self.timestep * self.timestep / self.mass], [self.timestep / self.mass]])
        self.c = np.array([(-9.81 / 2.) * self.timestep**2, -9.81 * self.timestep])
        self.x0 = np.array([0., -5.])
        self.wu = np.array([1e-4])
        self.wx = np.array([10., 10000.])

        # Costs
        self.xd = np.zeros((2,))
        self.ud = np.zeros((1,))
        self.M = np.identity(2)
        self.N = np.ones((1, 1))

        # Inequality constraints
        self.Gineq = np.ones((1, 1))
        self.hineq = np.array([200.])
        self.Eineq = np.zeros((1, 2))
        self.Eineq[0, 1] = 1
        self.fineq = np.zeros((1,))

        # Bound constraints
        self.uLower = np.zeros((1,))
        self.uUpper = np.zeros((1,))
        self.xLower = np.zeros((2,))
        self.xUpper = np.zeros((2,))
        self.uLower[0] = -float('Inf')
        self.uUpper[0] = 200
        self.xLower[0] = -float('Inf')
        self.xLower[1] = -float('Inf')
        self.xUpper[0] = float('Inf')
        self.xUpper[1] = 0

        # Equality constraints
        self.x0Eq = np.zeros((2,))
        self.xdEq = np.zeros((2,))
        self.Geq = np.ones((1, 1))
        self.heq = np.array([200.])
        self.Eeq = np.zeros((2, 2))
        self.Eeq[0, 0] = 1
        self.feq = self.x0Eq

    def test_lmpc_ineq(self):
        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.nbStep)

        controller = copra.LMPC(ps)
        xCost = copra.TargetCost(self.M, -self.xd)
        uCost = copra.ControlCost(self.N, -self.ud)
        trajConstr = copra.TrajectoryConstraint(self.Eineq, self.fineq)
        contConstr = copra.ControlConstraint(self.Gineq, self.hineq)
        xCost.weights(self.wx)
        uCost.weights(self.wu)

        controller.add_cost(xCost)
        controller.add_cost(uCost)
        controller.add_constraint(trajConstr)
        controller.add_constraint(contConstr)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = int(len(fullTraj) / 2)
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in range(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        self.assertLessEqual(np.amax(control), self.hineq[0])

        print("Test lmpc with inequalities")
        print(controller.solve_time(), "s")
        print(controller.solve_and_build_time(), "s")
        print()

    def test_lmpc_mixed(self):
        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.nbStep)

        controller = copra.LMPC(ps)
        xCost = copra.TargetCost(self.M, -self.xd)
        uCost = copra.ControlCost(self.N, -self.ud)
        mixedConstr = copra.MixedConstraint(self.Eineq, self.Gineq, self.hineq)
        xCost.weights(self.wx)
        uCost.weights(self.wu)

        controller.add_cost(xCost)
        controller.add_cost(uCost)
        controller.add_constraint(mixedConstr)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = int(len(fullTraj) / 2)
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in range(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        for i in range(self.Eineq.shape[0]):
            res = 0
            for j in range(self.Eineq.shape[1]):
                res += self.Eineq[i, j] * fullTraj[i * self.Eineq.shape[1] + j]
            for j in range(self.Gineq.shape[1]):
                res += self.Gineq[i, j] * control[i * self.Gineq.shape[1] + j]
            self.assertLessEqual(res, self.hineq[0])

        print("Test lmpc with inequalities")
        print(controller.solve_time(), "s")
        print(controller.solve_and_build_time(), "s")
        print()

    def test_lmpc_bound(self):
        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.nbStep)

        controller = copra.LMPC(ps)
        xCost = copra.TargetCost(self.M, -self.xd)
        uCost = copra.ControlCost(self.N, -self.ud)
        trajConstr = copra.TrajectoryBoundConstraint(self.xLower, self.xUpper)
        contConstr = copra.ControlBoundConstraint(self.uLower, self.uUpper)
        xCost.weights(self.wx)
        uCost.weights(self.wu)

        controller.add_cost(xCost)
        controller.add_cost(uCost)
        controller.add_constraint(trajConstr)
        controller.add_constraint(contConstr)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = int(len(fullTraj) / 2)
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in range(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        self.assertLessEqual(max(velTraj), self.xUpper[1] + 1e-6)
        self.assertLessEqual(max(control), self.uUpper[0] + 1e-6)

        print("Test lmpc with bounds")
        print(controller.solve_time(), "s")
        print(controller.solve_and_build_time(), "s")
        print()

    def test_lmpc_eq(self):
        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0Eq, self.nbStep)

        controller = copra.LMPC(ps)
        xCost = copra.TargetCost(self.M, -self.xdEq)
        uCost = copra.ControlCost(self.N, -self.ud)
        trajConstr = copra.TrajectoryConstraint(self.Eeq, self.feq, False)
        xCost.weights(self.wx)
        uCost.weights(self.wu)

        controller.add_cost(xCost)
        controller.add_cost(uCost)
        controller.add_constraint(trajConstr)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = int(len(fullTraj) / 2)
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in range(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0] + 1e-6)
        self.assertLessEqual(max(velTraj), self.feq[0] + 1e-6)

        print("Test lmpc with equalities")
        print(controller.solve_time(), "s")
        print(controller.solve_and_build_time(), "s")
        print()

    def test_constructors_initialisations(self):
        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.nbStep)

        controller = copra.LMPC(ps)

        copra.LMPC()
        copra.LMPC(copra.SolverFlag.QuadProgDense)
        copra.LMPC(ps, copra.SolverFlag.QuadProgDense)

        controller.initialize_controller(ps)

    @unittest2.expectedFailure
    def test_fail_construct_trajectory(self):
        constr = copra.TrajectoryConstraint()

    @unittest2.expectedFailure
    def test_fail_construct_control(self):
        constr = copra.ControlConstraint()

    @unittest2.expectedFailure
    def test_fail_construct_control(self):
        constr = copra.TrajectoryBoundConstraint()

    @unittest2.expectedFailure
    def test_fail_construct_control(self):
        constr = copra.ControlBoundConstraint()

    def test_constraint_and_cost_deletion(self):
        print("Testing 'test_constraint_deletion'.")
        print("In order to see the outputs, the copra must be installed under Debug mode.")

        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.nbStep)

        controller = copra.LMPC(ps)
        trajConstr = copra.TrajectoryConstraint(self.Eineq, self.fineq)
        contConstr = copra.ControlConstraint(self.Gineq, self.hineq)
        trajEqConstr = copra.TrajectoryConstraint(self.Eeq, self.feq, False)
        contEqConstr = copra.ControlConstraint(self.Geq, self.heq, False)
        trajBdConstr = copra.TrajectoryBoundConstraint(self.xLower, self.xUpper)
        contBdConstr = copra.ControlBoundConstraint(self.uLower, self.uUpper)
        targetCost = copra.TargetCost(self.M, -self.xd)
        trajectoryCost = copra.TrajectoryCost(self.M, -self.xd)
        controlCost = copra.ControlCost(self.N, -self.ud)
        M_mixed = np.ones((1, 2))
        mixedCost = copra.MixedCost(M_mixed, self.N, -self.ud)
        
        controller.add_constraint(trajConstr)
        controller.add_constraint(contConstr)
        controller.add_constraint(trajEqConstr)
        controller.add_constraint(contEqConstr)
        controller.add_constraint(trajBdConstr)
        controller.add_constraint(contBdConstr)
        controller.add_cost(targetCost)
        controller.add_cost(trajectoryCost)
        controller.add_cost(controlCost)
        controller.add_cost(mixedCost)

        del trajConstr

        targetCost.weights(self.wx)
        controlCost.weights(self.wu)

        del trajEqConstr
        del contEqConstr
        del trajBdConstr
        del contBdConstr
        del trajectoryCost
        del mixedCost

        self.assertFalse(controller.solve())
        self.assertTrue(controller.solve()) # Has kept the contConstr only

    def test_preview_systeme_still_exist(self):
        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.nbStep)

        controller = copra.LMPC(ps)
        del ps
        trajConstr = copra.TrajectoryConstraint(self.Eineq, self.fineq)
        contConstr = copra.ControlConstraint(self.Gineq, self.hineq)
        targetCost = copra.TargetCost(self.M, -self.xd)
        controlCost = copra.ControlCost(self.N, -self.ud)
        targetCost.weights(self.wx)
        controlCost.weights(self.wu)

        controller.add_constraint(trajConstr)
        controller.add_constraint(contConstr)
        controller.add_cost(targetCost)
        controller.add_cost(controlCost)

        self.assertTrue(controller.solve())
        control = controller.control()
        fullTraj = controller.trajectory()
        fTLen = int(len(fullTraj) / 2)
        posTraj = [0.] * fTLen
        velTraj = [0.] * fTLen
        for i in range(fTLen):
            posTraj[i] = fullTraj[2 * i]
            velTraj[i] = fullTraj[2 * i + 1]

        self.assertAlmostEqual(self.xd[1], velTraj[-1], places=3)
        self.assertLessEqual(max(posTraj), self.x0[0])
        self.assertLessEqual(np.amax(control), self.hineq[0])

    def test_throw_handler(self):
        ps = copra.PreviewSystem()
        ps.system(self.A, self.B, self.c, self.x0, self.nbStep)

        controller = copra.LMPC(ps)
        # Test trajectory constraint throws
        with self.assertRaises(RuntimeError):
            constr = copra.TrajectoryConstraint(np.identity(5), np.ones((2,)))
            controller.add_constraint(constr)

        # Test control constraint throws
        with self.assertRaises(RuntimeError):
            constr = copra.ControlConstraint(np.identity(5), np.ones((2,)))
            controller.add_constraint(constr)

        # Test mixed constraint throws
        with self.assertRaises(RuntimeError):
            constr = copra.MixedConstraint(np.identity(5), np.identity(5), np.ones((2,)))
            controller.add_constraint(constr)

        # Test trajectory bound constraint throws
        with self.assertRaises(RuntimeError):
            constr = copra.TrajectoryBoundConstraint(np.ones((3,)), np.ones((2,)))
            controller.add_constraint(constr)

        # Test control bound constraint throws
        with self.assertRaises(RuntimeError):
            constr = copra.ControlBoundConstraint(np.ones((3,)), np.ones((2,)))
            controller.add_constraint(constr)

    def test_dynamic_walk(self):
        A = np.array([
            [      1,      0,      0,0.11699999999999999,      0,      0],
            [      0,      1,      0,      0,0.11699999999999999,      0],
            [      0,      0,      1,      0,      0,0.11699999999999999],
            [      0,      0,      0,      1,      0,      0],
            [      0,      0,      0,      0,      1,      0],
            [      0,      0,      0,      0,      0,      1]])
        B = np.array([
	        [0.006844499999999999,      0,      0],
            [      0,0.006844499999999999,      0],
            [      0,      0,0.006844499999999999],
            [0.11699999999999999,      0,      0],
            [      0,0.11699999999999999,      0],
            [      0,      0,0.11699999999999999]])

        c = np.array([0., 0., 0., 0., 0., 0.])
        x_init = np.array([1.5842778860957882, 0.3422260214935311, 2.289067474385933, 0., 0., 0.])
        x_goal = np.array([1.627772868473883, 0.4156386515475985, 2.3984423755527136, 0.06745225960685897, 0.3882830795737303, 0.06845759848745198])
        nb_steps = 10
        G = np.array([
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

        h = np.array([47.76613348420254, 9.80665, 9.80665, 9.80665, 9.806649999999998, 49.86555936465245, 9.806649999999994, 47.76613348420254, 9.80665, 9.80665, 9.80665, 9.806649999999998, 49.86555936465245, 9.806649999999994, 47.76613348420254, 9.80665, 9.80665, 9.80665, 9.806649999999998, 49.86555936465245, 9.806649999999994, 47.76613348420254, 9.80665, 9.80665, 9.80665, 9.806649999999998, 49.86555936465245, 9.806649999999994, 47.76613348420254, 9.80665, 9.80665, 9.80665, 9.806649999999998, 49.86555936465245, 9.806649999999994, 47.76613348420254, 9.80665, 9.80665, 9.80665, 9.806649999999998, 49.86555936465245, 9.806649999999994, 9.806650000000007, 9.806650000000008, 9.80665, 9.806650000000001, 9.806650000000007, 9.806649999999996, 9.806650000000007, 9.806650000000008, 9.80665, 9.806650000000001, 9.806650000000007, 9.806649999999996, 9.806650000000007, 9.806650000000008, 9.80665, 9.806650000000001, 9.806650000000007, 9.806649999999996, 9.806650000000007, 9.806650000000008, 9.80665, 9.806650000000001, 9.806650000000007, 9.806649999999996])

        ps = copra.PreviewSystem()
        ps.system(A, B, c, x_init, nb_steps)

        controller = copra.LMPC(ps)
        contConstr = copra.ControlConstraint(G, h)
        M_cost = np.identity(6)
        targetCost = copra.TargetCost(M_cost, -x_goal)

        controller.add_constraint(contConstr)
        controller.add_cost(targetCost)

        controller.solve()
        print(controller.solve_time(), "s")

nb_steps = 10

if __name__ == '__main__':
    unittest2.main()
