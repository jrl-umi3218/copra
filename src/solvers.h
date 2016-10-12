#pragma once

#include <memory>
#include <Eigen/Core>

// QP solvers
#include <eigen-qld/QLD.h>
#include <eigen-QuadProg/QuadProg.h>
// #include <eigen-lssol/LSSOL.h>

namespace pc
{

/**
 * An interface to the quadratic solvers.
 * This interface is more like a pseudo-interface (it is not an abstract class).
 * This class allows to have a base class for all solvers. 
 * It provides all the necessary functions for using a qp solver.
 * In case of a QP does not have a corresponding function it sends a warning.
 */
class SolverInterface //TODO: Add warning for all functions
{
  public:
    /**
	 * Get the iteration information output as define by the solver documentation.
	 * @see QuadProgDenseSolver::SI_iter()
	 * @see QuadProgSparseSolver::SI_iter()
	 * @see QLDSolver::SI_iter()
	 * @return The iterations result.
	 */
    virtual const Eigen::VectorXi &SI_iter() const {};

    /**
	 * Get information of eventual fail's solver output as define by the solver documentation.
	 * @see QuadProgDenseSolver::SI_fail()
	 * @see QuadProgSparseSolver::SI_fail()
	 * @see QLDSolver::SI_fail()
	 * @return The fail number.
	 */
    virtual int SI_fail() const {};

    /**
	 * Get the solver's solution.
	 * @see QuadProgDenseSolver::SI_result()
	 * @see QuadProgSparseSolver::SI_result()
	 * @see QLDSolver::SI_result()
	 * @return The qp solver result.
	 */
    virtual const Eigen::VectorXd &SI_result() const {};

    /**
	 * Initialize the variables of the problem to solve.
	 * @see QuadProgDenseSolver::SI_problem()
	 * @see QuadProgSparseSolver::SI_problem()
	 * @see QLDSolver::SI_problem()
	 * @return The qp solver result.
	 */
    virtual void SI_problem(int nrVar, int nrEq, int nrInEq){};

    /**
	 * Solve the problem.
	 * @see QuadProgDenseSolver::SI_solve()
	 * @see QuadProgSparseSolver::SI_solve()
	 * @see QLDSolver::SI_solve()
	 * @return The qp solver result.
	 */
    virtual bool SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
			  const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
			  const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq){};
};

/**
 * QuadProg solver for dense matrix.
 */
class QuadProgDenseSolver : public SolverInterface
{
  public:
    QuadProgDenseSolver();
    QuadProgDenseSolver(int nrVar, int nrEq, int nrInEq);

    const Eigen::VectorXi &SI_iter() const override;
    int SI_fail() const override;

    const Eigen::VectorXd &SI_result() const override;

    /**
	 * Initialize the variables of the problem to solve.
	 * @param nrVar The number of decision variables
	 * @param nrEq The number of equality constrains
	 * @param nrInEq The number of inequality constrains
	 * @return The qp solver result.
	 */
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;

    /**
	 * Solve the problem.
	 * Solve the system 
	 * \f[
	 * 	 \left\{
		 	\begin{array}{cl}
			 {
				 \min_x & \frac{1}{2} x^T Q x + c^T x \\
				        & Aeq x = beq \\
						& AInEq x \leq bInEq
			 }
		 \right.
	 * \f]
	 * @param Q An N-by-N symmetric positive definite dense matrix.
	 * @param C An N-by-1 dense vector.
	 * @param Aeq An M-by-N dense matrix.
	 * @param beq An M-by-1 dense vector.
	 * @param Aineq An P-by-N dense matrix.
	 * @param Bineq An P-by-1 dense vector.
	 * @return The qp solver result.
	 */
    bool SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
		  const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
		  const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq) override;

  private:
    std::unique_ptr<QuadProgDense> solver_;
};

/**
 * QuadProg solver for sparse matrix.
 */
class QuadProgSparseSolver : public SolverInterface
{
  public:
    QuadProgSparseSolver();
    QuadProgSparseSolver(int nrVar, int nrEq, int nrInEq);

    const Eigen::VectorXi &SI_iter() const override;
    int SI_fail() const override;

    const Eigen::VectorXd &SI_result() const override;

    void SI_problem(int nrVar, int nrEq, int nrInEq) override;

    /**
	 * Solve the problem.
	 * Solve the system 
	 * \f[
	 * 	 \left\{
		 	\begin{array}{cl}
			 {
				 \min_x & \frac{1}{2} x^T Q x + c^T x \\
				        & Aeq x = beq \\
						& AInEq x \leq bInEq
			 }
		 \right.
	 * \f]
	 * @param Q An N-by-N symmetric positive definite sparse matrix.
	 * @param C An N-by-1 sparse vector.
	 * @param Aeq An M-by-N sparse matrix.
	 * @param beq An M-by-1 sparse vector.
	 * @param Aineq An P-by-N sparse matrix.
	 * @param Bineq An P-by-1 sparse vector.
	 * @return The qp solver result.
	 */
    bool SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
		  const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
		  const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq) override;

  private:
    std::unique_ptr<QuadProgSparse> solver_;
};

/**
 * QuadProg solver for both dense matrix.
 */
class QLDSolver : public SolverInterface //TODO: Enable sparse matrix
{
  public:
    QLDSolver();
    QLDSolver(int nrVar, int nrEq, int nrInEq);

    const Eigen::VectorXi &SI_iter() const override;
    int SI_fail() const override;

    const Eigen::VectorXd &SI_result() const override;

    void SI_problem(int nrVar, int nrEq, int nrInEq) override;
    bool SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
		  const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
		  const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq) override;

  private:
    std::unique_ptr<QLD> solver_;
};

} // namespace pc