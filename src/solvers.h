#pragma once

#include <memory>
#include <Eigen/Core>

// QP solvers
#include <eigen-qld/QLD.h>
#include <eigen-QuadProg/QuadProg.h>
// #include <eigen-lssol/LSSOL.h>

namespace pc
{

// This class act as an interface (but it is not an abstract class)
class SolverInterface
{
public:
    virtual const Eigen::VectorXi& SI_iter() const {};
	virtual int SI_fail() const {};

	virtual const Eigen::VectorXd& SI_result() const {};

	virtual void SI_problem(int nrVar, int nrEq, int nrInEq) {};
	virtual bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& C,
		const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& Bineq) {};
};

class QuadProgDenseSolver : public Solver
{
public:
    QuadProgDenseSolver();
    QuadProgDenseSolver(int nrVar, int nrEq, int nrInEq);

    const Eigen::VectorXi& SI_iter() const override;
    int SI_fail() const override;

	const Eigen::VectorXd& SI_result() const override;

	void SI_problem(int nrVar, int nrEq, int nrInEq) override;
	bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& C,
		const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& Bineq) override;


private:
    std::unique_ptr<QuadProgDense> solver_;
};

class QuadProgSparseSolver : public Solver
{
public:
    QuadProgSparseSolver();
    QuadProgSparseSolver(int nrVar, int nrEq, int nrInEq);

    const Eigen::VectorXi& SI_iter() const override;
    int SI_fail() const override;

	const Eigen::VectorXd& SI_result() const override;

	void SI_problem(int nrVar, int nrEq, int nrInEq) override;
	bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& C,
		const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& Bineq) override;


private:
    std::unique_ptr<QuadProgSparse> solver_;
};

class QLDSolver : public Solver
{
public:
    QLDSolver();
    QLDSolver(int nrVar, int nrEq, int nrInEq);

    const Eigen::VectorXi& SI_iter() const override;
    int SI_fail() const override;

	const Eigen::VectorXd& SI_result() const override;

	void SI_problem(int nrVar, int nrEq, int nrInEq) override;
	bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& C,
		const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& Bineq) override;


private:
    std::unique_ptr<QLD> solver_;
};

} // namespace pc