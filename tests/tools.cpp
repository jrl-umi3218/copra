#include "tools.h"

#include <boost/test/unit_test.hpp>
#include <sstream>

namespace tools {

Eigen::MatrixXd spanMatrix(const Eigen::MatrixXd& m, int size, int addCols)
{
    Eigen::MatrixXd mout = Eigen::MatrixXd::Zero(m.rows() * size, m.cols() * (size + addCols));
    for (int i = 0; i < size; ++i)
        mout.block(i * m.rows(), i * m.cols(), m.rows(), m.cols()) = m;
    return mout;
}

Eigen::VectorXd spanVector(const Eigen::VectorXd& v, int size)
{
    Eigen::VectorXd vout = Eigen::VectorXd::Zero(v.rows() * size);
    for (int i = 0; i < size; ++i)
        vout.segment(i * v.rows(), v.rows()) = v;
    return vout;
}

void printSortedTimers(SolverTimers& solT)
{
    using underlying_pair = std::pair<std::string, double>;
    std::sort(solT.st.begin(), solT.st.end(), [](const underlying_pair& lhs, const underlying_pair& rhs) { return lhs.second < rhs.second; });
    std::sort(solT.bt.begin(), solT.bt.end(), [](const underlying_pair& lhs, const underlying_pair& rhs) { return lhs.second < rhs.second; });
    std::sort(solT.ct.begin(), solT.ct.end(), [](const underlying_pair& lhs, const underlying_pair& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving speed: ";
    for (auto t : solT.st)
        ss << t.first << " (" + std::to_string(t.second) << "ms) > ";
    ss << "\nBuild speed:";
    for (auto t : solT.bt)
        ss << t.first << " (" + std::to_string(t.second) << "ms) > ";
    ss << "\nOverall speed:";
    for (auto t : solT.ct)
        ss << t.first << " (" + std::to_string(t.second) << "ms) > ";

    BOOST_TEST_MESSAGE(ss.str());
}

} // namespace tools