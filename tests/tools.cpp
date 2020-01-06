/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "tools.h"

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

std::string getSortedTimers(SolverTimers& solT)
{
    std::sort(solT.st.begin(), solT.st.end(), [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::sort(solT.bt.begin(), solT.bt.end(), [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::sort(solT.ct.begin(), solT.ct.end(), [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving speed: ";
    for (auto t : solT.st)
        ss << t.first << " (" + std::to_string(t.second) << "ms) > ";
    ss.seekp(-2, std::ios_base::end);
    ss << "\nBuild speed:";
    for (auto t : solT.bt)
        ss << t.first << " (" + std::to_string(t.second) << "ms) > ";
    ss.seekp(-2, std::ios_base::end);
    ss << "\nOverall speed:";
    for (auto t : solT.ct)
        ss << t.first << " (" + std::to_string(t.second) << "ms) > ";
    ss.seekp(-2, std::ios_base::end);
    ss.put(' ');

    return ss.str();
}

} // namespace tools
