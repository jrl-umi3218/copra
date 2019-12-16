set(OPTIONAL_EIGEN_QPS
    eigen-qld
    eigen-lssol
    eigen-gurobi
    eigen-osqp
)

macro(add_optional_qps PACKAGE_LIST)
    foreach(qp ${OPTIONAL_EIGEN_QPS})
        find_package(${qp} QUIET)
        if(${${qp}_FOUND})
            list(APPEND ${PACKAGE_LIST} ${qp})
        endif()
    endforeach()
endmacro()

macro(link_optional_qps project)
    foreach(qp ${OPTIONAL_EIGEN_QPS})
        if (${${qp}_FOUND})
            target_link_libraries(${project} PUBLIC ${qp}::${qp})
        endif()
    endforeach()
endmacro()