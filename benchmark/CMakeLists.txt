add_executable(solver_benchmark
    solver_benchmark.cc
    problem_generator.cc
)

target_link_libraries(solver_benchmark PRIVATE PoseLib::PoseLib Eigen3::Eigen)

install(
  TARGETS solver_benchmark
  RUNTIME DESTINATION bin
  COMPONENT "${INSTALL_BIN_DIR}"
)

# Compilation options
set_compile_flags(solver_benchmark)

# Optional paired benchmark; reference files are generated from an explicit Git
# revision by compare_mixed_pose.py and are never part of the installed library.
set(MIXED_POSE_REFERENCE_DIR "" CACHE PATH "Directory containing reference mixed pose solver sources")
if(MIXED_POSE_REFERENCE_DIR)
    add_executable(mixed_pose_benchmark
        mixed_pose_benchmark.cc
        problem_generator.cc
        ${MIXED_POSE_REFERENCE_DIR}/p1p2ll.cc
        ${MIXED_POSE_REFERENCE_DIR}/p2p1ll.cc
    )
    target_include_directories(mixed_pose_benchmark PRIVATE ${PROJECT_SOURCE_DIR}/PoseLib/solvers)
    target_link_libraries(mixed_pose_benchmark PRIVATE PoseLib::PoseLib Eigen3::Eigen)
    set_compile_flags(mixed_pose_benchmark)
endif()
