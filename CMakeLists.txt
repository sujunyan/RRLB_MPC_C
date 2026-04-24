cmake_minimum_required(VERSION 3.10)
project(RRLB_MPC_generated)

set(EXAMPLES toyExample robotArm helicopter springMass spacecraft quadcopter)
# set(EXAMPLES toyExample robotArm)

add_subdirectory(ALADIN)
add_subdirectory(osqp)
add_subdirectory(qpOASES)
add_subdirectory(hpipm)

# Custom targets for building corresponding executables grouped by example
# add_custom_target(toyExample DEPENDS ALADIN_toyExample osqp_toyExample tinyMPC_toyExample)
# add_custom_target(robotArm DEPENDS ALADIN_robotArm osqp_robotArm tinyMPC_robotArm)
# add_custom_target(helicopter DEPENDS ALADIN_helicopter osqp_helicopter tinyMPC_helicopter)