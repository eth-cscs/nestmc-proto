file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arbor)
if(ARB_WITH_ASSERTIONS)
    configure_file(arbor/assert_macro.hpp.enabled arbor/assert_macro.hpp COPYONLY)
else()
    configure_file(arbor/assert_macro.hpp.disabled arbor/assert_macro.hpp COPYONLY)
endif()

add_library(arbor-public-headers INTERFACE)

# At build time, public headers found in this directory and in
# binary directory (viz. assert.hpp, version.hpp). For consumers
# of an installed arbor package, the include files will be found
# in the 'include' subdirectory of the instalation path.

target_include_directories(arbor-public-headers INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Always make version.hpp-test, but only update version.hpp if it differs.

add_custom_command(
    OUTPUT _always_rebuild
    COMMAND true
)

set(arb_features)
if(ARB_WITH_ASSERTIONS)
    # define ARB_ASSERT_ENABLED in version.hpp
    list(APPEND arb_features ASSERT)
endif()
if(ARB_WITH_MPI)
    # define ARB_MPI_ENABLED in version.hpp
    list(APPEND arb_features MPI)
endif()
if(ARB_WITH_GPU)
    # define ARB_GPU_ENABLED in version.hpp
    list(APPEND arb_features GPU)
endif()
if(ARB_WITH_PROFILING)
    # define ARB_PROFILE_ENABLED in version.hpp
    list(APPEND arb_features PROFILE)
endif()

set(arb_constants)

add_custom_command(
    OUTPUT version.hpp-test
    DEPENDS _always_rebuild
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/git-source-id ${FULL_VERSION_STRING} ${arb_features} > version.hpp-test
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/git-source-id-constants ${arb_constants}  >> version.hpp-test
)

set(version_hpp_path arbor/version.hpp)
add_custom_command(
    OUTPUT ${version_hpp_path}
    DEPENDS version.hpp-test
    COMMAND ${CMAKE_COMMAND} -E copy_if_different version.hpp-test ${version_hpp_path}
)

add_custom_target(generate_version_hpp DEPENDS ${version_hpp_path})
add_dependencies(arbor-public-headers generate_version_hpp)

install(DIRECTORY arbor ${CMAKE_CURRENT_BINARY_DIR}/arbor
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp")

install(TARGETS arbor-public-headers EXPORT arbor-targets)
