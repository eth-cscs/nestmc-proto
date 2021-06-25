# Sources:

set(arbor_sources
    arbexcept.cpp
    assert.cpp
    backends/multicore/fvm.cpp
    backends/multicore/shared_state.cpp
    communication/communicator.cpp
    communication/dry_run_context.cpp
    benchmark_cell_group.cpp
    cable_cell.cpp
    cable_cell_param.cpp
    cell_group_factory.cpp
    common_types_io.cpp
    cv_policy.cpp
    execution_context.cpp
    gpu_context.cpp
    event_binner.cpp
    fvm_layout.cpp
    fvm_lowered_cell_impl.cpp
    hardware/memory.cpp
    hardware/power.cpp
    io/locked_ostream.cpp
    io/serialize_hex.cpp
    label_resolution.cpp
    lif_cell_group.cpp
    mc_cell_group.cpp
    mechcat.cpp
    mechanism.cpp
    mechinfo.cpp
    memory/gpu_wrappers.cpp
    memory/util.cpp
    morph/embed_pwlin.cpp
    morph/label_dict.cpp
    morph/label_parse.cpp
    morph/locset.cpp
    morph/morphexcept.cpp
    morph/morphology.cpp
    morph/mprovider.cpp
    morph/place_pwlin.cpp
    morph/primitives.cpp
    morph/region.cpp
    morph/segment_tree.cpp
    morph/stitch.cpp
    merge_events.cpp
    simulation.cpp
    partition_load_balance.cpp
    profile/clock.cpp
    profile/memory_meter.cpp
    profile/meter_manager.cpp
    profile/power_meter.cpp
    profile/profiler.cpp
    schedule.cpp
    spike_event_io.cpp
    spike_source_cell_group.cpp
    s_expr.cpp
    symmetric_recipe.cpp
    threading/threading.cpp
    thread_private_spike_store.cpp
    tree.cpp
    util/hostname.cpp
    util/unwind.cpp
    version.cpp
)

if(ARB_WITH_GPU)
    list(APPEND arbor_sources
        backends/gpu/fvm.cpp
        backends/gpu/fvm.cu	
        backends/gpu/shared_state.cpp
        backends/gpu/stimulus.cu
        backends/gpu/threshold_watcher.cu
        backends/gpu/matrix_assemble.cu
        backends/gpu/matrix_fine.cu
        backends/gpu/matrix_solve.cu
        backends/gpu/multi_event_stream.cpp
        backends/gpu/multi_event_stream.cu
        backends/gpu/shared_state.cu
        backends/gpu/forest.cpp
        backends/gpu/stimulus.cu
        backends/gpu/threshold_watcher.cu
        memory/fill.cu
    )
endif()

if(ARB_WITH_MPI)
    list(APPEND arbor_sources
        communication/mpi.cpp
        communication/mpi_error.cpp
        communication/mpi_context.cpp)
endif()

# Add special target for private include directory, for use by arbor target
# and arbor unit tests. Private headers are also used for the other binaries
# until the process of splitting our private and public headers is complete.

# Because we need to add this target to the EXPORT set, and it needs to be
# installed (despite being private to arbor), we have to qualify the include
# directory with a build-time only generator expression.

add_library(arbor-private-headers INTERFACE)
target_include_directories(arbor-private-headers INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")

if(ARB_WITH_NVCC OR ARB_WITH_CUDA_CLANG)
    target_include_directories(arbor-private-headers INTERFACE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
elseif(ARB_WITH_HIP_CLANG)
    target_include_directories(arbor-private-headers INTERFACE)
endif()

install(TARGETS arbor-private-headers EXPORT arbor-targets)

# Mechanisms, generated from .mod files; sets arbor_mechanism_sources
# variable, build_all_mods target. Note: CMake source file properties are
# directory-local.

add_subdirectory(../mechanisms "${CMAKE_BINARY_DIR}/mechanisms")
set_source_files_properties(${arbor_mechanism_sources} PROPERTIES GENERATED TRUE)

if(ARB_WITH_CUDA_CLANG OR ARB_WITH_HIP_CLANG)
    set_source_files_properties(${arbor_sources} PROPERTIES LANGUAGE CXX)
    set_source_files_properties(${arbor_mechanism_sources} PROPERTIES LANGUAGE CXX)
endif()

# Library target:

add_library(arbor ${arbor_sources} ${arbor_mechanism_sources})
add_dependencies(arbor build_all_mods)

target_link_libraries(arbor PRIVATE arbor-private-deps arbor-private-headers)
target_link_libraries(arbor PUBLIC arbor-public-deps arbor-public-headers)

if(ARB_WITH_CUDA_CLANG)
    target_link_libraries(arbor PRIVATE ${CUDA_LIBRARIES})
    target_link_libraries(arbor PUBLIC  ${CUDA_LIBRARIES})
endif()

set_target_properties(arbor PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

install(TARGETS arbor EXPORT arbor-targets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

