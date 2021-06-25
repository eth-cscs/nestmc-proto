include(CMakeParseArguments)

# If a MODCC executable is explicitly provided, don't make the in-tree modcc a dependency.

function(build_modules)
    cmake_parse_arguments(build_modules "" "MODCC;TARGET;SOURCE_DIR;DEST_DIR;MECH_SUFFIX" "MODCC_FLAGS;GENERATES" ${ARGN})

    if("${build_modules_SOURCE_DIR}" STREQUAL "")
        set(build_modules_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    endif()

    if("${build_modules_DEST_DIR}" STREQUAL "")
        set(build_modules_SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
    file(MAKE_DIRECTORY "${build_modules_DEST_DIR}")

    set(all_generated)
    foreach(mech ${build_modules_UNPARSED_ARGUMENTS})
        set(mod "${build_modules_SOURCE_DIR}/${mech}.mod")
        set(out "${build_modules_DEST_DIR}/${mech}")
        set(generated)
        foreach (suffix ${build_modules_GENERATES})
            list(APPEND generated ${out}${suffix})
        endforeach()

        set(depends "${mod}")
        if(build_modules_MODCC)
            set(modcc_bin ${build_modules_MODCC})
        else()
            list(APPEND depends modcc)
            set(modcc_bin $<TARGET_FILE:modcc>)
        endif()

        set(flags ${build_modules_MODCC_FLAGS} -o "${out}")
        if(build_modules_MECH_SUFFIX)
            list(APPEND flags -m "${mech}${build_modules_MECH_SUFFIX}")
        endif()

        add_custom_command(
            OUTPUT ${generated}
            DEPENDS ${depends}
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMAND ${modcc_bin} ${flags} ${mod}
            COMMENT "modcc generating: ${generated}"
        )
        set_source_files_properties(${generated} PROPERTIES GENERATED TRUE)
        list(APPEND all_generated ${generated})
    endforeach()

    # Fake target to always trigger .mod -> .hpp/.cu dependencies because CMake
    if (build_modules_TARGET)
        set(depends ${all_generated})
        if(NOT build_modules_MODCC)
            list(APPEND depends modcc)
        endif()
        add_custom_target(${build_modules_TARGET} DEPENDS ${depends})
    endif()
endfunction()

function("make_catalogue")
  cmake_parse_arguments(MK_CAT "" "NAME;SOURCES;OUTPUT;PREFIX;STANDALONE;VERBOSE" "CXX_FLAGS_TARGET;MECHS" ${ARGN})
  set(MK_CAT_OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/${MK_CAT_NAME}")

  # Need to set ARB_WITH_EXTERNAL_MODCC *and* modcc
  set(external_modcc)
  if(ARB_WITH_EXTERNAL_MODCC)
    set(external_modcc MODCC ${modcc})
  endif()

  if(MK_CAT_VERBOSE)
    message("Catalogue name:       ${MK_CAT_NAME}")
    message("Catalogue mechanisms: ${MK_CAT_MECHS}")
    message("Catalogue sources:    ${MK_CAT_SOURCES}")
    message("Catalogue output:     ${MK_CAT_OUT_DIR}")
    message("Build as standalone:  ${MK_CAT_STANDALONE}")
    message("Arbor cxx flags:      ${MK_CAT_CXX_FLAGS_TARGET}")
    message("Arbor cxx compiler:   ${ARB_CXX}")
    message("Script prefix:        ${MK_CAT_PREFIX}")
    message("Current cxx compiler: ${CMAKE_CXX_COMPILER}")
  endif()

  file(MAKE_DIRECTORY "${MK_CAT_OUT_DIR}")

  if (NOT TARGET build_all_mods)
    add_custom_target(build_all_mods)
  endif()

  build_modules(
    ${MK_CAT_MECHS}
    SOURCE_DIR "${MK_CAT_SOURCES}"
    DEST_DIR "${MK_CAT_OUT_DIR}"
    ${external_modcc} # NB: expands to 'MODCC <binary>' to add an optional argument
    MODCC_FLAGS -t cpu -t gpu ${ARB_MODCC_FLAGS} -N arb::${MK_CAT_NAME}_catalogue
    GENERATES .hpp _cpu.cpp _gpu.cpp _gpu.cu
    TARGET build_catalogue_${MK_CAT_NAME}_mods)

  set(catalogue_${MK_CAT_NAME}_source ${CMAKE_CURRENT_BINARY_DIR}/${MK_CAT_NAME}_catalogue.cpp)
  set(catalogue_${MK_CAT_NAME}_options -A arbor -I ${MK_CAT_OUT_DIR} -o ${catalogue_${MK_CAT_NAME}_source} -B multicore -C ${MK_CAT_NAME} -N arb::${MK_CAT_NAME}_catalogue)
  if(ARB_WITH_GPU)
    list(APPEND catalogue_${MK_CAT_NAME}_options -B gpu)
  endif()

  add_custom_command(
    OUTPUT  ${catalogue_${MK_CAT_NAME}_source}
    COMMAND ${MK_CAT_PREFIX}/generate_catalogue ${catalogue_${MK_CAT_NAME}_options} ${MK_CAT_MECHS}
    COMMENT "Building catalogue ${MK_CAT_NAME}"
    DEPENDS ${MK_CAT_PREFIX}/generate_catalogue)

  add_custom_target(${MK_CAT_NAME}_catalogue_cpp_target DEPENDS ${catalogue_${MK_CAT_NAME}_source})
  add_dependencies(build_catalogue_${MK_CAT_NAME}_mods ${MK_CAT_NAME}_catalogue_cpp_target)
  add_dependencies(build_all_mods build_catalogue_${MK_CAT_NAME}_mods)

  foreach(mech ${MK_CAT_MECHS})
    list(APPEND catalogue_${MK_CAT_NAME}_source ${MK_CAT_OUT_DIR}/${mech}_cpu.cpp)
    if(ARB_WITH_GPU)
      list(APPEND catalogue_${MK_CAT_NAME}_source ${MK_CAT_OUT_DIR}/${mech}_gpu.cpp ${MK_CAT_OUT_DIR}/${mech}_gpu.cu)
    endif()
  endforeach()
  set(${MK_CAT_OUTPUT} ${catalogue_${MK_CAT_NAME}_source} PARENT_SCOPE)

  if(${MK_CAT_STANDALONE})
    add_library(${MK_CAT_NAME}-catalogue SHARED ${catalogue_${MK_CAT_NAME}_source})
    target_compile_definitions(${MK_CAT_NAME}-catalogue PUBLIC STANDALONE=1)
    target_compile_options(${MK_CAT_NAME}-catalogue PUBLIC ${MK_CAT_CXX_FLAGS_TARGET})
    set_target_properties(${MK_CAT_NAME}-catalogue
      PROPERTIES
      SUFFIX ".so"
      PREFIX ""
      CXX_STANDARD 17)
    target_include_directories(${MK_CAT_NAME}-catalogue PUBLIC "${MK_CAT_ARBOR}/arbor/")

    if(TARGET arbor)
      target_link_libraries(${MK_CAT_NAME}-catalogue PRIVATE arbor)
    else()
      target_link_libraries(${MK_CAT_NAME}-catalogue PRIVATE arbor::arbor)
    endif()
  endif()
endfunction()
