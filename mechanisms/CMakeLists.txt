include(BuildModules.cmake)

# Define catalogues
make_catalogue(
  NAME bbp
  SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bbp"
  OUTPUT "CAT_BBP_SOURCES"
  MECHS CaDynamics_E2 Ca_HVA Ca_LVAst Ih Im K_Pst K_Tst Nap_Et2 NaTa_t NaTs2_t SK_E2 SKv3_1
  PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"
  CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
  STANDALONE FALSE
  VERBOSE ${ARB_CAT_VERBOSE})

make_catalogue(
  NAME allen
  SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/allen"
  OUTPUT "CAT_ALLEN_SOURCES"
  MECHS CaDynamics Ca_HVA Ca_LVA Ih Im Im_v2 K_P K_T Kd Kv2like Kv3_1 NaTa NaTs NaV Nap SK
  PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"
  CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
  STANDALONE FALSE
  VERBOSE ${ARB_CAT_VERBOSE})

make_catalogue(
  NAME default
  SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/default"
  OUTPUT "CAT_DEFAULT_SOURCES"
  MECHS exp2syn expsyn expsyn_stdp hh kamt kdrmt nax nernst pas
  PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"
  CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
  STANDALONE FALSE
  VERBOSE ${ARB_CAT_VERBOSE})

# Join sources
set(arbor_mechanism_sources
  ${CAT_BBP_SOURCES}
  ${CAT_ALLEN_SOURCES}
  ${CAT_DEFAULT_SOURCES}
  PARENT_SCOPE)

if(ARB_WITH_CUDA_CLANG OR ARB_WITH_HIP_CLANG)
    set_source_files_properties(${arbor_mechanism_sources} PROPERTIES LANGUAGE CXX)
endif()
