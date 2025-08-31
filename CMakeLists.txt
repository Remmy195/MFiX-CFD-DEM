cmake_minimum_required(VERSION 3.16)
project(MFiX LANGUAGES Fortran C)

# Options
option(ENABLE_MPI "Build with MPI (DMP) support" OFF)
option(ENABLE_OpenMP "Build with OpenMP (SMP) support" OFF)

set(_unused ${QUIET})
if(DEFINED QUIET)
  add_compile_definitions(QUIET=1)
endif()

# Read version from file
set(MFIX_SRC ${PROJECT_SOURCE_DIR}/model)
file(READ ${PROJECT_SOURCE_DIR}/VERSION MFIX_VERSION)

# Add compile options for Debug configuration globally
# Note we cannot use -fcheck=all,no-array-temps because
# this is not supported by gfortran 5.3.0 (Windows) so
# we need to use flags which are compatible with that.
add_compile_options(
    $<$<CONFIG:Debug>:-fcheck=bounds,do,mem,pointer,recursion>
)

# Set CMAKE_BUILD_TYPE (if undefined)
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_Fortran_FLAGS)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

set(SUFFIX "")
if(ENABLE_MPI)
  set(SUFFIX ${SUFFIX}_dmp)
endif()
if(ENABLE_OpenMP)
  set(SUFFIX ${SUFFIX}_smp)
endif()

if(ENABLE_POSTMFIX)
  set(NAME "postmfix")
  set(ENABLE_MPI OFF)
  set(ENABLE_OpenMP OFF)
  set(ENABLE_PYMFIX OFF)
  set(MFIX_LIB_DIR "") # TODO MFIX_LIB_DIR work for postmfix build
else()
  set(NAME "mfixsolver${SUFFIX}")
endif()

# This is a shortcut for running tests
if(MFIX_LIB_DIR) # Only supported for batch solver, for now
  unset(ENABLE_PYMFIX)
  set(ENABLE_PYMFIX OFF)
  set(LIB_NAME "mfix${SUFFIX}")
  set(LIB_FILE_NAME lib${LIB_NAME}.a)
endif()

set(TARGET ${NAME})

# build_mfixsolver sets UDF_DIR path
if(NOT DEFINED UDF_DIR)
  # if not build_mfixsolver, use current directory:
  # I would prefer not to do this but users expect this behavior :(
  set(UDF_DIR ${CMAKE_BINARY_DIR})
  message(STATUS "Setting UDF_DIR to ${UDF_DIR}")
endif()

# If UDF_DIR is relative, resolve it relative to current directory
get_filename_component(UDF_DIR "${UDF_DIR}" ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}")

# use "mod{SUFFIX}" subdir for compiled modules (-J flag)
set(CMAKE_Fortran_MODULE_DIRECTORY mod${SUFFIX})

# set include path (-I flag)
include_directories(${MFIX_SRC}/include)
include_directories(${UDF_DIR})

# Report build settings
message(STATUS "MFiX build settings summary: ")
message(STATUS "   MFiX version      = ${MFIX_VERSION}")
message(STATUS "   Build type        = ${CMAKE_BUILD_TYPE}")
message(STATUS "   CMake version     = ${CMAKE_VERSION}")
message(STATUS "   Fortran compiler  = ${CMAKE_Fortran_COMPILER}")
message(STATUS "   C compiler        = ${CMAKE_C_COMPILER}")
message(STATUS "   ENABLE_MPI        = ${ENABLE_MPI}")
message(STATUS "   ENABLE_OpenMP     = ${ENABLE_OpenMP}")
message(STATUS "   ENABLE_PYMFIX     = ${ENABLE_PYMFIX}")
if(CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} BUILDTYPE)
  message(STATUS "   Fortran flags     = ${CMAKE_Fortran_FLAGS_${BUILDTYPE}}")
else()
  message(STATUS "   Fortran flags     = ${CMAKE_Fortran_FLAGS}")
endif()

# Find MFiX project file
file(GLOB PROJECT_FILES ${UDF_DIR}/mfix.dat ${UDF_DIR}/*.mfx)
list(LENGTH PROJECT_FILES COUNT)
if(COUNT GREATER 1)
  message(FATAL_ERROR Multiple project files found: ${PROJECT_FILES})
elseif(COUNT EQUAL 1)
  list(GET PROJECT_FILES 0 PROJECT_FILE)
endif()

# Build species.inc from project file
# TODO check if species.inc is really needed
if(PROJECT_FILE)
  message(STATUS "Building custom solver for ${PROJECT_FILE}")
  add_custom_command(OUTPUT ${UDF_DIR}/species.inc
    COMMAND python ${MFIX_SRC}/rxn_preproc.py ${PROJECT_FILE}
    COMMENT "Generating species.inc"
    WORKING_DIRECTORY ${UDF_DIR})
endif()

# Avoid picking up cmake test files
function(cleanup SRC_LIST)
  set(TMP_LIST "")
  foreach(F ${${SRC_LIST}})
    if (NOT ${F} MATCHES /CMakeFiles/)
      list(APPEND TMP_LIST ${F})
    endif()
  endforeach()
  set(${SRC_LIST} ${TMP_LIST} PARENT_SCOPE)
endfunction()

# Find the sources
file(GLOB_RECURSE SRC_F90 ${MFIX_SRC}/[a-zA-Z]*.f ${MFIX_SRC}/[a-zA-Z]*.f90 ${MFIX_SRC}/[a-zA-Z]*.f95)
file(GLOB_RECURSE SRC_F77 ${MFIX_SRC}/[a-zA-Z]*.F)
file(GLOB_RECURSE SRC_C ${MFIX_SRC}/[a-zA-Z]*.c)
file(GLOB UDF_F90 ${UDF_DIR}/[a-zA-Z]*.f ${UDF_DIR}/[a-zA-Z]*.f90 ${UDF_DIR}/[a-zA-Z]*.f95)
file(GLOB UDF_F77 ${UDF_DIR}/[a-zA-Z]*.F)
file(GLOB UDF_C ${UDF_DIR}/[a-zA-Z]*.c)
cleanup(SRC_F90)
cleanup(SRC_F77)
cleanup(SRC_C)
cleanup(UDF_F90)
cleanup(UDF_F77)
cleanup(UDF_C)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  LIST(REMOVE_ITEM SRC_C ${MFIX_SRC}/gfortran_init.c)
endif()
if(APPLE)
  LIST(REMOVE_ITEM SRC_C ${MFIX_SRC}/xpow.c)
endif()
cleanup(SRC_F90)
cleanup(SRC_F77)
cleanup(SRC_C)
cleanup(UDF_F90)
cleanup(UDF_F77)
cleanup(UDF_C)

# Helper function to append to COMPILE_DEFINITIONS
function(append_source_compile_definitions source_file)
  get_source_file_property(defs ${source_file} COMPILE_DEFINITIONS)
  # Check if the property exists (it will be NOTFOUND if not previously set)
  if("${defs}" STREQUAL "NOTFOUND" OR "${defs}" STREQUAL "")
    set(defs "")
  endif()
  # Append all additional arguments as definitions
  foreach(def ${ARGN})
    list(APPEND defs "${def}")
  endforeach()
  # Set the updated property
  set_source_files_properties(${source_file}
    PROPERTIES COMPILE_DEFINITIONS "${defs}")
endfunction()


# Replace source files with modified copies in project directory (UDF)
macro(merge_udf SRC UDF)
  foreach(F ${${SRC}})
    cmake_path(GET F FILENAME F1)
    foreach(G ${${UDF}})
      cmake_path(GET G FILENAME G1)
      if(F1 STREQUAL G1)
        list(REMOVE_ITEM ${SRC} ${F})
        break()
      endif()
    endforeach()
  endforeach()
  list(APPEND ${SRC} ${${UDF}})
endmacro()

#message(STATUS "PRE-MERGE")
#message(STATUS "F90 sources: " ${SRC_F90})
#message(STATUS "F77 sources: " ${SRC_F77})
#message(STATUS "C sources: " ${SRC_C})

if (NOT ENABLE_POSTMFIX AND NOT MFIX_LIB_DIR)
  merge_udf(SRC_F90 UDF_F90)
  merge_udf(SRC_F77 UDF_F77)
  merge_udf(SRC_C UDF_C)
endif()

#message(STATUS "POST-MERGE")
#message(STATUS "F90 sources: " ${SRC_F90})
#message(STATUS "F77 sources: " ${SRC_F77})
#message(STATUS "C sources: " ${SRC_C})

set(UDFS "")
foreach(F ${UDF_F90} ${UDF_F77} ${UDF_C})
  file(RELATIVE_PATH REL ${UDF_DIR} ${F})
  list(APPEND UDFS ${REL})
endforeach()

if (NOT "${UDFS}" STREQUAL "")
  message(STATUS "Found UDFS: ${UDFS}")
endif()


# Ugh, glob is case-insensitive on MacOS and Windows
set(SRC_F77_TMP ${SRC_F77})
foreach(FILENAME ${SRC_F77_TMP})
  if(NOT ${FILENAME} MATCHES "\\.F")
    list(REMOVE_ITEM SRC_F77 ${FILENAME})
  endif()
endforeach()
set(SRC_F90_TMP ${SRC_F90})
foreach(FILENAME ${SRC_F90_TMP})
  if(${FILENAME} MATCHES "\\.F")
    list(REMOVE_ITEM SRC_F90 ${FILENAME})
  endif()
endforeach()

set(SRC_ALL ${SRC_F90} ${SRC_F77} ${SRC_C})

# TODO check if species.inc is really needed
if(PROJECT_FILE)
  list(APPEND SRC_ALL ${UDF_DIR}/species.inc)
endif()

# Set optimization flags
if (APPLE)
  set (MARCH "")
else()
  if (DEFINED ENV{CONDA_BUILD})
    set (MARCH "-march=haswell")
  else()
    set (MARCH "-march=native")
  endif()
endif()
if (DEFINED ENV{MARCH})
  set (MARCH "-march=$ENV{MARCH}")
endif()

# Trap FPE errors
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag(-ffpe-trap=invalid,zero,overflow ffpe_trap)
if (ffpe_trap)
  set(ffpe_trap "-ffpe-trap=invalid,zero,overflow")
endif()

# Don't spew irrelevant warnings at exit
# https://mfix.netl.doe.gov/forum/t/default-solver-not-available-for-mfix-gui-20-1-2/1641/12
if(NOT APPLE)
  check_fortran_compiler_flag(-ffpe-summary=none ffpe_summary)
  if(ffpe_summary)
    set(ffpe_summary -ffpe-summary=none)
  endif()
endif()

set(ffpe_flag "${ffpe_trap};${ffpe_summary}")

# Set bigendian
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(endian_flag "-fconvert=big-endian")
  # remove line length limit for gfortran
  set(line_length_flag "-ffree-line-length-0")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
  set(endian_flag "-Mextend;-byteswapio")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(endian_flag "-convert;big_endian;-assume;byterecl")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  set(endian_flag "-Mextend;-byteswapio")
endif()

#set (OPTIMIZE "-O3;-ggdb;${MARCH}")
set (FFLAGS "${ffpe_flag};${endian_flag};${line_length_flag}") # ;${OPTIMIZE}")

set_source_files_properties(${SRC_F90} PROPERTIES
  Fortran_FORMAT "FREE"
  Fortran_PREPROCESS ON
  COMPILE_OPTIONS "${FFLAGS}")

set_source_files_properties(${SRC_F77} PROPERTIES
  Fortran_FORMAT "FIXED"
  Fortran_PREPROCESS ON
  COMPILE_OPTIONS "-std=legacy;${FFLAGS}")

set_source_files_properties(${SRC_C} PROPERTIES
  COMPILE_OPTIONS "${CFLAGS}") #;${OPTIMIZE}")

set_source_files_properties(${MFIX_SRC}/xpow.c PROPERTIES
  COMPILE_OPTIONS "${CFLAGS};-O3;-ggdb;${MARCH}")

set_source_files_properties(${MFIX_SRC}/DGTSV.f PROPERTIES
  COMPILE_OPTIONS "${FFLAGS};-O3;-ggdb;${MARCH}")

# set_source_files_properties(${MFIX_SRC}/ODEPACK.F PROPERTIES
#   COMPILE_OPTIONS "${FFLAGS};-O3;-ggdb;${MARCH}")

# Is it safe to enable -O3 for all sources?  If not, we can do this
# for files where we know it matters (DGTSV.f, xpow.c, ODEPACK.F, ...)
#set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3 -ggdb ${MARCH}")
#set (CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -ggdb ${MARCH}")

if(ENABLE_PYMFIX)
  # We might be able to relax this now that we no longer use f2py
  if(NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    message(FATAL_ERROR "ENABLE_PYMFIX is only supported for GNU Fortran compiler")
  endif() #NOT GNU
  add_library(${TARGET} SHARED ${SRC_ALL})
  target_compile_definitions(${TARGET} PUBLIC PYMFIX)
  target_compile_definitions(${TARGET} PUBLIC QUIET)
  if (PROJECT_FILE)
     target_compile_definitions(${TARGET} PUBLIC SPECIES_INC)
   endif()
  set_property(TARGET ${TARGET} PROPERTY POSITION_INDEPENDENT_CODE ON)

  # Force needed symbols to be exported
  if (NOT APPLE)
    set_target_properties(${TARGET} PROPERTIES
      LINK_FLAGS "-u __main_MOD_print_flags -u __main_MOD_run_mfix0 -Wl,--wrap=pow")
  else(APPLE)
    # force all symbols  into mfixsolver.so
    target_link_libraries(${TARGET} PUBLIC "-undefined dynamic_lookup -all_load")
  endif(NOT APPLE)

  # prevent "lib" (mfixsolver.so not libmfixsolver.so)xs
  set_target_properties(${TARGET} PROPERTIES
    PREFIX ""
    OUTPUT_NAME ${TARGET}
    SUFFIX ".so")

elseif(ENABLE_POSTMFIX)
  # Postmfix isn't ready for the "source replacement" strategy so we
  #  have to revert to building a static lib and over-linking
  message(STATUS "Building postmfix")
  # TODO use lib_dir here
  add_library(mfix STATIC ${SRC_ALL})
  target_compile_definitions(mfix PUBLIC MFIX_VERSION=\"${MFIX_VERSION}\") # run_mod.f expects quotes
  if (MINGW)
    target_compile_definitions(mfix PUBLIC MINGW)
  endif()
  file(GLOB_RECURSE POSTMFIX_SRC ${MFIX_SRC}/../tools/postmfix/[a-zA-Z]*.f ${MFIX_SRC}/../tools/postmfix/[a-zA-Z]*.f90)
  add_executable(postmfix ${POSTMFIX_SRC})
  target_link_libraries(postmfix mfix)
  target_link_options(postmfix PUBLIC "-Wl,--allow-multiple-definition") # yuk
  set_source_files_properties(${POSTMFIX_SRC} PROPERTIES
    Fortran_FORMAT "FREE"
    Fortran_PREPROCESS ON
    COMPILE_OPTIONS "${FFLAGS}")

elseif(MFIX_LIB_DIR) # Special accelerated build for tests, do not use this
  # except for running tests
  add_executable(${TARGET} dummy.f ${UDF_F90} ${UDF_F77} ${UDF_C})
  if(PROJECT_FILE)
    target_sources(${TARGET} PUBLIC ${UDF_DIR}/species.inc)
    target_compile_definitions(${TARGET} PUBLIC SPECIES_INC)
  endif()
  target_link_libraries(${TARGET} PUBLIC ${LIB_NAME})
  target_link_options(${TARGET} PUBLIC -L ${MFIX_LIB_DIR})

  if(EXISTS ${MFIX_LIB_DIR}/${LIB_FILE_NAME})
    message(STATUS "Using prebuilt library " ${MFIX_LIB_DIR}/${LIB_FILE_NAME})
    # Dependency cannot be a directory, so we use 'run.mod' as placeholder
    # Copy modules from lib dir to local, in case any get overwritten
    add_custom_command(OUTPUT mod${SUFFIX}/run.mod
      COMMENT "Copying modules from ${MFIX_LIB_DIR}/mod${SUFFIX}"
      COMMAND ${CMAKE_COMMAND} -E copy_directory ${MFIX_LIB_DIR}/mod${SUFFIX} mod${SUFFIX})

    # Force rebuild of usr_read_namelist if usrnlst.inc present
    if(EXISTS ${UDF_DIR}/usrnlst.inc)
      add_custom_command(OUTPUT ${UDF_DIR}/usr_read_namelist.f
        COMMENT "Copying usr_read_namelist.f"
        COMMAND ${CMAKE_COMMAND} -E copy ${MFIX_SRC}/usr_read_namelist.f ${UDF_DIR})
      list(APPEND UDF_F90 ${UDF_DIR}/usr_read_namelist.f)
      target_sources(${TARGET} PUBLIC ${UDF_DIR}/usr_read_namelist.f)
    endif()

    if(PROJECT_FILE)
      target_sources(${TARGET} PUBLIC mod${SUFFIX}/run.mod)
    endif()
#   set(CMAKE_Fortran_MODULE_DIRECTORY ${MFIX_LIB_DIR}/mod${SUFFIX}
  else()
    message(STATUS "Building library " ${LIB_FILE_NAME})
    set(BUILD_MFIX_LIB TRUE)
    add_library(${LIB_NAME} STATIC ${SRC_ALL})
    add_custom_command(
      TARGET ${LIB_NAME}
      POST_BUILD
      COMMENT "Copying library ${LIB_NAME} to ${MFIX_LIB_DIR}"
      COMMAND ${CMAKE_COMMAND} -E copy ${LIB_FILE_NAME} ${MFIX_LIB_DIR})
    if(NOT EXISTS ${MFIX_LIB_DIR}/mod)
      add_custom_command(
        TARGET ${LIB_NAME}
        POST_BUILD
        COMMENT "Copying modules to ${MFIX_LIB_DIR}/mod${SUFFIX}"
        COMMAND ${CMAKE_COMMAND} -E copy_directory mod${SUFFIX} ${MFIX_LIB_DIR}/mod${SUFFIX})
      endif()
  endif()
  file(TOUCH dummy.f)

  # FIXME: DRY
  set_source_files_properties(${UDF_F90} PROPERTIES
    Fortran_FORMAT "FREE"
    Fortran_PREPROCESS ON
    COMPILE_OPTIONS "${FFLAGS}")

  set_source_files_properties(${UDF_F77} PROPERTIES
    Fortran_FORMAT "FIXED"
    Fortran_PREPROCESS ON
    COMPILE_OPTIONS "-std=legacy;${FFLAGS}")

  set_source_files_properties(${UDF_C} PROPERTIES
    COMPILE_OPTIONS "${CFLAGS}") #;${OPTIMIZE}")


else() # Default batch solver build
  add_executable(${TARGET} ${SRC_ALL})
  if (NOT APPLE)
    set_target_properties(${TARGET} PROPERTIES
      LINK_FLAGS "-Wl,--wrap=pow")
  endif()
  if (PROJECT_FILE)
    target_compile_definitions(${TARGET} PUBLIC SPECIES_INC)
  endif()
endif()

function (set_compile_definitions X)
  # Version
  target_compile_definitions(${X} PUBLIC MFIX_VERSION=\"${MFIX_VERSION}\")
  # usrnlst.inc
  if(EXISTS ${UDF_DIR}/usrnlst.inc)
    target_compile_definitions(${X} PUBLIC USR_NAMELIST)
  endif()

  # MinGW, ie Windows
  if(MINGW)
    target_compile_definitions(${X} PUBLIC MINGW)
  endif()

  # Define BURCAT.THR as path to Burcat database file
  target_compile_definitions(${X} PUBLIC
    BURCAT_THR=\"${PROJECT_SOURCE_DIR}/model/thermochemical/BURCAT.THR\" )
  # MPI
  if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_compile_options(${X} PUBLIC ${MPI_Fortran_COMPILE_FLAGS} )
    target_compile_definitions(${X} PUBLIC MPI )
    target_include_directories(${X} PUBLIC ${MPI_Fortran_INCLUDE_PATH} )
    target_link_libraries(${X} PUBLIC ${MPI_Fortran_LIBRARIES})
  endif()
  # OpenMP
  if(ENABLE_OpenMP)
    find_package(OpenMP REQUIRED)
    target_compile_options(${X} PUBLIC ${OpenMP_Fortran_FLAGS} )
    target_link_libraries(${X} PUBLIC ${OpenMP_Fortran_LIBRARIES} ${OpenMP_Fortran_FLAGS} )
  endif()
endfunction()

set_compile_definitions(${TARGET})
if(BUILD_MFIX_LIB)
  set_compile_definitions(${LIB_NAME})
endif()

if(ENABLE_PYMFIX)
  # NB, find_package(PythonInterp) finds /usr/bin/python, rather than the
  #  python in the active Conda environment
  #execute_process(COMMAND python -c "import sys; print(sys.executable)"
  #  OUTPUT_VARIABLE PYTHON_INTERP
  #  OUTPUT_STRIP_TRAILING_WHITESPACE)
  #message(STATUS "Found Python: " ${PYTHON_INTERP})

  if(DEFINED ENV{CONDA_BUILD})
    if(WIN32)
      set(Python3_EXECUTABLE "$ENV{PREFIX}/python")
    else()
      set(Python3_EXECUTABLE "$ENV{PREFIX}/bin/python")
    endif()
  elseif(DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}/bin/python")
    set(Python3_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")
  elseif(DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}/python.exe")
    set(Python3_EXECUTABLE "$ENV{CONDA_PREFIX}/python")
  else()
    find_package(Python3 "3.11" EXACT)
  endif()

  if(DEFINED ENV{CONDA_BUILD})
    set(SO_DIR "$ENV{PREFIX}")
    if(WIN32)
      set(SO_ABSPATH "${SO_DIR}/Scripts/${NAME}.so")
    else()
      set(SO_ABSPATH "${SO_DIR}/bin/${NAME}.so")
    endif()
  else()
    set(SO_DIR ${CMAKE_INSTALL_PREFIX})
    set(SO_ABSPATH "${SO_DIR}/${NAME}.so")
  endif()

  # Write the script file to start mfixsolver
  if(WIN32)
    set(MFIXSOLVER_SCRIPT "${NAME}.bat")
    file(WRITE ${CMAKE_BINARY_DIR}/${MFIXSOLVER_SCRIPT}
      "call ${Python3_EXECUTABLE} -m mfixgui.pymfix -m \"${SO_ABSPATH}\" %*
")
  else()
    set(MFIXSOLVER_SCRIPT "${NAME}.sh")
    file(WRITE ${CMAKE_BINARY_DIR}/${MFIXSOLVER_SCRIPT}
      "#!/bin/sh
${Python3_EXECUTABLE} -m mfixgui.pymfix -m \"${SO_ABSPATH}\" \"\$@\"
")
  endif(WIN32)
  # Ship it!
  install(TARGETS ${TARGET} DESTINATION "${CMAKE_INSTALL_PREFIX}")
  install(
    PROGRAMS ${CMAKE_BINARY_DIR}/${MFIXSOLVER_SCRIPT}
    DESTINATION ${CMAKE_INSTALL_PREFIX})
else()  # batch solver
  #      if(NOT ENABLE_POSTMFIX)
  install(TARGETS ${TARGET} DESTINATION ${CMAKE_INSTALL_PREFIX})
  #      endif()
endif()
