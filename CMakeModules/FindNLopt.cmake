
if(DEFINED PROJECT_NAME AND PROJECT_NAME)
  string(TOUPPER ${PROJECT_NAME} _UPPER_PROJECT_NAME)
  set(_PROJECT_DEPENDENCY_DIR ${_UPPER_PROJECT_NAME}_DEPENDENCY_DIR)
else()
  # Fallback if PROJECT_NAME is not set (shouldn't happen in normal usage)
  # This is a COPASI-specific Find module, so default to COPASI_DEPENDENCY_DIR
  set(_PROJECT_DEPENDENCY_DIR COPASI_DEPENDENCY_DIR)
endif()

message(STATUS "NLOPT_DIR: $ENV{NLOPT_DIR}")
message(STATUS "PROJECT_NAME: ${PROJECT_NAME}")
message(STATUS "_PROJECT_DEPENDENCY_DIR variable name: ${_PROJECT_DEPENDENCY_DIR}")
message(STATUS "PROJECT_DEPENDENCY_DIR: ${${_PROJECT_DEPENDENCY_DIR}}")

message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "CONAN_LIB_DIRS_NLOPT: ${CONAN_LIB_DIRS_NLOPT}")
message(STATUS "COPASI_DEPENDENCY_DIR: ${COPASI_DEPENDENCY_DIR}")
message(STATUS "PATHS: $ENV{NLOPT_DIR}/${CMAKE_INSTALL_LIBDIR}/cmake")
message(STATUS "PATHS: ${${_PROJECT_DEPENDENCY_DIR}}/${CMAKE_INSTALL_LIBDIR}/cmake")
message(STATUS "PATHS: ${${_PROJECT_DEPENDENCY_DIR}}/lib/cmake")
message(STATUS "PATHS: /usr/${CMAKE_INSTALL_LIBDIR}/cmake")


find_package(NLopt CONFIG QUIET
  CONFIGS NLoptConfig.cmake
  PATHS $ENV{NLOPT_DIR}/${CMAKE_INSTALL_LIBDIR}/cmake
        ${${_PROJECT_DEPENDENCY_DIR}}/${CMAKE_INSTALL_LIBDIR}/cmake
        ${${_PROJECT_DEPENDENCY_DIR}}/lib/cmake
        /usr/${CMAKE_INSTALL_LIBDIR}/cmake
        ${CONAN_LIB_DIRS_NLOPT}/cmake
    PATH_SUFFIXES nlopt
    NO_DEFAULT_PATH
)

find_package(NLopt CONFIG QUIET
  CONFIGS NLoptConfig.cmake
  PATHS $ENV{NLOPT_DIR}/${CMAKE_INSTALL_LIBDIR}/cmake
        ${${_PROJECT_DEPENDENCY_DIR}}/${CMAKE_INSTALL_LIBDIR}/cmake
        ${${_PROJECT_DEPENDENCY_DIR}}/lib/cmake
        /usr/${CMAKE_INSTALL_LIBDIR}/cmake
        ${CONAN_LIB_DIRS_NLOPT}/cmake
    PATH_SUFFIXES nlopt
)

if (NOT NLOPT_INCLUDE_DIRS)
  
find_path(NLOPT_INCLUDE_DIRS nlopt.hpp
PATHS ${${_PROJECT_DEPENDENCY_DIR}}
      ${${_PROJECT_DEPENDENCY_DIR}}/include
      /opt/include
      CMAKE_FIND_ROOT_PATH_BOTH
NO_DEFAULT_PATH)

if (NOT NLOPT_INCLUDE_DIRS)
find_path(NLOPT_INCLUDE_DIRS nlopt.hpp)
endif (NOT NLOPT_INCLUDE_DIRS)


if (NOT NLOPT_INCLUDE_DIRS)
message(FATAL_ERROR "NLOPT include dir not found not found!")
endif (NOT NLOPT_INCLUDE_DIRS)


find_library(NLOPT_LIBRARIES 
NAMES nlopt-static 
      nlopt
      libnlopt-static 
      libnlopt
PATHS ${${_PROJECT_DEPENDENCY_DIR}}
      ${${_PROJECT_DEPENDENCY_DIR}}
      ${${_PROJECT_DEPENDENCY_DIR}}/${CMAKE_INSTALL_LIBDIR}
      ${CONAN_LIB_DIRS_NLOPT}
      ~/Library/Frameworks
      /Library/Frameworks
      /sw/lib        # Fink
      /opt/local/lib # MacPorts
      /opt/csw/lib   # Blastwave
      /opt/lib
      /usr/freeware/lib64
      CMAKE_FIND_ROOT_PATH_BOTH
      NO_DEFAULT_PATH)

if (NOT NLOPT_LIBRARIES)
find_library(NLOPT_LIBRARIES 
    NAMES nlopt-static nlopt
)
endif (NOT NLOPT_LIBRARIES)

if (NOT NLOPT_LIBRARIES)
  message(FATAL_ERROR "NLOPT library not found!")
endif (NOT NLOPT_LIBRARIES)

add_library(NLopt::nlopt UNKNOWN IMPORTED)
set_target_properties(NLopt::nlopt 
  PROPERTIES
  IMPORTED_LOCATION ${NLOPT_LIBRARIES}
  INTERFACE_INCLUDE_DIRECTORIES ${NLOPT_INCLUDE_DIRS}
)

endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NLopt REQUIRED NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES)