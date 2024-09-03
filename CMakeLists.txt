# For Intel MKL, set MKLROOT=<MKL_path> when running cmake
#   e.g. MKLROOT=/opt/mkl/ cmake -S regenie_dir/ -B regenie_dir/build/
# For OpenBLAS, set OPENBLAS_ROOT=<OpenBLAS_path> when running cmake
#   note: it also requires lapacke library
# For static compilation on Linux systems, set STATIC=1 when running cmake
#   -> this excludes GLIBC


cmake_minimum_required(VERSION 3.13)

# detect OS architecture
execute_process(
  COMMAND uname -s
  OUTPUT_VARIABLE UNAME_S
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# Get Regenie version
file(STRINGS "VERSION" RG_VERSION)

project(regenie
  VERSION ${RG_VERSION}
  )

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Ensures -std=c++11

######################################
######## check input variables

# check BGEN_PATH
if("$ENV{BGEN_PATH}" STREQUAL "")
  message( FATAL_ERROR "Must specify path to BGEN library in 'BGEN_PATH'")
else()
  set(BGEN_PATH "$ENV{BGEN_PATH}" CACHE INTERNAL "Set BGEN library path")
  if (NOT EXISTS ${BGEN_PATH})
    message( FATAL_ERROR "Specified BGEN library directory '${BGEN_PATH}' does not exist")
  endif()
endif()

# check for static compilation
if($ENV{STATIC})
  set(BUILD_STATIC ON CACHE INTERNAL "Static compilation")
  message( STATUS "Static compilation mode")
endif()

# check Boost IOStreams
if($ENV{HAS_BOOST_IOSTREAM})
  set(HAS_BOOST_IOSTREAM 1 CACHE INTERNAL "Add Boost IO")
  set(RG_VERSION "${RG_VERSION}.gz" CACHE INTERNAL "Update version")
else()
  set(HAS_BOOST_IOSTREAM 0 CACHE INTERNAL "Skip Boost IO")
endif()

# check MKL
if(NOT "$ENV{MKLROOT}" STREQUAL "")
  set(MKLROOT "$ENV{MKLROOT}" CACHE INTERNAL "Set MKL library path")
  if (NOT EXISTS ${MKLROOT})
    message( FATAL_ERROR "Specified MKL library directory '${MKLROOT}' does not exist")
  endif()
  message( STATUS "Will compile with Intel MKL library")
endif()

# check HTSlib
if(NOT "$ENV{HTSLIB_PATH}" STREQUAL "")
  set(HTSLIB_PATH "$ENV{HTSLIB_PATH}" CACHE INTERNAL "Set HTSlib library path")
  if (NOT EXISTS ${HTSLIB_PATH})
    message( FATAL_ERROR "Specified HTSlib library directory '${HTSLIB_PATH}' does not exist")
  endif()
  message( STATUS "Will compile with HTSlib")
endif()

# check OpenBLAS
if(NOT "$ENV{OPENBLAS_ROOT}" STREQUAL "")
  set(OPENBLAS_ROOT "$ENV{OPENBLAS_ROOT}" CACHE INTERNAL "Set OpenBLAS library path")
  if (NOT EXISTS ${OPENBLAS_ROOT})
    message( FATAL_ERROR "Specified OpenBLAS library directory '${OPENBLAS_ROOT}' does not exist")
  endif()
  message( STATUS "Will compile with OpenBLAS library")
endif()

######################################
######## set flags and required libraries

set(BLA_STATIC               ${BUILD_STATIC})
set(Boost_USE_STATIC_LIBS    ${BUILD_STATIC})
set(Boost_USE_DEBUG_LIBS     OFF)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)

# list each file specifically
add_executable(regenie
  ${CMAKE_SOURCE_DIR}/src/Data.cpp
  ${CMAKE_SOURCE_DIR}/src/Files.cpp
  ${CMAKE_SOURCE_DIR}/src/Geno.cpp
  ${CMAKE_SOURCE_DIR}/src/HLM.cpp
  ${CMAKE_SOURCE_DIR}/src/Interaction.cpp
  ${CMAKE_SOURCE_DIR}/src/Joint_Tests.cpp
  ${CMAKE_SOURCE_DIR}/src/Masks.cpp
  ${CMAKE_SOURCE_DIR}/src/NNLS.cpp
  ${CMAKE_SOURCE_DIR}/src/Pheno.cpp
  ${CMAKE_SOURCE_DIR}/src/Regenie.cpp
  ${CMAKE_SOURCE_DIR}/src/SKAT.cpp
  ${CMAKE_SOURCE_DIR}/src/Step1_Models.cpp
  ${CMAKE_SOURCE_DIR}/src/Step2_Models.cpp
  ${CMAKE_SOURCE_DIR}/src/MultiTrait_Tests.cpp
  ${CMAKE_SOURCE_DIR}/src/MCC.cpp
  ${CMAKE_SOURCE_DIR}/src/Ordinal.cpp
  ${CMAKE_SOURCE_DIR}/src/survival_data.cpp
  ${CMAKE_SOURCE_DIR}/src/cox_ridge.cpp
  ${CMAKE_SOURCE_DIR}/src/cox_score.cpp
  ${CMAKE_SOURCE_DIR}/src/cox_firth.cpp
  )
target_include_directories(regenie PRIVATE ${CMAKE_SOURCE_DIR}/src)

set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic -ffast-math -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wno-long-long -Wno-c11-extensions -fPIC")
add_definitions(-DVERSION_NUMBER="${RG_VERSION}")

if("${UNAME_S}" STREQUAL "Linux")
  find_package(OpenMP REQUIRED)
  target_link_libraries(regenie PRIVATE OpenMP::OpenMP_CXX)
  if(${BUILD_STATIC})
    target_link_options(regenie BEFORE PRIVATE -static-libgcc PRIVATE -static-libstdc++)
  endif()
elseif("${UNAME_S}" STREQUAL "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

set(EXTERN_LIBS_PATH "${CMAKE_SOURCE_DIR}/external_libs")
target_include_directories(regenie PRIVATE ${EXTERN_LIBS_PATH}/)

# BGEN library and its dependencies
find_library(ZSTD_LIBRARY libzstd.a HINTS "${BGEN_PATH}/build/3rd_party/zstd-1.1.0" REQUIRED)
find_library(DB_LIBRARY libdb.a HINTS "${BGEN_PATH}/build/db" REQUIRED)
find_library(SQLITE3_LIBRARY libsqlite3.a HINTS "${BGEN_PATH}/build/3rd_party/sqlite3" REQUIRED)
find_library(Boost_LIBRARY libboost.a HINTS "${BGEN_PATH}/build/3rd_party/boost_1_55_0" REQUIRED)
find_library(BGEN_LIBRARY libbgen.a HINTS "${BGEN_PATH}/build" REQUIRED)
target_link_libraries(regenie PRIVATE ${ZSTD_LIBRARY} ${BGEN_LIBRARY} ${DB_LIBRARY} ${SQLITE3_LIBRARY} ${Boost_LIBRARY})
target_include_directories(regenie PRIVATE ${BGEN_PATH} ${BGEN_PATH}/genfile/include/ ${BGEN_PATH}/3rd_party/boost_1_55_0/ ${BGEN_PATH}/3rd_party/zstd-1.1.0/lib ${BGEN_PATH}/db/include/ ${BGEN_PATH}/3rd_party/sqlite3)

# MVTNorm library
set(MVTN_PATH "${EXTERN_LIBS_PATH}/mvtnorm")
add_custom_target(
   libMvtnorm
   COMMAND make
   WORKING_DIRECTORY ${MVTN_PATH}
)
target_link_libraries(regenie PRIVATE ${MVTN_PATH}/libMvtnorm.a)
add_dependencies(regenie libMvtnorm)

# QF library
set(QF_PATH "${EXTERN_LIBS_PATH}/qf")
add_custom_target(
   libqf
   COMMAND make
   WORKING_DIRECTORY ${QF_PATH}
)
target_link_libraries(regenie PRIVATE ${QF_PATH}/qf.a)
add_dependencies(regenie libqf)

# Quadpack library
set(QUAD_PATH "${EXTERN_LIBS_PATH}/quadpack")
add_custom_target(
   libquad
   COMMAND make
   WORKING_DIRECTORY ${QUAD_PATH}
)
target_link_libraries(regenie PRIVATE ${QUAD_PATH}/libquad.a)
add_dependencies(regenie libquad)

# PGEN library
set(PGEN_PATH "${EXTERN_LIBS_PATH}/pgenlib")
add_custom_target(
   pgenlib
   COMMAND make
   WORKING_DIRECTORY ${PGEN_PATH}
)
target_link_libraries(regenie PRIVATE ${PGEN_PATH}/pgenlib.a)
target_include_directories(regenie PRIVATE ${PGEN_PATH} ${PGEN_PATH}/simde/ ${PGEN_PATH}/include/)
add_dependencies(regenie pgenlib)

# REMETA library
if(EXISTS ${HTSLIB_PATH})
  set(REMETA_PATH "${EXTERN_LIBS_PATH}/remeta")
  add_custom_target(
     remeta
     COMMAND make HTSLIB_PATH=${HTSLIB_PATH}
     WORKING_DIRECTORY ${REMETA_PATH}
  )
  target_link_libraries(regenie PUBLIC ${REMETA_PATH}/remeta.a)
  target_include_directories(regenie PUBLIC ${REMETA_PATH})
  add_dependencies(regenie remeta)

  add_definitions(-DWITH_HTSLIB)
  find_library(HTSLIB libhts.a HINTS ${HTSLIB_PATH})
  find_library(BZ2_LIB bz2 REQUIRED)
  find_library(LZMA_LIB lzma REQUIRED)
  find_library(CURL_LIB curl REQUIRED)
  find_library(CRYPTO_LIB crypto REQUIRED)
  target_link_libraries(regenie PUBLIC 
    ${HTSLIB}
    ${BZ2_LIB} ${LZMA_LIB} ${CURL_LIB} ${CRYPTO_LIB} 
  )
endif()

# Intel MKL
if(EXISTS ${MKLROOT})
  add_definitions(-DWITH_MKL -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
  target_include_directories(regenie PRIVATE ${MKLROOT}/include/)
  if(${BUILD_STATIC}) # specify static libs
    find_library(MKL_LP64_LIB libmkl_intel_lp64.a 
      HINTS "${MKLROOT}/lib/intel64"
      "${MKLROOT}/lib" 
      REQUIRED)
    find_library(MKL_THREAD_LIB libmkl_gnu_thread.a 
      HINTS "${MKLROOT}/lib/intel64"
      "${MKLROOT}/lib" 
      REQUIRED)
    find_library(MKL_CORE_LIB libmkl_core.a 
      HINTS "${MKLROOT}/lib/intel64"
      "${MKLROOT}/lib" 
      REQUIRED)
    target_link_libraries(regenie PRIVATE "-Wl,--start-group" ${MKL_LP64_LIB} ${MKL_THREAD_LIB} ${MKL_CORE_LIB} "-Wl,--end-group" -lgomp)
  else() # use dynamic libs
    find_library(MKL_LP64_LIB mkl_intel_lp64 
      PATHS "${MKLROOT}/lib/intel64"
      "${MKLROOT}/lib" 
      REQUIRED)
    find_library(MKL_THREAD_LIB mkl_gnu_thread 
      PATHS "${MKLROOT}/lib/intel64"
      "${MKLROOT}/lib" 
      REQUIRED)
    find_library(MKL_CORE_LIB mkl_core 
      PATHS "${MKLROOT}/lib/intel64"
      "${MKLROOT}/lib" 
      REQUIRED)
    target_link_libraries(regenie PRIVATE "-Wl,--no-as-needed" ${MKL_LP64_LIB} ${MKL_THREAD_LIB} ${MKL_CORE_LIB} -lgomp)
  endif()
elseif(EXISTS ${OPENBLAS_ROOT}) # OpenBLAS
  add_definitions(-DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
  target_include_directories(regenie PRIVATE ${OPENBLAS_ROOT}/include/)
  find_library(LAPACK_LIB lapack REQUIRED)
  find_library(BLAS_LIB openblas HINTS "${OPENBLAS_ROOT}/lib/" REQUIRED)
  target_link_libraries(regenie PRIVATE ${LAPACK_LIB} -llapacke ${BLAS_LIB})
endif()

# cxxopts (header-only)
target_include_directories(regenie PRIVATE ${EXTERN_LIBS_PATH}/cxxopts/include/)

# LBFGS (header-only)
target_include_directories(regenie PRIVATE ${EXTERN_LIBS_PATH}/LBFGSpp/include/)

# Eigen (header-only)
target_include_directories(regenie PRIVATE ${EXTERN_LIBS_PATH}/eigen-3.4.0/)

# Boost IO
if(${HAS_BOOST_IOSTREAM})
  if("${UNAME_S}" STREQUAL "Darwin")
    find_library(BOOST_LIB_IO libboost_iostreams libboost_iostreams.a REQUIRED)
    target_link_libraries(regenie PRIVATE ${BOOST_LIB_IO})
  elseif(${BUILD_STATIC})
    find_library(BOOST_LIB_IO libboost_iostreams.a REQUIRED)
    target_link_libraries(regenie PRIVATE ${BOOST_LIB_IO})
  else()
    target_link_libraries(regenie PRIVATE -lboost_iostreams)
  endif()
  add_definitions(-DHAS_BOOST_IOSTREAM)
  message( STATUS "Will compile with Boost Iostreams library")
endif()

# Other libraries
find_library(ZLIB_LIBRARY libz.a z REQUIRED)
find_library(M_LIB m REQUIRED)
find_library(DL_LIB dl REQUIRED)
if("${UNAME_S}" STREQUAL "Linux")
  set(GFORTRAN_LIBRARY "-lgfortran")
elseif("${UNAME_S}" STREQUAL "Darwin")
  find_library(GFORTRAN_LIBRARY gfortran REQUIRED)
endif()
target_link_libraries(
  regenie PRIVATE
  ${ZLIB_LIBRARY} ${M_LIB} ${DL_LIB} ${PTHREAD_LIB}
  ${GFORTRAN_LIBRARY}
)

install(TARGETS regenie RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
message( STATUS "REGENIE v" ${RG_VERSION})

add_custom_target(full-clean
  COMMAND cd "${MVTN_PATH}" && make clean
  COMMAND cd "${QF_PATH}" && make clean
  COMMAND cd "${QUAD_PATH}" && make clean
  COMMAND cd "${PGEN_PATH}" && make clean
  )
