SET(MAJOR_VERSION 1)
SET(MINOR_VERSION 2)
SET(PATCH_VERSION 0)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(PolyMaker)
enable_language(CXX)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBRARIES OFF)
SET(CMAKE_EXE_LINKER_FLAGS "")
SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")



#set(CUDA_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS})
set(CUDA_NVCC_FLAGS  ${CUDA_NVCC_FLAGS}; -O3 -gencode arch=compute_20,code=sm_20)
set(CUDA_CUDART_LIBRARY /usr/local/cuda-7.5/lib64/libcudart.so)
find_package(CUDA REQUIRED)

if (CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "-- Could not find CUDA. Program_options library")
endif()


find_package(Eigen3 3.0.0 REQUIRED)

if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(FATAL_ERROR "-- Could not find any EIGEN3 library")
endif()

find_package(Boost REQUIRED)
#find_package(Boost 1.34.0 COMPONENTS Program_options ...)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "-- Could not find BOOST.Program_options library")
endif()


cuda_add_executable(poly
    main.cu
    classes.h
    input.cpp
    input.h
    output.cpp
    output.h
    )

target_link_libraries(poly boost_program_options)

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -lm -lcuda -std=c++11 ")

IF(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_SET_DESTDIR "on")
SET(CPACK_PACKAGING_INSTALL_PREFIX "/tmp")
SET(CPACK_GENERATOR "DEB")

SET(CPACK_PACKAGE_DESCRIPTION "PolyMaker utility.")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Ploymaker utility creates polycrystalline configuration files for molecular dynamics.")
#SET(CPACK_PACKAGE_VENDOR "Vendor")
SET(CPACK_PACKAGE_CONTACT "Tymofiy Lukinov ")
SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
SET(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${MAJOR_VERSION}.${MINOR_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${MAJOR_VERSION}.${MINOR_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}")

#SET(CPACK_DEBIAN_PACKAGE_DEPENDS "kdebase-runtime (>= 4:4.2.1), kdelibs5 (>= 4:4.2.1), libc6 (>= 2.1.3), libgcc1 (>= 1:4.1.1), libplasma3, libqt4-dbus (>= 4.5.0), libqtcore4 (>= 4.5.0), libqtgui4 (>= 4.5.0), libstdc++6 (>= 4.2.1)")

SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
SET(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

SET(CPACK_COMPONENTS_ALL Libraries ApplicationData)
INCLUDE(CPack)

ENDIF(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
