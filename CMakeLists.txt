cmake_minimum_required(VERSION 3.21)
project(metric_field)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Eigen3 QUIET)
if (NOT EIGEN3_FOUND)
    message("Eigen3 not found, downloading...")
    include(FetchContent)
    FetchContent_Declare(
            eigen3
            URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
    )
    FetchContent_Populate(eigen3)

    add_library(eigen INTERFACE)
    add_library(Eigen3::Eigen ALIAS eigen)
    target_include_directories(eigen INTERFACE
            $<BUILD_INTERFACE:${eigen3_SOURCE_DIR}>
            )
endif ()

# OVM
include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

find_package(OpenVolumeMesh::OpenVolumeMesh QUIET)
if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    FetchContent_Declare(ovm
        GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git
        GIT_TAG v3.3
        SOURCE_DIR "${EXTERNAL_DIR}/OpenVolumeMesh"
    )
    FetchContent_MakeAvailable(ovm)
    
endif()


# Gaussian Processes
include(FetchContent)
FetchContent_Declare(
gauss_process
GIT_REPOSITORY https://github.com/andreacasalino/GaussianProcesses
GIT_TAG        main
)
FetchContent_MakeAvailable(gauss_process)

# CGAL for making bounding and enclosing mesh

find_package(CGAL REQUIRED)


# Set warnings
set(MSVC_WARNINGS
        /W4 # Baseline reasonable warnings
        /w14242 # 'identifier': conversion from 'type1' to 'type1', possible loss
        # of data
        /w14254 # 'operator': conversion from 'type1:field_bits' to
        # 'type2:field_bits', possible loss of data
        /w14263 # 'function': member function does not override any base class
        # virtual member function
        /w14265 # 'classname': class has virtual functions, but destructor is not
        # virtual instances of this class may not be destructed correctly
        /w14287 # 'operator': unsigned/negative constant mismatch
        /we4289 # nonstandard extension used: 'variable': loop control variable
        # declared in the for-loop is used outside the for-loop scope
        /w14296 # 'operator': expression is always 'boolean_value'
        /w14311 # 'variable': pointer truncation from 'type1' to 'type2'
        /w14545 # expression before comma evaluates to a function which is missing
        # an argument list
        /w14546 # function call before comma missing argument list
        /w14547 # 'operator': operator before comma has no effect; expected
        # operator with side-effect
        /w14549 # 'operator': operator before comma has no effect; did you intend
        # 'operator'?
        /w14555 # expression has no effect; expected expression with side- effect
        /w14619 # pragma warning: there is no warning number 'number'
        /w14640 # Enable warning on thread un-safe static member initialization
        /w14826 # Conversion from 'type1' to 'type_2' is sign-extended. This may
        # cause unexpected runtime behavior.
        /w14905 # wide string literal cast to 'LPSTR'
        /w14906 # string literal cast to 'LPWSTR'
        /w14928 # illegal copy-initialization; more than one user-defined
        # conversion has been implicitly applied
        )

set(CLANG_WARNINGS
        -Wall
        -Wextra # reasonable and standard
        -Wshadow # warn the user if a variable declaration shadows one from a
        # parent context
        -Wnon-virtual-dtor # warn the user if a class with virtual functions has a
        # non-virtual destructor. This helps catch hard to
        # track down memory errors
        #-Wold-style-cast # warn for c-style casts
        -Wcast-align # warn for potential performance problem casts
        -Wunused # warn on anything being unused
        -Woverloaded-virtual # warn if you overload (not override) a virtual
        # function
        -Wpedantic # warn if non-standard C++ is used
        -Wconversion # warn on type conversions that may lose data
        -Wsign-conversion # warn on sign conversions
        -Wnull-dereference # warn if a null dereference is detected
        -Wdouble-promotion # warn if float is implicit promoted to double
        -Wformat=2 # warn on security issues around functions that format output
        # (ie printf)
        )

set(GCC_WARNINGS
        ${CLANG_WARNINGS}
        -Wmisleading-indentation # warn if indentation implies blocks where blocks
        # do not exist
        -Wduplicated-cond # warn if if / else chain has duplicated conditions
        -Wduplicated-branches # warn if if / else branches have duplicated code
        -Wlogical-op # warn about logical operations being used where bitwise were
        # probably wanted
        -Wuseless-cast # warn if you perform a cast to the same type
        )

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    message("Building with MSVC warnings")
    set(PROJECT_WARNINGS ${MSVC_WARNINGS})
elseif (CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    message("Building with Clang warnings")
    set(PROJECT_WARNINGS ${CLANG_WARNINGS})
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("Building with GCC warnings")
    set(PROJECT_WARNINGS ${GCC_WARNINGS})
else()
    message("Could not determine compiler for warnings")
endif ()

set(EXECUTABLE_FIELD metric_field)
set(LIB_FIELD metric_field_lib)
set(INPUT_FIELD input_handler)

add_executable(${INPUT_FIELD}
        src/input_handler.cpp)

add_executable(${EXECUTABLE_FIELD}
        src/metric_field.cpp
        src/other_fields.cpp
        src/predicates.c
        src/tests.cpp
        src/main.cpp
        )


set(SOURCES
        src/wrapper.cpp
        src/metric_field.cpp
        src/other_fields.cpp
        src/predicates.c)

add_library(${LIB_FIELD} STATIC ${SOURCES} )

target_compile_options(${EXECUTABLE_FIELD}
        PRIVATE ${PROJECT_WARNINGS})

# for shared lib
target_compile_options(${LIB_FIELD}
        PUBLIC ${PROJECT_WARNINGS})

target_compile_options(${INPUT_FIELD}
        PRIVATE ${PROJECT_WARNINGS})

target_include_directories(${EXECUTABLE_FIELD}
        PRIVATE include)

target_include_directories(${INPUT_FIELD}
        PRIVATE include)

# for shared lib
target_include_directories(${LIB_FIELD}
        PUBLIC include)

target_link_libraries(${EXECUTABLE_FIELD}
        PRIVATE Eigen3::Eigen
        PRIVATE OpenVolumeMesh::OpenVolumeMesh
        # PRIVATE GaussianProcess
        PRIVATE CGAL::CGAL)

target_link_libraries(${INPUT_FIELD}
        PRIVATE Eigen3::Eigen
        PRIVATE OpenVolumeMesh::OpenVolumeMesh
        PRIVATE CGAL::CGAL)

# for shared lib
target_link_libraries(${LIB_FIELD}
        PUBLIC Eigen3::Eigen
        PUBLIC OpenVolumeMesh::OpenVolumeMesh
        # PUBLIC GaussianProcess
        PUBLIC CGAL::CGAL)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Using libc++")
    target_compile_options(${EXECUTABLE_FIELD}
            PRIVATE -stdlib=libc++)
    target_link_libraries(${EXECUTABLE_FIELD}
            PRIVATE c++
            PRIVATE c++abi)
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Using libc++")
    target_compile_options(${LIB_FIELD}
            PUBLIC -stdlib=libc++)
    target_link_libraries(${LIB_FIELD}
            PUBLIC c++
            PUBLIC c++abi)
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Using libc++")
    target_compile_options(${INPUT_FIELD}
            PRIVATE -stdlib=libc++)
    target_link_libraries(${INPUT_FIELD}
            PRIVATE c++
            PRIVATE c++abi)
endif ()

# ########################################################