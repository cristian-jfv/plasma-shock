cmake_minimum_required(VERSION 3.25)

project(numerical_exercise Fortran)

enable_language(Fortran)

set(GENERAL_SOURCES bookkeeping.f90 model.f90 mesh.f90 scheme.f90)
add_library(general_modules OBJECT ${GENERAL_SOURCES})
target_compile_options(general_modules PRIVATE -fdefault-real-8 -fopenmp -fimplicit-none)

# First Order Arithmetic Mean Matrix
set(FOAMM_SOURCES foamm_main.f90 foamm_scheme.f90)
add_executable(foamm_solver $<TARGET_OBJECTS:general_modules> ${FOAMM_SOURCES})
target_compile_options(foamm_solver PRIVATE -fdefault-real-8 -fopenmp -fimplicit-none)

# Lax-Friedrich
set(LXF_SOURCES lxf_main.f90 lxf_scheme.f90)
add_executable(lxf_solver $<TARGET_OBJECTS:general_modules> ${LXF_SOURCES})
target_compile_options(lxf_solver PRIVATE -fdefault-real-8 -fopenmp -fimplicit-none)

# TVD Lax-Friedrich
set(TVDLXF_SOURCES tvdlxf_main.f90 tvdlxf_scheme.f90)
add_executable(tvdlxf_solver $<TARGET_OBJECTS:general_modules> ${TVDLXF_SOURCES})
target_compile_options(tvdlxf_solver PRIVATE -fdefault-real-8 -fopenmp -fimplicit-none)

set(CMAKE_VERBOSE_MAKEFILE ON)
