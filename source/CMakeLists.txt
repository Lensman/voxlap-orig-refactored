# -------------------------------------
# Voxlap5 cmake config
# -------------------------------------

cmake_minimum_required (VERSION 2.6)
project(Voxlap)

set (Voxlap_VERSION_MAJOR 5)
set (Voxlap_VERSION_MINOR 2)

# Currently supports win or sdl

set ( GFX sdl )

# -------------------------------------
# These will propagate to voxlap5.h 
# -------------------------------------

set (USEPHYSICS 1)
set (USEV5ASM 1)
set (UZEKZ 1)
set (USESCREEN 1)
set (FILEIO 1)

set (Voxlap_VSID 1024)
set (Voxlap_MAXXDIM 1024)
set (Voxlap_MAXYDIM 768)
set (Voxlap_MAXZDIM 256)

# -------------------------------------
# Set global paths (used by subfolders)
# -------------------------------------

set (LIBRARY_PATH "../libraries" )
set (CMAKE_ASM_CREATE_STATIC_LIBRARY 1)
set (INCLUDE_DIRECTORY "../include" )

# -------------------------------------
# Locate packages (SDL / DDRAW)
# -------------------------------------

if ( GFX STREQUAL "sdl" )
  Find_Package ( SDL REQUIRED )
elseif( GFX STREQUAL "win" )
  include(FindDirectInput.cmake)
endif()

exec_program(
    "git"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "describe --all --long"
    OUTPUT_VARIABLE VERSION 
)

string( REGEX MATCH "-g.*$" VERSION_SHA1 ${VERSION} )
string( REGEX REPLACE "[-g]" "" VERSION_SHA1 ${VERSION_SHA1} )

set (Voxlap_GIT_HASH "${VERSION_SHA1}")

# -------------------------------------
# Add git verion number to voxlap5.h
# -------------------------------------
configure_file (
  "../include/voxlap5.h.in"
  "../include/voxlap5.h"
)

if(MSVC)

  # -------------------------------------
  # Get correct assembler definitions
  # -------------------------------------

  INCLUDE(CMakeDetermineMASMDialect.cmake)
  INCLUDE(CMakeDetermineASMCompiler.cmake)
  add_definitions(/Fx /Gy /GF /Gs /USEKZ /MD /W0 /Ob2 /Oi /Fdbinaries)

  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:SSE2 /fp:fast")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /fp:fast")
else()
  set (CMAKE_C_FLAGS ${CMAKE_C_FLAGS "-funsigned-char -m32 -mfpmath=sse -msse2 -m3dnow -ffast-math")
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-msse2")
endif()


# -------------------------------------
# Set definitions
# -------------------------------------

IF(USEV5ASM)
  add_definitions(-DUSEV5ASM)
ENDIF(USEV5ASM)

add_definitions(-DGS -DUSEKZ=${USEKZ} -DFILEIO=${FILEIO} -DUSEPHYSICS=${USEPHYSICS} -DUSESCREEN=${USESCREEN} -DMODELLING=1 )

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../include)

# -------------------------------------
# Build Modules
# -------------------------------------

add_subdirectory (modules)

# -----------------------------------
# Build executables
# -----------------------------------

add_subdirectory (apps)
