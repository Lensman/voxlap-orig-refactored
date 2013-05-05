# -------------------------------------
# Voxlap5 cmake config
# -------------------------------------

cmake_minimum_required (VERSION 2.6)
project(Voxlap)

set (Voxlap_VERSION_MAJOR 5)
set (Voxlap_VERSION_MINOR 2)

# Currently supports win or sdl

set ( GFX sdl )

# -------------------------------------------
# These will propagate to voxlap5.h / v5.?asm
# -------------------------------------------

set (USEPHYSICS 1)
set (USEV5ASM 1)
set (UZEKZ 1)
set (USESCREEN 1)
set (FILEIO 1)
set (USEZBUFFER 1)

set (Voxlap_VSID 1024)

# vsid 1024 = prec = 4096 2048 = 2048
# for mipmap in v5.?asm log2(VSID) 512 = 9, 1024 = 10, 2048 = 11

if ( Voxlap_VSID EQUAL 4096 )
  set (Voxlap_LVSID 12)     
  set (Voxlap_PREC_TABLE 2048)  
elseif ( Voxlap_VSID EQUAL 2048 )
  set (Voxlap_LVSID 11)     
  set (Voxlap_PREC_TABLE 2048)  
elseif( Voxlap_VSID EQUAL 1024 )
  set (Voxlap_LVSID 10)
  set (Voxlap_PREC_TABLE 4096)
elseif( Voxlap_VSID EQUAL 512 )
  set (Voxlap_LVSID 9)
  set (Voxlap_PREC_TABLE 4096)
endif()

set (Voxlap_USEZBUFFER ${USEZBUFFER} ) 
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

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/v5.asm.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/v5.asm"
)

if(MSVC)

  # -------------------------------------
  # Get correct assembler definitions
  # -------------------------------------

  INCLUDE(CMakeDetermineMASMDialect.cmake)
  INCLUDE(CMakeDetermineASMCompiler.cmake)
  add_definitions(/Fx /Gy /GF /Gs /USEKZ /MD /Ob2 /Oi /Fdbinaries)

  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:SSE2 /fp:fast")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /fp:fast")
else()
  set (CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-funsigned-char -m32 -mfpmath=sse -msse2 -m3dnow -ffast-math")
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
