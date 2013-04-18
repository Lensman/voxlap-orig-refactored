cmake_minimum_required (VERSION 2.6)
project(Voxlap)
set (INCLUDE_DIRECTORY "../../../include" )
# Have to use custom commands for the voxlap core 

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/v5.obj
  COMMAND ml /c /coff ${CMAKE_CURRENT_SOURCE_DIR}/v5.asm
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/kplib.obj
  COMMAND cl /c /J /TP ${CMAKE_CURRENT_SOURCE_DIR}/kplib.cpp /Ox /Ob2 /GFy /Gs /MD
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/voxlap5.obj
  COMMAND cl /c /J /TP ${CMAKE_CURRENT_SOURCE_DIR}/voxlap5.cpp /Ox /Ob2 /GFy /Gs /MD
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cpu_detect.obj
  COMMAND cl /c /J /TP ${CMAKE_CURRENT_SOURCE_DIR}/cpu_detect.cpp /Ox /Ob2 /GFy /Gs /MD
)

set(v5_object_file ${CMAKE_CURRENT_BINARY_DIR}/v5.obj)
set(voxlap5_object_file ${CMAKE_CURRENT_BINARY_DIR}/voxlap5.obj)
set(cpu_detect_object_file ${CMAKE_CURRENT_BINARY_DIR}/cpu_detect.obj)
set(kplib_object_file ${CMAKE_CURRENT_BINARY_DIR}/kplib.obj)

set_source_files_properties (
  ${v5_object_file}
  PROPERTIES
  EXTERNAL_OBJECT false # to say that "this is actually an object file,so it should not be compiled, only linked"
  GENERATED true        # to say that "it is OK that the obj-files do not exist before build time"
)

set(VOXLAP_CORE
  ${kplib_object_file}
  ${voxlap5_object_file}
  ${v5_object_file}
  ${cpu_detect_object_file}
  "${CMAKE_CURRENT_SOURCE_DIR}/../../include/sysmain.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../include/cpu_detect.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../include/kplib.h"
)


function(compile_executable APP_NAME SYSTEM )
  set( APPLICATION_FULL_NAME ${GFX}_${APP_NAME} )
  add_executable ( 
    ${APPLICATION_FULL_NAME}
    ${SYSTEM}
    "${CMAKE_CURRENT_SOURCE_DIR}/../${APP_NAME}.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../${GFX}main.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/sysmain.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/cpu_detect.h"
    #"${CMAKE_CURRENT_SOURCE_DIR}/../../include/kplib.h"

    ${VOXLAP_CORE} 
  )
  if ( GFX STREQUAL "sdl" )
    target_link_libraries(${APPLICATION_FULL_NAME} ${SDL_LIBRARY})
    target_link_libraries(${APPLICATION_FULL_NAME} ${SDLMAIN_LIBRARY})
  elseif ( GFX STREQUAL "win" )
    target_link_libraries(${APPLICATION_FULL_NAME} ${DIRECTINPUT_DINPUT_LIBRARY})
    target_link_libraries(${APPLICATION_FULL_NAME} ${DIRECTINPUT_DDRAW_LIBRARY})
    target_link_libraries(${APPLICATION_FULL_NAME} ${DIRECTINPUT_DXGUID_LIBRARY})
  endif()
    
  message(STATUS "making ${SYSTEM} exe: ${APP_NAME} using ${VOXLAP_CORE}" )
endfunction()

if ( NOT NO_EXECUTABLES )
  compile_executable(simple WIN32)
  compile_executable(voxed WIN32)
endif()

