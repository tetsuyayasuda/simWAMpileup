cmake_minimum_required(VERSION 2.8)

set(CMAKE_INSTALL_PREFIX ./)
set(app simWAMpileup)

### SFITS ###
set(SFITS_LIB sfitsio sllib z bz2 readline curses)

### ROOT ###
set(ROOTSYS $ENV{ROOTSYS})
set(ROOTDIR "")
set(ROOT_INC_DIR ${ROOTSYS}/include${ROOTDIR})
set(ROOT_LIB_DIR ${ROOTSYS}/lib${ROOTDIR})
set(ROOT_LIB
 Core RIO Net Hist Graf Graf3d Gpad Tree
 Rint Postscript Matrix Physics MathCore Thread)

include_directories(include ${ROOT_INC_DIR})
link_directories(${ROOT_LIB_DIR})

add_executable(${app}
  ${app}.cc
  src/powCalc.cc
  src/simPileup.cc
  src/readText.cc
  src/readRespWAM.cc
  src/readSpecWAM.cc
)

target_link_libraries(${app} ${SFITS_LIB} ${ROOT_LIB}
)
