cmake_minimum_required(VERSION 3.16)
project(hitmap)

find_package(ROOT 6.26 CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/source/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

find_package(yaml-cpp)
include_directories(${YAML_CPP_INCLUDE_DIRS})

#add_executable(hitmap main.cc ${headers} ${sources})
add_executable(hitmap hitmap.cc ${headers} ${sources})
target_link_libraries(hitmap PUBLIC ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS} ${YAML_CPP_LIBRARIES})