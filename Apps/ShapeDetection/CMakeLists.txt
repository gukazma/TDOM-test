set(TARGET_NAME ShapeDetection)


file(GLOB_RECURSE CPP_HEAD_FILES  ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB_RECURSE CPP_SOUCE_FILES  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

source_group(TREE ${CMAKE_CURRENT_LIST_DIR} FILES ${CPP_HEAD_FILES} ${CPP_SOUCE_FILES})

add_executable(${TARGET_NAME} ${CPP_HEAD_FILES} ${CPP_SOUCE_FILES})
set_target_properties(${TARGET_NAME} PROPERTIES FOLDER Apps)
find_package(tinyobjloader CONFIG REQUIRED)
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${TARGET_NAME} PRIVATE
CGAL::CGAL
Eigen3::Eigen
tinyobjloader::tinyobjloader
)