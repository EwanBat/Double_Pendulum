cmake_minimum_required(VERSION 3.10)
project(DoublePendulum VERSION 1.0)

# Standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Recherche des packages (décommentez selon vos besoins)
# find_package(OpenGL REQUIRED)
# find_package(GLFW3 REQUIRED)
# find_package(SDL2 REQUIRED)

# Fichiers sources (ajustez selon vos fichiers)
set(SOURCES
    src/main.cpp
    src/double_pendulum.cpp
    src/simple_pendulum.cpp
)

# Fichiers headers (optionnel)
set(HEADERS
    include/double_pendulum.hpp
    include/simple_pendulum.hpp
)

# Répertoires d'inclusion
include_directories(include)

# Création de l'exécutable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Liaison avec les bibliothèques (décommentez selon vos besoins)
# target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL)
# target_link_libraries(${PROJECT_NAME} PRIVATE glfw)
# target_link_libraries(${PROJECT_NAME} PRIVATE SDL2::SDL2)

# Options de compilation
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()
