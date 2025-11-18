cmake_minimum_required(VERSION 3.12)
project(nrPlotter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include FetchContent for managing all dependencies
include(FetchContent)


message(STATUS "Setting up CLI11 via FetchContent...")
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.5.0  # Nowsza wersja: wymaga CMake 3.9, brak deprecjacji z <3.5
)
FetchContent_MakeAvailable(cli11)


# --- STB_IMAGE Dependency (Image Loading) ---
message(STATUS "Fetching stb_image.h...")
FetchContent_Declare(
    stb
    GIT_REPOSITORY https://github.com/nothings/stb.git
    GIT_TAG        master # A stable header-only library
)
FetchContent_GetProperties(stb)
if(NOT stb_POPULATED)
    FetchContent_MakeAvailable(stb)
endif()

# --- GLAD Dependency (OpenGL Loader) ---
add_library(glad STATIC external/glad/src/glad.c)
target_include_directories(glad PUBLIC external/glad/include)

# --- GLFW Dependency ---
find_package(glfw3 3.3 QUIET)
if(glfw3_FOUND)
  message(STATUS "Found system-installed glfw3.")
  set(GLFW_TARGET glfw3::glfw)
else()
  message(STATUS "glfw3 not found. Fetching from GitHub...")
  FetchContent_Declare(
      glfw
      GIT_REPOSITORY https://github.com/glfw/glfw.git
      GIT_TAG        3.3.8
  )
  FetchContent_MakeAvailable(glfw)
  set(GLFW_TARGET glfw)
endif()

# --- pbPlots Dependency ---
# This is the new, automated way to get pbPlots.
message(STATUS "Fetching pbPlots...")
FetchContent_Declare(
    pbplots_fetch # Renamed to avoid conflict with the library target name
    GIT_REPOSITORY https://github.com/InductiveComputerScience/pbPlots.git
    # Using 'main' branch. For production, a specific GIT_TAG (commit hash) is safer.
    GIT_TAG        master
)
FetchContent_GetProperties(pbplots_fetch)
if(NOT pbplots_fetch_POPULATED)
    FetchContent_MakeAvailable(pbplots_fetch)
endif()

# Since pbPlots has no CMakeLists.txt, we create our own library target for it.
add_library(pbplots STATIC
    ${pbplots_fetch_SOURCE_DIR}/Cpp/pbPlots.cpp
    ${pbplots_fetch_SOURCE_DIR}/Cpp/supportLib.cpp
)
# Make the pbplots headers available to any target that links it.
target_include_directories(pbplots PUBLIC
    ${pbplots_fetch_SOURCE_DIR}/Cpp
)

# --- Main Executable ---
add_executable(main
    src/main.cpp
    src/plotter.cpp
    src/plotter.h
	src/overlay.cpp
	src/overlay.h
)

# Add include directories for ALL our header-only/fetched dependencies
target_include_directories(main PRIVATE
    ${stb_SOURCE_DIR}
    ${glad_SOURCE_DIR}/include
	${cli11_SOURCE_DIR}/include
)

# Link the executable against all our library targets
target_link_libraries(main PRIVATE
    glad
    ${GLFW_TARGET}
    pbplots
)

# --- Platform specific link libraries ---
if(UNIX AND NOT APPLE AND NOT glfw3_FOUND)
    find_package(Threads REQUIRED)
    target_link_libraries(main PRIVATE Threads::Threads m dl)
endif()

if(UNIX AND NOT APPLE)
    target_link_libraries(main PRIVATE dl)
endif()