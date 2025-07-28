nrPlotter - numerica recipes plotter

![alt text](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)

![alt text](https://img.shields.io/badge/Platform-Linux%20%7C%20Windows%20%7C%20macOS-lightgrey.svg)

![alt text](https://img.shields.io/badge/Build-CMake-green.svg)

nrPlotter is a cross-platform C++ application for visualizing data, designed as a template for scientific computing and numerical recipes. It uses a modern CMake build system that automatically manages all dependencies, allowing for a seamless build experience.

The application can generate plots from hardcoded data, user-defined functions, or interactive user input, and allows for detailed configuration via command-line options.
Dependencies

This project uses CMake's FetchContent module to download and configure all required libraries automatically. You do not need to install any of them manually.

    CMake (>= 3.12): The build system generator used to configure and build the project.

    GLFW: A multi-platform library used for creating the window, handling the OpenGL context, and receiving keyboard/mouse input.

    GLAD: An OpenGL Loading Library that generates the code needed to use modern OpenGL functions.

    pbPlots: A simple C++ library used to generate the plot data and render it to a PNG image.

    stb_image: A single-header C library for loading the generated PNG image as a texture in OpenGL.

    CLI11: A powerful, modern, and header-only C++ library used for parsing all command-line arguments.

Installation and Build

Ensure you have the necessary development tools installed on your system before proceeding.
Prerequisites

    Git: For cloning the repository.

    CMake (version 3.12 or newer): For configuring the build.

    A C++17 compliant compiler:

        Linux: g++ or clang++ (install with sudo apt-get install build-essential).

        macOS: Apple Clang (install with xcode-select --install).

        Windows: Visual Studio 2019 or newer (with the "Desktop development with C++" workload).

Build Steps

1. Clone the Repository

Open a terminal or command prompt and clone the project.
Generated bash

      
git clone https://github.com/SWENG-UKSW/nr_protter.git
cd nr_protter


2. Configure the Build with CMake

Create a build directory and run CMake from within it. This step will automatically fetch all dependencies.
Generated bash

      
mkdir build
cd build
cmake ..

**Note:** The first time you run `cmake ..`, it may take a moment as it downloads the required libraries. Subsequent runs will be instantaneous.

**3. Compile the Project**

Once configuration is complete, build the executable.

cmake --build .


The compiled executable (nrPlotter or nrPlotter.exe) will be located inside the build directory.
Features

The application supports several modes of interaction for generating and manipulating plots.
Key	Action
a	Plot from Function: Generates a new plot from a pre-defined mathematical function (e.g., sine wave).
s	Plot Hardcoded Data: Generates a plot from a hardcoded set of (x,y) data points.
d	Plot from Clicks: Generates a new plot using the coordinates of all the points the user has added by clicking on the window.
x	Collect X-Range: Enters a mode where the next two mouse clicks define the minimum and maximum X-axis range for future plots.
y	Collect Y-Range: Enters a mode where the next two mouse clicks define the minimum and maximum Y-axis range for future plots.
c	Clear User Input: Clears all stored mouse clicks and resets the custom X/Y ranges.
ESC	Exit: Closes the application.
Command-Line Options

The plot dimensions and viewing ranges can be configured at launch using the following options.
Option	Long Option	Description
-x	--plot-width	Sets the plot width in pixels.
-y	--plot-height	Sets the plot height in pixels.
-a	--pad-width	Sets the plot padding width in pixels.
-s	--pad-height	Sets the plot padding height in pixels.
-z	--range-minx	Sets the minimum value of the X-axis.
-c	--range-maxx	Sets the maximum value of the X-axis.
-t	--range-miny	Sets the minimum value of the Y-axis.
-u	--range-maxy	Sets the maximum value of the Y-axis.
-h	--help	Displays the help message with all available options.
Example Usage

To run the application and generate a plot of size 1280x720 with a specific X-axis range:

./nrPlotter --plot-width 1280 --plot-height 720 --range-minx -10.0 --range-maxx 10.0

    