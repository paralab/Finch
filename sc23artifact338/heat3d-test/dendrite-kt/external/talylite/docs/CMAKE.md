Q. Why did we switch to CMake?

A. CMake is the standard for cross-platform C++ projects. CMake projects contain a file named `CMakeLists.txt` that defines how the project is built, similar to a [Makefile](https://en.wikipedia.org/wiki/Makefile). The major selling point of CMake is that CMake can use this file to generate [many other types of project files](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html). This include project files for [IDEs](https://en.wikipedia.org/wiki/Integrated_development_environment) like Visual Studio, Xcode, and [CLion](https://www.jetbrains.com/clion/), as well as command-line build systems like make. "Out-of-source builds" are another useful feature in CMake which keeps all the files generated during the build process in a different directory than the source code. This helps keep the source code directory nice and clean.

Also, it generates a Makefile with pretty colors.

---

Q. How do I install CMake?

A. It may already be installed - try `cmake --version` to check. TalyFEM requires version 3.5 or above. If you don't have it, `sudo apt-get install cmake3` or `sudo apt-get install cmake`. If that fails, you can download and install it manually [from the CMake website](https://cmake.org/download/).

---

Q. How do I use CMake?

A. From inside a CMake project folder (which contains a `CMakeLists.txt` file):

```bash
mkdir build
cd build
cmake ..
make -j4
```

You can pass options to CMake by adding `-D[OPTION]=[VALUE]` to your CMake command. For TalyFEM projects, you need to add `-Dtalyfem_DIR=/path/to/talyfem/build`.

---

Q. Do I need to re-run CMake every time I make a change to my code?

A. No! You only need to re-run CMake if your CMakeLists.txt file changes. Even then, the Makefile that CMake generates is programmed to automatically re-run CMake if it detects a CMakeLists.txt change.

---

Q. I made a copy of my project and edited the code, but when I run `make` it keeps building the old version!

A. The files in the build directory that CMake generates are not relative. They are hard-coded to use the absolute path to the project they were created for. You must regenerate the build folder from scratch.

```bash
rm -rf build/* && cd build && cmake .. -Dtalyfem_DIR=/path/to/talyfem/build && make
```

---

Q. I keep getting configuration errors! How do I debug CMake problems?

A. If you're on a desktop or have X forwarding set up, try using `cmake-gui`. It works like the `cmake` command, but it opens up a nice user interface that lets you review available CMake options instead of immediately generating your project files.

---

Q. CMake is configuring OK, but when I run `make` I get odd errors. How can I debug?

A. If you think CMake has generated bad build files, try `make VERBOSE=1`. This will print out the compiler and linker commands `make` is invoking as it runs.

---

Q. How do I control what compiler CMake uses?

A. Add something like `-DCMAKE_CXX_COMPILER=mpicxx` to your CMake command line arguments.

---

Q. How do I control compiler flags when using CMake? How do I turn on optimizations?

A. Typically the reason you want to fiddle with compilation flags is to toggle optimizations and debugging information. CMake has a concept of different "build types" that tries to automatically manage this for you across different compilers (where the command-line options may change).

* `-DCMAKE_BUILD_TYPE=Release`: `-O3 -DNDEBUG` (max optimizations, asserts disabled, no debug info)
* `-DCMAKE_BUILD_TYPE=RelWithDebInfo`: `-O2 -g` (good optimizations, asserts enabled, debug info)
* `-DCMAKE_BUILD_TYPE=Debug`: `-O0 -g` (no optimization, asserts enabled, debug info)

You should develop with the `Debug` build type, test at scale with `RelWithDebInfo`, and do production runs with `Release`. The default build type (if none is specified) is `Debug`.

---

Q. No, seriously, I know what I'm doing, how do I control compiler flags?

A. `-DCMAKE_CXX_FLAGS="-g -O3"` will add custom flags, regardless of build type.

---

Q. I changed some arguments or environment variables, but CMake seems to be ignoring them!

A. CMake caches a lot of variables to make re-running CMake fast - but sometimes it is too aggressive. Try clearing out your build directory (`rm -rf build/*`) and re-running CMake.

---

Q. How does CMake find where PETSc is installed?

A. It looks at the `$PETSC_DIR` and `$PETSC_ARCH` environmental variables - you should set both before running CMake.

---

Q. How does CMake find where libconfig is installed?

A. First, it checks the "usual" system directories (like `/usr`). This means CMake will automatically pick up the install directory if you installed libconfig with the package manager or used `make install` with the default install prefix. If that fails, CMake will check the `$LIBCONFIG_DIR` environmental variable.

---

Q. After I compile TalyFEM, do I still need to set `PETSC_DIR`, `PETSC_ARCH`, and `LIBCONFIG_DIR` when I'm compiling my own project?

A. No. The PETSc and libconfig directories are baked into the `talyfemConfig.cmake` file that CMake generated when you configured TalyFEM. These environment variables will be ignored (as long as you don't try to re-configure TalyFEM).

---

Q. I have an old TalyFEM project that uses a Makefile. How do I upgrade to this new CMake version of TalyFEM?

1. Back up your project.

2. Delete your `makefile` (or `Makefile`).

3. Create a new file named `CMakeLists.txt` in the root of your project directory. Put this in it:

```cmake
cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

# This tells CMake what the name of our project is (optional).
project(MyProject)

# This declares an executable named 'hello' and lists the .cpp files that need
# to be compiled to make the executable.
add_executable(hello
  # list source files to compile
  src/main.cpp
  src/another_file.cpp
)

# Next, we add the ./include directory to the include search path for the program.
# This is what lets us "#include" any header files we put in the "include/" folder.
target_include_directories(hello PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Now we tell CMake to find TalyFEM.
find_package(talyfem REQUIRED)

# Finally, we link our excecutable defined above with TalyFEM. This will
# also link our program with TalyFEM's dependencies (PETSc, libconfig).
target_link_libraries(hello talyfem)
```

4. Adjust the list of .cpp files in "add_executable" according to your project.

5. Run `mkdir build && cd build`

6. Make sure you have CMake 3.5 or above installed. You can use `cmake --version` to check.

7. Run `cmake .. -Dtalyfem_DIR=/path/to/talyfem/build` (Don't forget the `/build`! This folder should contain a file named `talyfemConfig.cmake` that was  generated by CMake during configuration.)

8. Run `make` to compile your code (or `make -j4` to compile it in parallel).

---

Q. How do I generate Doxygen documentation for my project?

A. TalyFEM includes some CMake scripts to make this easy. In your CMakeLists.txt file, somewhere after `find_package(talyfem REQUIRED)`, add:

```cmake
include(TalyDoxygen)
add_doxygen(TARGET myproject INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
```

(update myproject accordingly)

Then, type `make docs` to build your documentation. The main page should be generated at `docs/html/index.html`.

---

Q. How do I add reproducibility information to my project?

A. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md).