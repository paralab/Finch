# What is it?

TalyFEM implements some build reproducibility measures by generating a header file containing build-time information using the `build_info_gen.py` Python script. This information can be written out into a libconfig-formatted file at run-time using the `Repro` class. Currently, the information recorded includes:

* The current time and date
* The build directory
* The hostname of the build machine and user who ran the build
* The values of all environment variables at the time of compilation
* The compiler ID and version
* The contents of CMakeCache.txt (contains CMake configuration information for this build)
* The git commit hash this build is based on
* A diff patch containing any uncommitted changes to git-tracked files

It also records some run-time information:

* Command-line arguments the program was run with
* The contents of config.txt
* The number of MPI processes
* Rank 0's hostname (for identifying which cluster the program was run on)

# How do I use it?

```cpp
int main(int argc, char** argv) {
  PetscInitialize(...);

  Repro r(argc, argv);
  r.write();

  // the rest of your program...
}
```

This will write a file named `repro.cfg` that contains information about how TalyFEM was built and how the program was run. This will only write information about how TalyFEM (the library) was built - it won't contain any information about how your project executable was built.

To write information about your project, you'll need to generate your own "build info" header and register it with the `Repro` class. To generate the header, add the following to your CMakeLists.txt:

```cmake
if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
  include(TalyReproducibility)
endif()

# generate reproducibility information header (build_info_mytarget.h)
add_taly_reproducibility(TARGET mytarget)
```

(if you copied the CMakeLists.txt from the library, you may already have this)

Then, in your code, add:

```cpp
#define BUILDINFO_MYTARGET_IMPL  // new!
#include <build_info_mytarget.h>  // new!

int main(int argc, char** argv) {
  PetscInitialize(...);

  Repro r(argc, argv);
  r.add_build_info<BuildInfo_mytarget>();  // new!
  r.write();
}
```

(adjust define, include, and add_build_info<> according to your target's name)

# Advanced Usage

You can also include input files (like meshes or initial conditions) in `repro.cfg`. Be warned that this will increase the size of `input.cfg` according to the size of the input files.

```cpp
Repro r(argc, argv);
r.add_input_file("myfile.plt");
r.write();
```

Files are encoded in base64, so you can include any type of file you like.

---

If you want to add your own custom information to the `repro.cfg` file, use the `Repro::add_handler(...)` method:

```cpp
  /**
   * Type for a repro "handler".
   * Handlers write to a setting in a libconfig file.
   */
  typedef std::function<void(libconfig::Setting&)> repro_func_t;

  /**
   * Register a new section for this Repro object.
   * @param name Name of section. Multiple sections can be specified with
   *             '.' - ex: 'runtime.hostname'. It is strongly recommended
   *             that all handlers are prefixed with 'runtime.*'.
   * @param type Type of libconfig setting this section is. Common types:
   *             libconfig::Setting::TypeString, TypeInt, TypeBoolean,
   *             TypeFloat, TypeList, TypeGroup, TypeArray
   *             See the libconfig documentation for more.
   * @param f    Function that will write data to the setting. Typically a
   *             C++ lambda function.
   *             (See: http://en.cppreference.com/w/cpp/language/lambda)
   */
  void add_handler(const std::string& name, libconfig::Setting::Type type,
                   const repro_func_t& f);

```

See `src/utils/reproducibility.cpp` for examples of handlers.