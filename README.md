# Spectra
A C library for spectral estimation and signal processing.

The main goal of Spectra is to provide a C library that can be easily used
in scientific programming to perform spectral estimation and other operations
such as spectrum sensing.

The spectral estimation functions try to imitate the behavior of SciPy, hence
the source code structure of this part is very similar to that of SciPy.

I started developing this library while working on a project funded by the
[European Space Agency (ESA)](https://www.esa.int/) at
[Dapcom Data Services](https://www.dapcom.es/) (ESA Contract No. 4000137290) .

## Project information
### Features
Spectral estimation methods:  
* Periodogram.
* Averaged periodogram (Welch's method).

Supported window types:  
* Rectangular.
* Hann.
* Hamming.
* Blackman.

Other features:  
* Noise power estimation.
* Energy detector.

If you are missing some features first make sure that they are not listed below.
Otherwise check how can you contribute.

### Work in progress
Here there is a list of some features that I know that are missing:  
* [ ] Support for real FFT (input data type must be checked).
* [ ] Padding support to avoid discarding the last frame if it it smaller.
* [ ] Support for Cross Spectral Density (CSD) estimation.
* [ ] Support for median averaging in the Welch's method.

### Documentation
The project documentation is managed with Doxygen. It can be generated by running
`doxygen` on the root of the repository. This will create the directory `docs/` with
all the documentation in HTML.

### Contributing
If you find a bug or you would like to suggest a new feature, feel free to open an issue.

Pull requests are also welcome but please create an issue first if it's a new feature.

### Acknowledgements
I would like to thank the [European Space Agency (ESA)](https://www.esa.int/) and
[Dapcom Data Services](https://www.dapcom.es/) for the funding that allowed
the development of this library. I would also like to thank Jordi Portell
for his contributions in the code.

### License
This project is licensed under the BSD 3-clause license.  
The KissFFT library is also licensed under the BSD 3-clause license.  
The GSL library (only used in the energy detector) is licensed under the GPL-3.

## Build
The compilation of this project is managed by CMake. You can either build the library or
the different examples provided.

First we generate the build system:
```
cmake -B build .
```

To build all the targets we can execute:
```
cmake --build build
```

To build only one target:
```
cmake --build build --target <target_name>
```

The available target names are:  
* `spectra`: builds the shared library.
* `examples`: builds all the examples (and also the library).
* `detect_signal_aic_ed`: builds only the example (and the library).

## Troubleshooting
### KISS FFT library not found
```
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.
Please set them or make sure they are set and tested correctly in the CMake files:
KISSFFT_LIB
    linked by target "spectra"
```
1. Make sure that the library is installed or available in the `lib` project directory.
2. CMake searches for `kissfft-float` and `kiss_fft_float` in the `lib` directory of the
project and also in the system library path. On UNIX, make sure that the files `libkissfft-float.so`
or `libkiss_fft_float.so` exist. They must not have the version number at the end.
If they do, you can create a symlink.
