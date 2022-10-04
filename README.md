# SPECTRUM
A C library for spectral estimation and signal processing.

The main goal of SPECTRUM is to provide a C library that can be easily used
in scientific programming to perform spectral estimation and other operations
such as spectrum sensing.

I started developing this library while working on a project funded by the
[European Space Agency (ESA)](https://www.esa.int/) at
[Dapcom Data Services](https://www.dapcom.es/).

## Features
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

If you are missing some features first make sure that
they are not listed below. Otherwise check how can
you contribute.

## Work in progress
Here there is a list of some features that I know that are missing:  
* [ ] Support for real FFT (input data type must be checked).
* [ ] Padding support to avoid discarding the last frame if it it smaller.
* [ ] Support for Cross Spectral Density (CSD) estimation.

## Contributing
If you find a bug or you would like to suggest a new feature, feel free to open an issue.

Pull requests are also welcome but please create an issue first if it's a new feature.

## Acknowledgements
I would like to thank the European Space Agency (ESA) and Dapcom Data Services for the funding that allowed
the development of this library. I would also like to
thank Jordi Portell for his contributions in the code.
