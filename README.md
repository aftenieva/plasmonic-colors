# plasmonic-colors
Color calcuator for plasmonic gold nanoparticles
=========

`spectrum_to_rgb` is a pure Python module to calculate the RGB colors of nanoparticles, percieved by the observer.
The absorption, scattering and extinction spectra of a gold sphere in various environments are calculated with the help of Mie-theory-based functions from the 'miepython' package, developed by Scott Pharl (https://github.com/scottprahl/miepython).
After providing the sphere radius in nanometers, and selecting the surrounding environment (solvent) as a manual input, the output spectra will be displayed on the screen.
The corresponding scattered, absorbed and transmitted RGB colors on the scale from 0 to 255 (the latter being the complementary to the color, calculated from extinction spectra) will be printed in the console.

The described code uses the data from the accompanying files 'n k.txt' and 'x y z D65' that contain the optical properties of gold (real and imaginary part of the refractive index) and CIE color matching functions together with the standard illuminant spectrum, respectively.
In principle, materials other than gold can be used, as well as other illuminants or color matching functions, with an only requirement of the matching wavelength span. Here, the whole visible spectrum (360 - 830 nm) is covered with a step of 1 nm. 


License
-------

`spectrum_to_rgb` is licensed under the terms of the MIT license.
