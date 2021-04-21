Computer-aided calculation of PKF dynamics with Maxima (CAC-PKF-M v0.1)
=======================================================================


O. Pannekoucke (olivier.pannekoucke@meteo.fr)

date:
 * 2017 developpment of the v0.1
 * 2021 v0.1 of the code is made available in open access

CeCILL-B FREE SOFTWARE LICENSE AGREEMENT

Introduction
------------

This package is a part of the preliminary symbolic computation in [Maxima][], of the parametric Kalman filter (PKF) dynamics applied to the Burgers' equation detailed in [Pannekoucke et al. (2018)][]. The motivation is to provide a first attempt to validate the calculation done by hand, and to think about a systematic derivation of the PKF equations. This notebook relies on the compagnion library `./lib/libpkf.mc` which implements the expectation operator used to facilitate the calculations. 

The code has been developped in 2017, and runs with the 5.43.2 version of Maxima (Jan. 2020).

A Python implementation of the symbolic calculation of the PKF dynamics, [SymPKF][] is also available (see [Pannekoucke and Arbogast (2021)][]).

Details of the repository
-------------------------

 * `pkf-burgers.wxm` is the [wxMaxima][] notebook used to the calculation of the diffusion part of the PKF applied to the Burgers' equation. 
 * `pkf-burgers.html` is the HTML saving of the output of the notebook [pkf-burgers.wxm](./pkf-burgers.html)
 * `./lib/libpkf.mc` contains the library which implements the expectation operator used for the PKF computation.


References
----------

Maxima, a Computer Algebra System. Version 5.43.2 : https://maxima.sourceforge.io/

O. Pannekoucke, M. Bocquet, and R. Ménard, “Parametric covariance dynamics for the nonlinear diffusive Burgers’ equation,” Nonlinear Processes in Geophysics, vol. 2018, pp. 1–21, 2018, doi: https://doi.org/10.5194/npg-2018-10.

O. Pannekoucke and P. Arbogast, “SymPKF: a symbolic and computational toolbox for the design of univariate parametric Kalman filter dynamics. https://arxiv.org/abs/2103.09226

wxMaxima : a document based interface for the computer algebra system Maxima. https://wxmaxima-developers.github.io/wxmaxima/

[Maxima]: https://maxima.sourceforge.io/ "Maxima, a Computer Algebra System. Version 5.43.2"
[wxMaxima]: https://wxmaxima-developers.github.io/wxmaxima/ "wxMaxima : a document based interface for the computer algebra system Maxima."
[Pannekoucke et al. (2018)]: https://doi.org/10.5194/npg-2018-10 'O. Pannekoucke, M. Bocquet, and R. Ménard, “Parametric covariance dynamics for the nonlinear diffusive Burgers’ equation,” Nonlinear Processes in Geophysics, vol. 2018, pp. 1–21, 2018, doi: https://doi.org/10.5194/npg-2018-10.'
[SymPKF]: https://github.com/opannekoucke/sympkf 'O. Pannekoucke, “SymPKF: a symbolic and computational toolbox for the design of parametric Kalman filter dynamics.” Zenodo, 2021, doi: 10.5281/ZENODO.4608514.'
[Pannekoucke and Arbogast (2021)]: https://arxiv.org/abs/2103.09226 'O. Pannekoucke and P. Arbogast, “SymPKF: a symbolic and computational toolbox for the design of univariate parametric Kalman filter dynamics. https://arxiv.org/abs/2103.09226'