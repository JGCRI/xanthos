---
title: 'Xanthos 2.0: A Python package for iterative unfolding'
tags:
  - Python
  - global hydrology
  - human-Earth systems
  - modeling
authors:
  - name: Chris R. Vernon
    orcid: 0000-0003-2164-7789
    affiliation: "1"
    email: james.bourbeau@icecube.wisc.edu
  - name: Zigfried Hampel-Arias
    orcid: 0000-0003-0253-9117
    affiliation: "1, 2"
    email: zhampel@wipac.wisc.edu
affiliations:
 - name: University of Wisconsin - Madison, Wisconsin, USA
   index: 1
 - name: IIHE, Université Libre de Bruxelles, Bruxelles, Belgium
   index: 2
date: 15 May 2018
bibliography: paper.bib
---

# Summary

In an ideal world, experimentalists would have access to the perfect detector:
an apparatus that makes no error in measuring a desired quantity. However,
in real life scenarios, detectors have:

- Finite resolutions

- Characteristic biases that cannot be eliminated

- Less than 100% detection efficiencies

- Statistical uncertainties

By building a matrix that encodes a detector’s smearing of the desired true
quantity into the measured observable(s), a deconvolution can be performed
that provides an estimate of the true variable. This deconvolution process is
known as unfolding.


PyUnfold is an extensible framework for the unfolding of discrete probability
distributions via the iterative unfolding method described in [@agostini].
This method falls into the general class of inverse problems, and is especially powerful
when no analytic form of the inverse function is explicitly available, instead requiring
an estimate (e.g. via a finite amount of simulation) of the response function.
Given that measured data comprise a finite sample, PyUnfold also implements the uncertainty
contributions stated in [@adye2].


The unfolding method itself is data-agnostic, referring to the measurement process
as the smearing of a set of true causes into a set of detectable effects.
For example one could define as causes the true energy of a particle and the effects
the measured energy of that particle in a detector.
Another example might be a set of diseases (causes) and possible clinical symptoms (effects).
So long as it is possible to encode estimable resolutions and biases connecting causes to
effects in a binned response matrix, one can perform a deconvolution with PyUnfold.


The primary purpose of PyUnfold is to provide an unfolding toolkit for members of all
scientific disciplines in an easy-to-use package.
For example, unfolding methods are commonly used in the analysis pipeline of the
high-energy physics (HEP) community.
However, the deconvolution packages used in HEP maintain a strong dependence on the
ROOT data analysis framework [@root], which is almost exclusively used in HEP.
Instead, PyUnfold is built on top of the Python scientific computing stack (i.e. NumPy,
SciPy, and pandas), thus broadening its scope to a general scientific audience.


PyUnfold has been designed to be both easy to use for first-time users as well as
flexible enough for fine-tuning an analysis and seamlessly testing the robustness
of results. It provides support for the following:

- Custom, user defined initial prior probability distributions, the default being
the uniform prior. The non-informative Jeffreys' prior [@jeffreys] is accessible
as a utility function.

- Unfolding stopping criteria based on test statistic calculations comparing unfolded
distributions from one iteration to the next. These include
Kolmogorov-Smirnov [@kolmogorov][@smirnov], $\chi^2$, relative difference,
and Bayes factor [@pfendner] tests.

- Tunable spline regularization as a means of ensuring that unfolded distributions do not
suffer from growing fluctuations potentially arising from the finite binning of the
response matrix.

- Option to choose between Poisson or multinomial forms of the covariance matrices
for both the data and response contributions to the uncertainty calculation.

- Multivariate unfolding via definitions of subsets of causes, which are regularized
in their respective blocks or groups.


Further mathematical details regarding the iterative unfolding procedure, including complete
derivations of the statistical and systematic uncertainty propagation can be found in the
online documentation.


PyUnfold has been applied successfully for the measurement of the cosmic-ray energy spectrum
by the HAWC Observatory [@hawc-crspectrum] and is currently being used in an analysis by
members of the IceCube Neutrino Observatory.


# Acknowledgements

The authors acknowledge support from the Wisconsin IceCube Particle Astrophysics Center
at the UW-Madison, and especially for the guidance of Professor Stefan Westerhoff.
We also acknowledge the financial support provided by the National Science Foundation,
the Belgian American Educational Foundational Fellowship, and Wallonie-Bruxelles International.

# References