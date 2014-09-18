nucleusParametersTuner
======================

This small program finds the best configuration set of parameters for the
woods-saxon nucleus density distribution for various nuclei 

We assume that the density distribution of individual nucleon is a gaussian

\rho(r) = 1/(2 \pi width)^{1.5} exp(-r^2/(2 width^2))

Since the nucleus charge density is measured at low collision energy, the
fit should also be done at low collision energy. When the collision energy
increases, the size of the individual nucleon grows with the inelastic 
nucleon-nucleon cross section and so does the size of the nucleus. 

We will perform our fit at \sqrt{s} = 23.5 GeV according to the reference, 
[http://arxiv.org/pdf/1108.5379.pdf],
where the gaussian width for the nucleon is 0.473 fm.
