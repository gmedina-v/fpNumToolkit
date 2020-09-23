# fpNumToolkit
This toolkit written in the Wolfram language allows computation of fixed points and scaling dimensions in the polynomial approximation of the gravitational action with non-interacting matter in the framework of the functional renormalisation group.
The gravitational theory is assummed to be the one where higher curvature interactions consist of a linear combination of the squared Ricci scalar, and Ricci and Riemann tensors evaluated on a spherical background. See http://arxiv.org/abs/1801.00162 and http://arxiv.org/abs/2008.09181 for details.

The main input is a list containing the coefficients of the series expansion of the flow equation around some point.
The couplings must be named \\[Lambda][n][t] where n is an index that can vary from 0 to N-1, with N the order of the polynomial approximation.

Expressions for initialising the flow equations are not included. These can be obtained from existing studies in the literature.
