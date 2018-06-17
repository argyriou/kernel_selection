### Learning the kernel hyperparameters continuously

This algorithm learns convex combinations of continuously parameterized kernels. For example, we can select a kernel among the set of Gaussian kernels whose parameter lies in a given interval.

Like [multiple kernel learning](http://www.jmlr.org/papers/volume5/lanckriet04a/lanckriet04a.pdf), this approach is an alternative to standard hyperparameter tuning by (cross)validation. Thus it requires only the training set and no validation set. In this setting, tuning the kernel parameters is viewed as learning a convex combination of *infinite, continuously parameterized kernels*. The algorithm is detailed in [Learning Convex Combinations of Continuously Parameterized Basic Kernels](http://www0.cs.ucl.ac.uk/staff/M.Pontil/reading/colt05.pdf) and [A DC-Programming Algorithm for Kernel Selection](http://ttic.uchicago.edu/~argyriou/papers/dc-prog.pdf).

The resulting optimization problem is not convex, in general. However, in certain cases of interest (such as Gaussian kernels), a *DC<sup>1</sup> decomposition* of the problem can be readily obtained. The algorithm implemented comes from the area of DC optimization and is tractable for a limited number of kernel parameters.

To run the experiments in [A DC-Programming Algorithm for Kernel Selection](http://ttic.uchicago.edu/~argyriou/papers/dc-prog.pdf), execute the `digit_runs*.m` scripts.

<sup>1</sup>Difference of convex functions.
