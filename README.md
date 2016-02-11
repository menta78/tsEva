# tsEva #

This MATLAB toolbox contains an implementation of the Transformed-Stationary (TS) methodology for non-stationary EVA as described in Mentaschi et al. (2016). In synthesis this approach consists in (i) transforming a non-stationary time series into a stationary one to which the stationary EVA theory can be applied; and (ii) reverse-transforming the result into a non-stationary extreme value distribution.

This toolbox is free of external dependencies, and contains calls for both the statistical analysis and the graphical rendering of the results.

A version working on MATLAB R2014b is available at the branch 0.1_R2014b of this project, and can be downloaded [here](https://bitbucket.org/menta78/tseva/get/0.1_R2014b.zip).
The subdirectory "test" in the toolbox contains a sample script which illustrates the features of the toolbox.

## References ##
Mentaschi L, Vousdoukas M, Voukouvalas E, Sartini L, Feyen L, Besio G and Alfieri L. Non-stationary Extreme Value Analysis: a simplified approach for Earth science applications. Submitted to Hydrology and Earth System Sciences (HESS); 2016