## Overview

This repository provides a set of MATLAB codes for simulations used in the following paper:

Goo, J., Sakhanenko, L. & Zhu, D.C. A chi-square type test for time-invariant fiber pathways of the brain. Stat Inference Stoch Process 25, 449–469 (2022). [https://doi.org/10.1007/s11203-022-09268-6](https://doi.org/10.1007/s11203-022-09268-6).

There are two main codes for running simulations.

1. **simulation_example_H0.m** This script simulates an artificial fiber tract whose orientation does not change over time.

2. **simulation_example_H1.m** This script simulates an artificial fiber tract whose orientation changes over time.

The computational time is slower at the following line:

> dnhat_array=dnhat_Xi(DD_array(:,1),n,nx,ny,nz,nt,h,parform); % To get dnhat(Xi) for all i=1,...,n. 

The code is for $\hat{D}_n(U_i),i=1,\dots,n$ in the residual tensor.

This issue has been resolved using a new noise tensor estimator in our second paper. 

The paper has not yet been published, but the corresponding code can be found at the **pLDTI** repository.

More details on the numerical implementation can be found at the aforementioned paper.
 
