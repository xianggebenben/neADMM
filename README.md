# neADMM: Nonconvex generalization of Alternating Direction Method of Multipliers for nonlinear equality constrained problems
This is the implmentation  of the nonlinear-equality Alternating Direction Method of Multipliers (neADMM).

## How to Use

The codes of this paper are in two folders:

1. Numerical examples folder contains two source files example1.m and example2.m.

2. Application folder contains source codes for two applications: 

    (i). For 1-bit compressive sensing, run the file main.m. 
  

    (ii). For vaccine adverse event detection, run the file main.m.

## Cite

@article{WANG2021100009,
title = {Nonconvex generalization of Alternating Direction Method of Multipliers for nonlinear equality constrained problems},
journal = {Results in Control and Optimization},
pages = {100009},
year = {2021},
issn = {2666-7207},
doi = {https://doi.org/10.1016/j.rico.2021.100009},
url = {https://www.sciencedirect.com/science/article/pii/S2666720721000035},
author = {Junxiang Wang and Liang Zhao},
keywords = {Nonconvex ADMM, Nonlinear equality constraints, Spherical constraints, Multi-instance learning},

abstract = {The classic Alternating Direction Method of Multipliers (ADMM) is a popular framework to solve linear-equality constrained problems. In this paper, we extend the ADMM naturally to nonlinear equality-constrained problems, called neADMM. The difficulty of neADMM is to solve nonconvex subproblems. We provide globally optimal solutions to them in two important applications. Experiments on synthetic and real-world datasets demonstrate excellent performance and scalability of our proposed neADMM over existing state-of-the-start methods.}
}
