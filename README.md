Data Assimilation toolbox
=========================

Introduction
============

The Data Assimilation toolbox is a set of fortran programs with
the aims to facilitate data assimilation for everyone.
The driver defines the main program and the user provides the code
of his/her own model which we call the solver, because it solves the
model equations.
Here is the principle, the user specifies the action in a namelist file
and driver takes care of everything. Wait a minute, is it that simple?
Of course not, the user can not just specifies the action, cross his/her 
fingers
and have the data assimilation run on an imaginary model and data.
The goal is to be as close as possible to: have the user do
the least possible. The least possible includes:
- providing the solver that can accomplish some predefined tasks.
- providing an interface between the driver and the solver.

The interface between the driver and the solver is made up of two
routines: an initialization routine and a run routine.
When the user provides those prerequisite, the driver takes care
of the well defined data assimilation process. An additional advantage
for a user starting its project after deciding to use this toolbox is that
there is a substantial numbers of data structues and functions that
perform common tasks relevant to numerical models and data assimilation.
The user can then focus on the only think that is relevant
to his/her problem. For example, the toolbox contains supports for:
- multivariate random generator (Gaussian only for now).
- covariance matrix and associated functions.
- explicit localization in covariance matrices.
- observation data structure and manipulation functions.
- gmres solver.
- RBCG solver.
- CSR sparse matrix.
- basic regularization.
- generalized diffusion projection.
- many data structures to keep track of what is going on.
- etc.

Interface between the driver and the solver
===========================================

Initialization routine
----------------------

The initialization routine is call be the driver to initialize
the environmenet for the model:
- allocate dynamically allocated variables
- read model parameters
- etc.

Run routine
-----------
The run routine is the interface routine that performs all the actions
required by the driver. For example, to run a 4DVAR simulation,
the basic action performed by the driver is to request the cost function
and its grandient. The run routine provided by the user must be able to
compute the cost function as well as its gradient at the request of the driver.

Communcation between the driver and the solver
----------------------------------------------

The communications between the driver and the solver are made easy by an
exchange data structure. It contains all the data that must be exchanged between 
the driver and the solver. It is even possible for the solver to be an 
independant program. The only requirement is that it be linked to the 
communication module that defines the exchange data 
structure. The solver can be an external independant program for many reasons:
- the user is a third party user and would like to upgrade to new versions
  as it becomes available.
- the solver is a generated program and can be regenerated at any time.
- the user would like to keep everything separated to make it easy to debug.
- etc.

In any case, the user has the control. He can rewrite the interface as it
please the needs to accomodate the constraints. We had a case where the solver
was a huge independant program. We design the interface program in two parts:
- one part included with the driver
- the other part included in the solver

What is possible
================

For the time being, we have the infrastrucres available for the following:
- 4DVAR with Newton type minimization
- Implicit Particle filter
- Incremental 4DVAR
- Ensemble Kalman Filter

4DVAR with Newton type minimization
-----------------------------------
this was the first subdriver included in the toolbox and has been used
in many projects. It is based on [M1QN3], a solver of large-scale unconstrained 
minimization problems.

This subdriver is fully automated, all what the user need to do is to
implement an interface with a run routine that can compute the
cost function and its gradient when requested.

Implicit Particle filter
------------------------
The implicit particle filter is an additional step on top of the 4DVAR
For now, this subdriver uses the subdriver fdvar_newtom. It has been tested
in two projects. It is also fully automated, all what the user need to do is to
implement an interface with a run routine that can compute the
cost function and its gradient when requested.

Incremental 4DVAR
-----------------
This is a newly added subdriver; it has not yet been tested in the toolbox,
althought all the core functions have been tested in other projects.
The minimization algorithm impleted for this is th RBCG.

Ensemble Kalman Filter
----------------------
The EnKF subdriver is not automated for now. The driver supplies the
library for the EnKF analysis and it is the responsability of the user
to build the ensemble matrix, the ensemble observation and to to call
the analysis function from the library.


The user provides the code for his/her own model, and the tool
The entry point of the toolbox is 
the driver that runs everything. The driver
The user specifies in a namelist

Code organization
=================
Source code
-----------
The code is organized into subdirectories, each subdirectory
contains closely related modules.
For the moment being, the following subdirectories are provided:
- DA contains everything related to the data assimilation.
- tool contains general tools used by the program.
- obs contains everything related to observation and observations processing
- GD is an additional set of tools provided for generalized diffusion
  regularization and projection.
- chkp is a begining of checkpointing module.

External libraries
------------------
This code depends on lapack library for matrix computation.

Compilation
-----------
The compilation is made easy by using the autotools suites


Release notes
=============
The code is used only internally for development and will be released later
when it is ready for distribution. If you are interested in trying it, let me 
know, I will get you started and ready to use the Data Assimilation toolbox
