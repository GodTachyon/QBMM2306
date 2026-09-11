# QBMM2306
This is a modified version of [OpenQBMM2306](https://github.com/OpenQBMM/OpenQBMM) software. 

# Objective
The objective of this dependecy is to develop a univariate solver PBM for two phase gas - liquid system. At its fundamental level this implies changes of the diameter calculation to be based on the Sauter mean diameter, d32 (3rd moment/2nd moment) along with using the gas phase velocity for transporting the moments in physical space. The d32 is also to be used in evaluating the momentum forces on the phases

# Plan of action
The univarite solver is linked directly withing the "application" folder following simimlar code structure of exisiting mass based NDF solver for two phase flow called "polydisperseBubbleFoam".
The univariate solver exists as a phaseModel selected via the "phaseProperties" file and linked to a new solver called "pbeBubbleFoam"

## Current status
Currently, the d32 & gas phase velocity has been linked but initial testing is exhibiting "unrealistic" behaviour for existing test cases (bubble columns).
To check -
1. This might be linked to the base class which is calling the diameter. Needs a different place to store and evaluate it. The "d()" member function exists and could be linked to d32. This should also resolve momentum force calculation as there are member functions which utilises the "d()" variable.
2. Double check the gas phase flux called for moment advection is correct.
3. Ensure only necessary source term models are called to prevent run time mismatch between existing solvers.
4. Quality of life improvements since code is duplicated in some places.

# Acknowledgements
Credit for origional code - [OpenQBMM2306](https://github.com/OpenQBMM/OpenQBMM) solver maintained by Alberto Passalacqua.
