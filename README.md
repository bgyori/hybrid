# Approximate probabilistic verification of hybrid systems
## Introduction
This is a MATLAB implementation of an approximate probabilisic verification method for non-linear hybrid systems.
The system dynamics is described by a combination of discrete modes and continuous ordinary differential 
equations governing the temporal evolution of the state variables in each mode. The transitions between modes are 
limited by guards, which are conditioned on the current value of the state variables. 

To verify whether the model satisfies a bounded linear time temporal logic (BLTL) property, we use a probabilistic approximation
technique, and construct a statistical mode checking (SMC) scheme relying on repeated trajectory simulations according to 
Algorithm 1 of [1]. Whether the system satisfies the property is posed as a hypothesis test and decided according to 
Algorithm 2 of [1].

## Model construction
A model is represented as a structure with the following fields. 

* nstates: the number of state variables of the model
* nmodes: the number of modes of the model
* mode0: the index of the initial mode (between 1 and nmodes)
* x0: a matrix of size (nstates,2), each row of which contains the lower and upper bound for the initial condition of a state variable
* modes: a vector of mode structures with nmodes elements
	* modes(i).ode: the function handle corresponding to the ODE equations governing the i-th mode
	* modes(i).guards: a vector of guards whose source is the i-th mode
		* modes(i).guards(j).target: the index of the target mode of the j-th guard
		* modes(i).guards(j).formula: the function handle corresponding to the evaluation of the guard condition. The j-th guard is enabled if the formula evaluates to true. 

## Property construction
A property is represented as a structure with the following fields

* formula: a function handle corresponding to the evaluation a trace with respect to the property
* alpha: the precision of the hypothesis test in the SMC procedure
* delta: the probability threshold parameter for the hypothesis test in the SMC procedure
* dt: the length of the time step for the trajectory simulation
* nt: the number of time steps to simulate
* J: the number of guard evaluations within a single time step

[1] BM Gyori, B Liu, S Paul, R Ramanathan and PS Thiagarajan, Approximate probabilstic verification of hybrid systems <http://arxiv.org/abs/1412.6953>