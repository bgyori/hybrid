%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014

addpath models/heart
clear variables


% Create the model
model = createHeart();

% Time-step governing the controller (one guard check within each
% time-step)
property.dt = 0.5;
% Number of time-steps to simulate
property.nt = 1000;
% Number of time points to check guards inside each time step
property.J = 10;
% Parameters for checking the property
property.r = 0.9;
property.alpha = 0.01;
property.beta = property.alpha;
property.delta = 0.01;

% For a property, the formula is a function that receives the mode
% sequence and returns false/true
% C1: The AP mode is reached: F<=500(~[Rest mode])
property.formula = @(q,x) any(q~=1);

% C2: After successfully generating an AP with transmembrane potential above $1.2$, the cardiac cell
% should return to transmembrane potential below $0.006$.
% $\mathbf{F}^{\leq 500}([1.2 \leq u]) \wedge \mathbf{F}^{\leq 500}(\mathbf{G}^{\leq 100}([u \leq
% 0.006]))$.
%property.formula = @quantpropertyC2;

tf = SMC(model,property);
