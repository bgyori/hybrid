%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014

addpath models/circadian
clear variables


% Time-step governing the controller (J guard checks within each
% time-step)
property.dt = 0.05;
% Number of time-steps to simulate
property.nt = 2000;
% Number of time points to check guards inside each time step
property.J = 10;
% Parameters for checking the property
property.alpha = 0.01;
property.delta = 0.01;

% %% Case 1
% model = createCircadian_wt();
% property.formula = @quantpropertyC1;
% %% Case 2
% model = createCircadian_crymut();
% property.formula = @quantpropertyC1;
% %% Case 3
% model = createCircadian_reverbmut();
% property.formula = @quantpropertyC1;
% %% Case 4
% model = createCircadian_wt();
% property.formula = @quantpropertyC2;
% %% Case 5
% model = createCircadian_nopercrydep();
% property.formula = @quantpropertyC2;
% %% Case 6
 model = createCircadian_nopercrydep();
 property.formula = @quantpropertyC1;


tf = SMC(model,property);
