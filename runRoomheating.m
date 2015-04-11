%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014

addpath models/roomheating
clear variables

% Time-step governing the controller (one guard check within each
% time-step)
property.dt = 0.025;
% Number of time-steps to simulate
property.nt = 200;
% Number of time points to check guards inside each time step
property.J = 10;
% Parameters for checking the property
property.r = 0.9;
property.alpha = 0.01;
property.beta = property.alpha;
property.delta = 0.01;

% Create the model, 3 rooms, 2 heaters
model = createRoomHeating(3,2);

%% PROPERTIES
% - Property R1 - Initially heaters are in room 1 and 2, and no heater will be moved within 5 days
% $\mathbf{G}^{\leq 5}([\text{Heater in R1}] \wedge \text{[Heater in R2]})$
property.formula = @(q,x) (all((q==1)|(q==2)|(q==3)|q(==4)));


% - Property R2- The temperature in all three rooms stabilizes within 1 day
% between 18 and 22 degrees, and stays in that range for 4 days.
% $\mathbf{F}^{\leq 1}(\mathbf{G}^{\leq 4}([18 < x_1 < 22] \wedge [18 < x_2 < 22] \wedge [18 < x_3 < 22]))$
% property.formula = @quantpropertyR2;


tf = SMC(model,property);
