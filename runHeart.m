%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014

addpath models/heart
clear variables

% Time-step governing the controller (one guard check within each
% time-step)
property.dt = 0.1;
% Number of time-steps to simulate
property.nt = 5000;
% Number of time points to check guards inside each time step
property.J = 10;
% Parameters for checking the property
property.alpha = 0.01;
property.delta = 0.01;

% For a property, the formula is a function that receives the mode
% sequence and returns false/true
% C1: The AP mode is reached: F<=500(~[Rest mode])
%property.formula = @(q,x) any(q~=1);

% C2: After successfully generating an AP, the cardiac cell
% should return to a low transmembrane potential resting mode.
% $\mathbf{F}^{\leq 500}([q==4]) \wedge \mathbf{F}^{\leq 500}(\mathbf{G}^{\leq 100}([q==1]))$.
%property.formula = @propertyC2;

% C3: the epicardial cells has the 'spike-and-dome' AP morphology. Its AP stays above 1.4 for 1 ms 
% and then decrease to [0.8,1.15] (spike). After that it resume increase and stay above 1.15 for 100 ms (dome)  
% F<=500(G<=1([u>=1.4]) && F<=500([0.8<=u<=1.15] && F<=500(G<=100(u>=1.15)))))
% property.formula = @quantpropertyC3;

%% Case 1
model = createHeart_epi_healthy_trans();
property.formula = @(q,x) any(q~=1);
% %% Case 2
% model = createHeart_endo_healthy_trans();
% property.formula = @(q,x) any(q~=1);
% %% Case 3
% model = createHeart_mid_healthy_trans();
% property.formula = @(q,x) any(q~=1);
% %% Case 4
% model = createHeart_epi_diseased_trans();
% property.formula = @(q,x) any(q~=1);
% %% Case 5
% model = createHeart_endo_diseased_trans();
% property.formula = @(q,x) any(q~=1);
% %% Case 6
% model = createHeart_mid_diseased_trans();
% property.formula = @(q,x) any(q~=1);
% %% Case 7
% model = createHeart_epi_healthy_trans();
% property.formula = @propertyC2;
% %% Case 8
% model = createHeart_endo_healthy_trans();
% property.formula = @propertyC2;
% %% Case 9
% model = createHeart_mid_healthy_trans();
% property.formula = @propertyC2;
% %% Case 10
% model = createHeart_epi_healthy_sust();
% property.formula = @propertyC2;
% %% Case 11
% model = createHeart_endo_healthy_sust();
% property.formula = @propertyC2;
% %% Case 12
% model = createHeart_mid_healthy_sust();
% property.formula = @propertyC2;
% %% Case 13
% model = createHeart_epi_healthy_trans();
% property.formula = @quantpropertyC3;
% %% Case 14
% model = createHeart_epi_healthy_trans_ts2();
% property.formula = @quantpropertyC3;
% %% Case 15
% model = createHeart_endo_healthy_sust();
% property.formula = @quantpropertyC3;
% %% Case 16
% model = createHeart_mid_healthy_sust();
% property.formula = @quantpropertyC3;

tf = SMC(model,property);
