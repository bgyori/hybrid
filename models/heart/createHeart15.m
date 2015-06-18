%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014
%Model is based on:
%A. Bueno-Orovio et al., Minimal model for human ventricular action
%potentials in tissue, J. Theor. Biol. 253:544-560, 2008

% the stimulus epsilon last for 1 millionsecond, 
% 4 modes with stimulus + 4 modes without stimuls

function model = createHeart15()

	% Number of modes
	model.nmodes = 8;
	% Number of state variables
	model.nstates = 5;
	% Initial mode
	model.mode0 = 1;
	% Initial state distribution
	model.x0(:,1) = [0.0;0.001];     %[low;high]
	model.x0(:,2) = [0.99;1.01];	 %[low;high]
    model.x0(:,3) = [0.99;1.01];     %[low;high]
	model.x0(:,4) = [0;0.001];       %[low;high]
    model.x0(:,5) = [0.0;0.0];       %[low;high]
	
  	% Model parameters
        ENDO_TVP   =     1.4506; 
        ENDO_TV1M  =    75.;  
        ENDO_TV2M  =  10.; 
        ENDO_TWP   =   280.0;   
        ENDO_TW1M  =    6.0;     
        ENDO_TW2M  =    140.;  
        ENDO_TS1   =    2.7342;  
        ENDO_TS2   =   2.;   
        ENDO_TFI   =    0.1;  
        ENDO_TO1   =  470.; 
	%ENDO_TO1   =  0.004; diseased state
        ENDO_TO2   =    6.  ; 
        ENDO_TSO1  =   40.0;
        ENDO_TSO2  =    1.2;
        ENDO_TSI   =    2.9013; 
        ENDO_KSO    =  2.;
        ENDO_USO    =  0.65;    
        ENDO_TWINF  =  0.0273; 
        ENDO_THV    =  0.3;  
        ENDO_KWM    =  200.;  
        ENDO_KS     =  2.0994; 
        ENDO_UWM    =  0.016;  
        ENDO_US     =  0.9087;  
        ENDO_UU     =  1.56;
        ENDO_WINF   = 0.78; 
    
	%%%%%%%%%%%%%%%%%%%%%
	% Resting mode
	% ODE definition, return column vector dx [u; w; v; s; t]
	model.modes(1).ode = @(t,x) ([ (1 - 0.0 - x(1)/ENDO_TO1 - 0.0)  ; ((1.0 -(x(1)/ENDO_TWINF) - x(2))/(ENDO_TW1M + (ENDO_TW2M - ENDO_TW1M) * 0.5 * (1.+tanh(ENDO_KWM*(x(1)-ENDO_UWM))))) ; ((1.0-x(3))/ENDO_TV1M) ; ((((1.+tanh(ENDO_KS*(x(1) - ENDO_US))) * 0.5) - x(4))/ENDO_TS1) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(1).guards(1).formula = @(x) (x(1)>=0.006);
	% Target mode of the guard
	model.modes(1).guards(1).target = 2;
	%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%
	% Gate v closing mode
	% ODE definition, return column vector dx
	model.modes(2).ode = @(t,x) ([ (1 - 0.0 - x(1)/ENDO_TO2 - 0.0) ; ((ENDO_WINF-x(2))/(ENDO_TW1M + (ENDO_TW2M - ENDO_TW1M) * 0.5 * (1.+tanh(ENDO_KWM*(x(1)-ENDO_UWM))))) ; (-x(3)/ENDO_TV2M) ; ((((1.+tanh(ENDO_KS*(x(1)-ENDO_US))) * 0.5) - x(4))/ENDO_TS1) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(2).guards(1).formula = @(x) (x(1)>=0.13);
    model.modes(2).guards(2).formula = @(x) (x(1)<0.006);
	% Target mode of the guard
	model.modes(2).guards(1).target = 3;
    model.modes(2).guards(2).target = 1;    
	%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%%%%%%%%
	% Gate w closing mode
	% ODE definition, return column vector dx
	model.modes(3).ode = @(t,x) ([ (1 - 0.0 - 1./(ENDO_TSO1+((ENDO_TSO2-ENDO_TSO1)*(1.+tanh(ENDO_KSO*(x(1) - ENDO_USO)))) * 0.5) - (-x(2) * x(4)/ENDO_TSI)) ; (-x(2)/ENDO_TWP) ; (-x(3)/ENDO_TV2M) ; ((((1.+tanh(ENDO_KS*(x(1)-ENDO_US))) * 0.5) - x(4))/ENDO_TS2) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(3).guards(1).formula = @(x) (x(1)>=0.3);
    model.modes(3).guards(2).formula = @(x) (x(1)<0.13);
	% Target mode of the guard
	model.modes(3).guards(1).target = 4;
    model.modes(3).guards(2).target = 2;
	%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%
	% AP initiation mode
	% ODE definition, return column vector dx
	model.modes(4).ode = @(t,x) ([ (1 - (-x(3) * (x(1) - ENDO_THV) * (ENDO_UU - x(1))/ENDO_TFI) - 1./(ENDO_TSO1+((ENDO_TSO2 - ENDO_TSO1)*(1.+tanh(ENDO_KSO*(x(1) - ENDO_USO)))) * 0.5) - (-x(2) * x(4)/ENDO_TSI)) ; (-x(2)/ENDO_TWP) ; (-x(3)/ENDO_TVP) ; ((((1.+tanh(ENDO_KS*(x(1) - ENDO_US))) * 0.5) - x(4))/ENDO_TS2); 1  ]);
	% Array of guards
	% Formula of the guard
	model.modes(4).guards(1).formula = @(x) (x(1)<0.3);
	% TRANSIENT STIMULATION
	model.modes(4).guards(2).formula = @(x) (x(5)>1);
	% SUSTAINED STIMULATION
	% model.modes(4).guards(2).formula = @(x) (x(5)>500);
	% Target mode of the guard
	model.modes(4).guards(1).target = 3;
	model.modes(4).guards(2).target = 5;
	%%%%%%%%%%%%%%%%%%%    
	
	%%%%%%%%%%%%%%%%%%%%%
	% AP initiation mode w/o stim
	% ODE definition, return column vector dx
	model.modes(5).ode = @(t,x) ([(0 - (-x(3) * (x(1) - ENDO_THV) * (ENDO_UU - x(1))/ENDO_TFI) - 1./(ENDO_TSO1+((ENDO_TSO2 - ENDO_TSO1)*(1.+tanh(ENDO_KSO*(x(1) - ENDO_USO)))) * 0.5) - (-x(2) * x(4)/ENDO_TSI)) ; (-x(2)/ENDO_TWP) ; (-x(3)/ENDO_TVP) ; ((((1.+tanh(ENDO_KS*(x(1) - ENDO_US))) * 0.5) - x(4))/ENDO_TS2); 1  ]);
	% Array of guards
	% Formula of the guard
	model.modes(5).guards(1).formula = @(x) (x(1)<0.3);
	% Target mode of the guard
	model.modes(5).guards(1).target = 6;
	%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%
	% Gate w closing mode w/o stim
	% ODE definition, return column vector dx
	model.modes(6).ode = @(t,x) ([(0 - 0.0 - 1./(ENDO_TSO1+((ENDO_TSO2-ENDO_TSO1)*(1.+tanh(ENDO_KSO*(x(1) - ENDO_USO)))) * 0.5) - (-x(2) * x(4)/ENDO_TSI)) ; (-x(2)/ENDO_TWP) ; (-x(3)/ENDO_TV2M) ; ((((1.+tanh(ENDO_KS*(x(1)-ENDO_US))) * 0.5) - x(4))/ENDO_TS2) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(6).guards(1).formula = @(x) (x(1)>=0.3);
    model.modes(6).guards(2).formula = @(x) (x(1)<0.13);
	% Target mode of the guard
	model.modes(6).guards(1).target = 5;
    model.modes(6).guards(2).target = 7;
	%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%
	% Gate v closing mode w/o stim
	% ODE definition, return column vector dx
	model.modes(7).ode = @(t,x) ([ (0 - 0.0 - x(1)/ENDO_TO2 - 0.0) ; ((ENDO_WINF-x(2))/(ENDO_TW1M + (ENDO_TW2M - ENDO_TW1M) * 0.5 * (1.+tanh(ENDO_KWM*(x(1)-ENDO_UWM))))) ; (-x(3)/ENDO_TV2M) ; ((((1.+tanh(ENDO_KS*(x(1)-ENDO_US))) * 0.5) - x(4))/ENDO_TS1) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(7).guards(1).formula = @(x) (x(1)>=0.13);
    model.modes(7).guards(2).formula = @(x) (x(1)<0.006);
	% Target mode of the guard
	model.modes(7).guards(1).target = 6;
    model.modes(7).guards(2).target = 8;    
	%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%
	% Resting mode w/o stim
	% ODE definition, return column vector dx [u; w; v; s]
	model.modes(8).ode = @(t,x) ([ (0 - 0.0 - x(1)/ENDO_TO1 - 0.0)  ; ((1.0 -(x(1)/ENDO_TWINF) - x(2))/(ENDO_TW1M + (ENDO_TW2M - ENDO_TW1M) * 0.5 * (1.+tanh(ENDO_KWM*(x(1)-ENDO_UWM))))) ; ((1.0-x(3))/ENDO_TV1M) ; ((((1.+tanh(ENDO_KS*(x(1) - ENDO_US))) * 0.5) - x(4))/ENDO_TS1) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(8).guards(1).formula = @(x) (x(1)>=0.006);
	% Target mode of the guard
	model.modes(8).guards(1).target = 7;
	%%%%%%%%%%%%%%%%%%%
	
	
	
	% ODE solver to simulate the model
	model.solver = @ode15s;
end
