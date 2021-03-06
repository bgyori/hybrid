%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014
%Model is based on:
%A. Bueno-Orovio et al., Minimal model for human ventricular action
%potentials in tissue, J. Theor. Biol. 253:544-560, 2008

% the stimulus epsilon last for 1 millionsecond, 
% 4 modes with stimulus + 4 modes without stimuls

function model = createHeart_epi_healthy_sust()

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
    EPI_TVP   =     1.4506; 
    EPI_TV1M  =    60.;  
    EPI_TV2M  =  1150.; 
    EPI_TWP   =   200.0;   
    EPI_TW1M  =    60.0;     
    EPI_TW2M  =    15.;          
    EPI_TS1   =    2.7342;  
    EPI_TS2   =   16.;
    EPI_TFI   =    0.11;
    EPI_TO1   =  400.;
	%EPI_TO1   =  0.004; diseased state
    EPI_TO2   =    6.;
    EPI_TSO1  =   30.0181;
    EPI_TSO2  =    0.9957;
    EPI_TSI   =    1.8875;
    EPI_KSO    =  2.0458;
    EPI_USO    =  0.65;        
    EPI_TWINF  =  0.07;
    EPI_THV    =  0.3; 
    EPI_KWM    =  65.;
    EPI_KS     =  2.0994;                  
    EPI_UWM    =  0.03;
    EPI_US     =  0.9087;
    EPI_UU     =  1.55;
    EPI_WINF   = 0.94;	
    
    
    
	%%%%%%%%%%%%%%%%%%%%%
	% Resting mode
	% ODE definition, return column vector dx [u; w; v; s; t]
	model.modes(1).ode = @(t,x) ([ (1 - 0.0 - x(1)/EPI_TO1 - 0.0)  ; ((1.0 -(x(1)/EPI_TWINF) - x(2))/(EPI_TW1M + (EPI_TW2M - EPI_TW1M) * 0.5 * (1.+tanh(EPI_KWM*(x(1)-EPI_UWM))))) ; ((1.0-x(3))/EPI_TV1M) ; ((((1.+tanh(EPI_KS*(x(1) - EPI_US))) * 0.5) - x(4))/EPI_TS1) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(1).guards(1).formula = @(x) (x(1)>=0.006);
	% Target mode of the guard
	model.modes(1).guards(1).target = 2;
	%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%
	% Gate v closing mode
	% ODE definition, return column vector dx
	model.modes(2).ode = @(t,x) ([ (1 - 0.0 - x(1)/EPI_TO2 - 0.0) ; ((EPI_WINF-x(2))/(EPI_TW1M + (EPI_TW2M - EPI_TW1M) * 0.5 * (1.+tanh(EPI_KWM*(x(1)-EPI_UWM))))) ; (-x(3)/EPI_TV2M) ; ((((1.+tanh(EPI_KS*(x(1)-EPI_US))) * 0.5) - x(4))/EPI_TS1) ; 1 ]);
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
	model.modes(3).ode = @(t,x) ([ (1 - 0.0 - 1./(EPI_TSO1+((EPI_TSO2-EPI_TSO1)*(1.+tanh(EPI_KSO*(x(1) - EPI_USO)))) * 0.5) - (-x(2) * x(4)/EPI_TSI)) ; (-x(2)/EPI_TWP) ; (-x(3)/EPI_TV2M) ; ((((1.+tanh(EPI_KS*(x(1)-EPI_US))) * 0.5) - x(4))/EPI_TS2) ; 1 ]);
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
	model.modes(4).ode = @(t,x) ([ (1 - (-x(3) * (x(1) - EPI_THV) * (EPI_UU - x(1))/EPI_TFI) - 1./(EPI_TSO1+((EPI_TSO2 - EPI_TSO1)*(1.+tanh(EPI_KSO*(x(1) - EPI_USO)))) * 0.5) - (-x(2) * x(4)/EPI_TSI)) ; (-x(2)/EPI_TWP) ; (-x(3)/EPI_TVP) ; ((((1.+tanh(EPI_KS*(x(1) - EPI_US))) * 0.5) - x(4))/EPI_TS2); 1  ]);
	% Array of guards
	% Formula of the guard
	model.modes(4).guards(1).formula = @(x) (x(1)<0.3);
	% TRANSIENT STIMULATION
	%model.modes(4).guards(2).formula = @(x) (x(5)>1);
	% SUSTAINED STIMULATION
	model.modes(4).guards(2).formula = @(x) (x(5)>500);
	% Target mode of the guard
	model.modes(4).guards(1).target = 3;
	model.modes(4).guards(2).target = 5;
	%%%%%%%%%%%%%%%%%%%    
	
	%%%%%%%%%%%%%%%%%%%%%
	% AP initiation mode w/o stim
	% ODE definition, return column vector dx
	model.modes(5).ode = @(t,x) ([(0 - (-x(3) * (x(1) - EPI_THV) * (EPI_UU - x(1))/EPI_TFI) - 1./(EPI_TSO1+((EPI_TSO2 - EPI_TSO1)*(1.+tanh(EPI_KSO*(x(1) - EPI_USO)))) * 0.5) - (-x(2) * x(4)/EPI_TSI)) ; (-x(2)/EPI_TWP) ; (-x(3)/EPI_TVP) ; ((((1.+tanh(EPI_KS*(x(1) - EPI_US))) * 0.5) - x(4))/EPI_TS2); 1  ]);
	% Array of guards
	% Formula of the guard
	model.modes(5).guards(1).formula = @(x) (x(1)<0.3);
	% Target mode of the guard
	model.modes(5).guards(1).target = 6;
	%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%
	% Gate w closing mode w/o stim
	% ODE definition, return column vector dx
	model.modes(6).ode = @(t,x) ([(0 - 0.0 - 1./(EPI_TSO1+((EPI_TSO2-EPI_TSO1)*(1.+tanh(EPI_KSO*(x(1) - EPI_USO)))) * 0.5) - (-x(2) * x(4)/EPI_TSI)) ; (-x(2)/EPI_TWP) ; (-x(3)/EPI_TV2M) ; ((((1.+tanh(EPI_KS*(x(1)-EPI_US))) * 0.5) - x(4))/EPI_TS2) ; 1 ]);
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
	model.modes(7).ode = @(t,x) ([ (0 - 0.0 - x(1)/EPI_TO2 - 0.0) ; ((EPI_WINF-x(2))/(EPI_TW1M + (EPI_TW2M - EPI_TW1M) * 0.5 * (1.+tanh(EPI_KWM*(x(1)-EPI_UWM))))) ; (-x(3)/EPI_TV2M) ; ((((1.+tanh(EPI_KS*(x(1)-EPI_US))) * 0.5) - x(4))/EPI_TS1) ; 1 ]);
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
	model.modes(8).ode = @(t,x) ([ (0 - 0.0 - x(1)/EPI_TO1 - 0.0)  ; ((1.0 -(x(1)/EPI_TWINF) - x(2))/(EPI_TW1M + (EPI_TW2M - EPI_TW1M) * 0.5 * (1.+tanh(EPI_KWM*(x(1)-EPI_UWM))))) ; ((1.0-x(3))/EPI_TV1M) ; ((((1.+tanh(EPI_KS*(x(1) - EPI_US))) * 0.5) - x(4))/EPI_TS1) ; 1 ]);
	% Array of guards
	% Formula of the guard
	model.modes(8).guards(1).formula = @(x) (x(1)>=0.006);
	% Target mode of the guard
	model.modes(8).guards(1).target = 7;
	%%%%%%%%%%%%%%%%%%%
	
	
	
	% ODE solver to simulate the model
	model.solver = @ode15s;
end
