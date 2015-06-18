% ===============================================================
% ==   Minimal Resistor Model (4 state variables)              ==
% ==                                                           ==
% ==   Supplement for the CMACS Workshop                       ==
% ==                                                           ==
% ==                                                           ==
% ==   Author:                                                 ==
% ==                                                           ==
% ==     E. Bartocci                                           ==
% ==                                                           == 
% ==   Date:  11/05/10                                         ==
% ==                                                           ==
% ==   Free distribution with authors permission               ==
% ==                                                           ==
% ==   SUNY Stony Brook, Stony Brook, NY                       ==
% ==                                                           == 
% ===============================================================  

% The following are the parameters that you can find in the paper 
% A. Bueno-Orovio, M. Cherry, and F. Fenton, ?Minimal model for 
% human ventricular action potentials in tissue,? Journal of 
% Theoretical Biology, no. 253, pp. 544?560, 2008. 


function [out] = single_cell_mrm4V ()

         global ENDO_TVP; 
         global ENDO_TV1M;
         global ENDO_TV2M; 

         global ENDO_TWP;   

         global ENDO_TW1M; %190    
         global ENDO_TW2M;
         
         global ENDO_TS1;
         global ENDO_TS2;
         global ENDO_TFI;
         global ENDO_TO1;
         global ENDO_TO2;
         global ENDO_TSO1;
         global ENDO_TSO2;

         global ENDO_TSI;


         global ENDO_TWINF;
         global ENDO_THV;

         global ENDO_KWM;
         global ENDO_KS;
         global ENDO_KSO;  
         global ENDO_UWM;    
         global ENDO_US;    
         global ENDO_UU;  
         global ENDO_USO;   
         global ENDO_WINF;

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
         

         ENDO_THVM   =  0.2; % <--
         ENDO_THVINF =  0.006; 
         ENDO_THW    =  0.13; 
         ENDO_THWINF =  0.006;  
         ENDO_THSO   =  0.13; 
         ENDO_THSI   =  0.13; 
         ENDO_THO    =  0.006; 
         ENDO_UO     =  0.;  

         
         ut1     = [];
         vt1     = [];
         wt1     = [];
         st1     = [];
         stimt1  = [];
      
         u = 0.0;
         v = 1.0;
         w = 1.0;
         s = 0.0;
         t1 = 0;
         t2 = 0;

         for i=1:20000   
             [u,v,w,s, t1, t2, stim] = nextStepNonlinear (u,v,w,s, t1, t2, 0.05);
         %for i=1:20000   
             %[u,v,w,s, t1, t2, stim] = nextStepNonlinear (u,v,w,s, t1, t2, 0.0001);
             ut1(i) =    u;
             vt1(i) =    v;
             wt1(i) =    w;
             st1(i) =    s;
             stimt1(i)=  stim;     
         end
         
         
         ph = plot(linspace(0,1000, 20000), [ut1; vt1; wt1; st1; stimt1]);
         %%ph = plot(linspace(0,10, 10000), [ut1; vt1; wt1; st1; stimt1]);
         %ph = plot(linspace(0,2, 20000), [ut1]);
         legend('u', 'v', 'w','s', 'stimulus \epsilon');
         %ylabel('u');
         xlabel('Time (milliseconds)');
         %legend('u state variable', 'v state variable', 'w state variable','s state variable', 'stimulus');
         title ('Minimal Resistor Model');
end


function [u,v,w,s, t1, t2, stim] = nextStepNonlinear (u,v,w,s, t1, t2, dt)
         global ENDO_TVP; 
         global ENDO_TV1M;  
         global ENDO_TV2M; 

         global ENDO_TWP;   

         global ENDO_TW1M; %190    
         global ENDO_TW2M;
         
         global ENDO_TS1;
         global ENDO_TS2;
         global ENDO_TFI;
         global ENDO_TO1;
         global ENDO_TO2;
         global ENDO_TSO1;
         global ENDO_TSO2;

         global ENDO_TSI;


         global ENDO_TWINF;
         global ENDO_THV;

         global ENDO_KWM;
         global ENDO_KS;
         global ENDO_KSO;  
         global ENDO_UWM;    
         global ENDO_US;    
         global ENDO_UU;  
         global ENDO_USO;   
         global ENDO_WINF;


         %% Stimulating the cell at 0, 300, 700 milliseconds
         t1 = t1 + dt;
         t2 = t2 + dt;
         %stim = heaviside(t1 - 0)*(1 - heaviside(t2 - 1)) + heaviside(t1 - 300)*(1 - heaviside(t2 - 301)) + heaviside(t1 - 700)*(1 - heaviside(t2 - 701));
         stim = heaviside(t1 - 0)*(1 - heaviside(t2 - 1)) + heaviside(t1 - 500)*(1 - heaviside(t2 - 501));
         
         
         if u < 0.006
              w = w + ((1.0 -(u/ENDO_TWINF) - w)/(ENDO_TW1M + (ENDO_TW2M - ENDO_TW1M) * 0.5 * (1.+tanh(ENDO_KWM*(u-ENDO_UWM)))))*dt;     
              v = v + ((1.0-v)/ENDO_TV1M)*dt;
              s = s + ((((1.+tanh(ENDO_KS*(u - ENDO_US))) * 0.5) - s)/ENDO_TS1)*dt;
              jfi = 0.0;
              jso = u/ENDO_TO1;
              jsi = 0.0;  
              
         elseif u < 0.13
              w = w + ((ENDO_WINF-w)/(ENDO_TW1M + (ENDO_TW2M - ENDO_TW1M) * 0.5 * (1.+tanh(ENDO_KWM*(u-ENDO_UWM)))))*dt;
              v = v + (-v/ENDO_TV2M)*dt;
              %v = v + ((1.0-v)/ENDO_TV1M)*dt;
              s = s +((((1.+tanh(ENDO_KS*(u-ENDO_US))) * 0.5) - s)/ENDO_TS1)*dt;
              jfi = 0.0;
              jso = u/ENDO_TO2;
              jsi = 0.0;

         elseif u < 0.2
              w = w + (-w/ENDO_TWP)*dt;
              v = v + (-v/ENDO_TV2M)*dt;
              %v = v + ((1.0-v)/ENDO_TV1M)*dt;
              s = s + ((((1.+tanh(ENDO_KS*(u-ENDO_US))) * 0.5) - s)/ENDO_TS2)*dt;
              jfi = 0.0;
              jso = 1./(ENDO_TSO1+((ENDO_TSO2-ENDO_TSO1)*(1.+tanh(ENDO_KSO*(u - ENDO_USO)))) * 0.5);
              jsi = -w * s/ENDO_TSI;            
              
              
              
         elseif u < 0.3
              w = w + (-w/ENDO_TWP)*dt;
              v = v + (-v/ENDO_TV2M)*dt;
              s = s + ((((1.+tanh(ENDO_KS*(u-ENDO_US))) * 0.5) - s)/ENDO_TS2)*dt;
              jfi = 0.0;
              jso = 1./(ENDO_TSO1+((ENDO_TSO2-ENDO_TSO1)*(1.+tanh(ENDO_KSO*(u - ENDO_USO)))) * 0.5);
              jsi = -w * s/ENDO_TSI;
         else
              w  = w + (-w/ENDO_TWP)*dt;
              v  = v + (-v/ENDO_TVP)*dt;
              s  = s +((((1.+tanh(ENDO_KS*(u - ENDO_US))) * 0.5) - s)/ENDO_TS2)*dt;
                       
              jfi = -v * (u - ENDO_THV) * (ENDO_UU - u)/ENDO_TFI;
              jso = 1./(ENDO_TSO1+((ENDO_TSO2 - ENDO_TSO1)*(1.+tanh(ENDO_KSO*(u - ENDO_USO)))) * 0.5);
              jsi = -w * s/ENDO_TSI;  
         end

         u = u  - (jfi+jso+jsi-stim)*dt;
                 
end


