%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014

% SMC: Perform statistical model checking on a hybrid system model
% with respect to a property
function tf = SMC(model,property,onesided)
	if nargin < 3
		onesided = true;
	end

    tf = true;
	tStart = tic;
    N = 0;
	nTrue = 0;

    if onesided
        N_max = log(property.alpha)/log(1-property.delta);
		batchSize = 4;
		batchNum = ceil(N_max / batchSize);
		N = 0;
		for j = 1:batchNum
			tf = false(1,batchSize);
			parfor i = 1:batchSize
				% Simulate system
				[q,x] = simTraj(model,property.dt,property.nt,property.J);
				% Verify trajectory
				tf(i) = property.formula(q,x);
				N = N + 1;
			end
			nTrue = nTrue + sum(tf);
			if any(~tf)
				tf = false; 
				break;
			end
			fprintf('N: %d, t:%.1f\n',N,toc(tStart));
		end
    else
    	while true
    		% Simulate system
    		[q,x] = simTraj(model,property.dt,property.nt,property.J);
    		% Verify trajectory
    		tf1 = property.formula(q,x);
            N = N + 1;
            nTrue = nTrue + tf1;
            % Check SPRT condition
    		tf = SPRT(N, nTrue, property.alpha, property.beta, property.delta, property.r);
            if tf ~= -1
                break;
            end
    	end
    end
	tTotal = toc(tStart);
	
	fprintf('\nNumber of samples/true: %d/%d, Runtime: %f\n',N,nTrue,tTotal);
	if tf
		fprintf('H0 accepted\n');
	else
		fprintf('H1 accepted\n');
	end
end

function result = SPRT(nTotal, nTrue, alpha, beta, delta, r)
% SPRT: Sequential probability ratio test
% NTOTAL: Total number of samples evaluated so far
% NTRUE: Number of samples evaluated as true
% ALPHA, BETA, DELTA: SPRT parameters
% R: Probability threshold
% RESULT: 0=Accept H0, 1=Accept H1, -1=No decision

	acceptance_number = floor((log(beta/(1-alpha)) + ...
						nTotal*(log((1-(r+delta))/(1-(r-delta))))) / ...
						((log((r-delta)/(r+delta))- log((1-(r-delta))/(1-(r+delta))))));
    rejection_number = ceil((log((1-beta)/alpha) + ...
						nTotal*(log((1-(r+delta))/(1-(r-delta)))))/ ...
						((log((r-delta)/(r+delta))- log((1-(r-delta))/(1-(r+delta))))));

    result = -1;

    if (nTrue>=acceptance_number)
        result = 1;
    elseif (nTrue<=rejection_number)
        result = 0;
    end

end


% SIMTRAJ: Simulates a trajectory of a hybrid systems model
function [q,x] = simTraj(model,dt,nt,J)
% MODEL: The model to simulate
% DT: time-step governing the controller, one mode change is done within
% each time-step
% NT: number of time-steps to simulate
% J: number of time points to check guards at within each DT

	% Reserve state and mode output
	x = nan(nt+1,model.nstates);
	q = nan(nt+1,1);
	
	% Generate random initial state
	x(1,:) = model.x0(1,:)+rand(1,size(model.x0,2)).*...
						   (model.x0(2,:)-model.x0(1,:));
	% Set initial mode
	q(1) = model.mode0;
	
	% Absolute time
	abst = 0;
	for t=1:nt
		modecur = q(t);
		
		% Generate time points to check guards
		ts = sort(dt*rand(1,J));
		ts = [0 ts dt];
		% Simulate system
		[~,xout] = model.solver(model.modes(modecur).ode,abst+ts,x(t,:));

		% Check which guards are enabled
		g = checkGuards(model,modecur,xout);
		% Find times where none of the guards are enabled
		g1 = ~any(g,2);
		% Add the guard for self-transition when nothing else enabled
		g = [g1 g];
		% Number of time points where each guard is satisfied
		gr = sum(g);

		% Sample index of guard to take
		r = samplemult(gr);
		
		
		if r == 1 % Self-transition
			modenext = modecur;
			x(t+1,:) = xout(end,:);
		else % Not self-transition
			modenext = model.modes(modecur).guards(r-1).target;
			% Find time points where guard is enabled
			tidx = find(g(:,r));
			% Sample time index of transition
			r2 = ceil(rand*length(tidx));
			tchange = tidx(r2);
			if tchange==length(ts) % Change at last time point
				x(t+1,:) = xout(end,:);
			else % Change before last time point
				% Simulate again
				[~,xout] = model.solver(model.modes(modenext).ode,...
						abst+ts(tchange:end)-min(ts(tchange:end)),xout(tchange,:));
				x(t+1,:) = xout(end,:);
			end
		end
		q(t+1) = modenext;
		abst = abst + dt;
	end
end

% CHECKGUARDS: Check which guards in the model's current mode are satisfied
function guards = checkGuards(model,mode,x)
% MODEL: Hybrid ODE model
% MODE: The current mode the model is in
% X: The current state of the model
% GUARDS: A vector of guards that are satisfied in the current mode
	% Number of guards associated with mode
	nguards = length(model.modes(mode).guards);
	
	% Number of states in x
	ns = size(x,1);
	
	% Initialize guard state vector
	guards = zeros(ns,nguards);
	
	for i=1:nguards
		% Evaluate the guard formula
		guards(:,i) = model.modes(mode).guards(i).formula(x);
	end
end

% SAMPLEMULT: multinomial sampling on the distribution d
function didx = samplemult(d)
% D: a vector encoding samples in a multinomial distribution
% DIDX: the index in d of the element chosen randomly
	s = sum(d);
	if s==0
		error('All zero vector');
	end
	dnorm = d/sum(d);
	cd = cumsum(dnorm);
	u = rand;
	didx = binsearch(u,cd);
end

% Binary search used by samplemult
function idx = binsearch(u,cw)
	N = length(cw);
	a = 1;
	b = N;
	while(1)
		i = floor((a+b)/2);
		if(cw(i)<u)
			a = i+1;
		elseif(cw(i)>=u)
			b = i;
		end
		if(a==b), break; end
	end
	idx = a;
end
