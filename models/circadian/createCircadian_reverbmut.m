%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2015
%Model is based on:
% H. Matsuno et al., A new regulatory interaction suggested by simulations
% for circadian genetic control mechanism in mammals, 
% Journal of Bioinformatics and Computational Biology 4 1:139-153, 2006


function model = createCircadian_reverbmut()
	switchSpecies = [5,7,12];
	% PER_CRY thresholds
	%switchThresh{1} = [1.4, 1.5, 2.2];
    switchThresh{1} = [1.4, 1.5, 2.2];
	% REV_ERB thresholds
	switchThresh{2} = [1.1];
	% CLOCK_BMAL thresholds
	switchThresh{3} = [1.0];

	levels = cellfun(@length,switchThresh)+1;
	
	
	% Number of state variables
	model.nstates = 12;
	% Initial state distribution
	model.x0(:,1:12) = repmat([0.0;0.001],1,12);
	model.x0(:,5) = [1.999;2.001];
    model.x0(:,8) = [0.99;1.01];
	model.x0(:,12) = [1.99;2.01];
	% Number of modes
	model.nmodes = prod(levels);
	% Initial mode
	model.mode0 = idxToMode(xToIdx(mean(model.x0),switchSpecies,switchThresh),levels);
	
  	% Model parameters
    model.p = [1/5, 1/7, 1/5, 1/10, 1/15, 1/5, 1/10, 1/5, 1/7, 1/5 ...
		1/10, 1/15, 1, 0.05, 1/5, 1/10, 1, 0.05, 1/5, 1 ...
        0.05, 1/10, 1/2, 1/5, 1/10, 1, 0.05, 1/5, 1.5];
    
	for m = 1:model.nmodes
		idx = modeToIdx(m,levels);
		g{m} = getModeGuard(idx,switchSpecies,switchThresh);
	end
	
	for m = 1:model.nmodes
		idx = modeToIdx(m,levels);
		neighbors = getNeighborModes(idx,levels);
		for j=1:length(neighbors)
			model.modes(m).guards(j).target = neighbors(j);
			model.modes(m).guards(j).formula = g{neighbors(j)};
		end
		model.modes(m).ode = @(t,x) circ_ode(t,x,model.p,idx);
	end
	
	
	
	% ODE solver to simulate the model
	model.solver = @ode45;
end

function dx = circ_ode(t,x,p,mode_idx)
	[clockt,dc] = clock(t);

	switch1 = (mode_idx(1) <= 1);
	switch3 = (mode_idx(1) <= 1);
	switch5 = (mode_idx(1) < 1);
	switch6 = (mode_idx(1) >= 3);
	
	switch7 = (mode_idx(2) < 1);
	
	switch2 = (mode_idx(3) >= 1);
	switch4 = (mode_idx(3) >= 1);
	switch8 = (mode_idx(3) >= 1);
	
	%a = [switch1,switch2,switch3,switch4,switch5,switch6,switch7,switch8]
	dx = zeros(12,1);
	% per
    dx(1) = - p(1)*x(1) + p(13)*switch1*switch2 + p(14);
    % PER
    dx(2) = - p(2)*x(2) + p(15)*x(1) - p(16)*x(2)*x(4);
    % cry
    dx(3) = - p(3)*x(3) + p(17)*switch3*switch4 + p(18);
    % CRY
    dx(4) = - p(4)*x(4) + p(19)*x(3) - p(16)*x(2)*x(4);
    % PER/CRY
    dx(5) = - p(5)*x(5) + p(16)*x(2)*x(4);
    % rev-erb
    dx(6) = - p(6)*x(6) + p(20)*switch5*switch8 + p(21);
    % REV-ERB
    dx(7) = - p(7)*x(7);% + p(22)*x(6);
    % clock
    % dx(8) = x(8) - p(8)*x(8) + p(23);
	dx(8) = dc;
    % CLOCK
    dx(9) = - p(9)*x(9) + p(24)*clockt - p(25)*x(9)*x(11);
    % bmal
    dx(10) = - p(10)*x(10) +  p(26)*switch6*switch7 + p(27);
    % BMAL
    dx(11) = - p(11)*x(11) + p(28)*x(10) - p(25)*x(9)*x(11);
    % CLOCK/BMAL
    dx(12) = - p(12)*x(12) + p(25)*x(9)*x(11);
end

function [c,dc] = clock(t)
	clockt = [0,3.35,26.30,54.22,78.79,114.49,131.40,159.32,188.99,231.40,260.13,283.22,316.10,326.30,376.24,397.44,439.06,463.22,505.90,548.32,570.33,584.69,610.20,639.06,661.87,681.34,697.44,720];
	clockx = [1,1.08,1.1616,1.0837,0.8960,0.7023,0.6414,0.7023,0.9565,1.3063,1.4047,1.3438,1.2023,1.1991,1.2600,1.1991,0.9766,0.8988,0.8292,0.7165,0.7454,0.8324,1.0402,1.2715,1.3438,1.2834,1.1066,1.1];
	clocknt = 28;
	scale = 1.0/10.0;
	t = mod(t,clockt(end)*scale);
	for i=2:clocknt
		if (clockt(i)*scale >= t)
			dt = (clockt(i)-clockt(i-1))*scale;
			dx = clockx(i)-clockx(i-1);
			c = clockx(i) + dx*(t-clockt(i)*scale)/dt;
			dc = dx/dt;
			return;
		end
	end
end

function idx = xToIdx(x,switchSpecies,switchThresh)
	for i=1:length(switchSpecies)
		idx(i) = find(x(switchSpecies(i)) <= [switchThresh{i}, Inf],1,'first')-1;
	end
end

function mode = idxToMode(idx,levels)
	w = fliplr(cumprod(fliplr([levels(2:end),1])));
	mode = dot(idx,w)+1;
end

function idx = modeToIdx(mode,levels)
	w = fliplr(cumprod(fliplr([levels(2:end),1])));
	a = mode-1;
	for i=1:length(levels)
		idx(i) = floor(a / w(i));
		a = mod(a,w(i));
	end
end

function g = getModeGuard(idx,switchSpecies,switchThresh)
	gstr = '1 ';
	for i=1:length(idx)
		lims = [-Inf, switchThresh{i}, Inf];
		
		gstr = [gstr, ...
			sprintf('& (x(:,%d)>%f) & (x(:,%d)<%f)',...
					switchSpecies(i),lims(idx(i)+1),...
					switchSpecies(i),lims(idx(i)+2))];
	end
	evalc(sprintf('g = @(x) (%s)',gstr));
end

function neighbors = getNeighborModes(idx,levels)
	neighbors = [];
	for i=1:length(idx)
		idx_tmp = idx;
		if idx(i)-1 >= 0
			idx_tmp(i) = idx(i)-1;
			m = idxToMode(idx_tmp,levels);
			neighbors = [neighbors m];
		end
		if idx(i)+1 < levels(i)
			idx_tmp(i) = idx(i)+1;
			m = idxToMode(idx_tmp,levels);
			neighbors = [neighbors m];
		end
	end
end

