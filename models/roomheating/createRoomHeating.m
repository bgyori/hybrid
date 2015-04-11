%Approximate probabilistic verification of hybrid systems
%Benjamin M. Gyori, 2014
%Model is based on: 
%Ansgar Fehnker, Franjo Ivancic, Benchmarks for Hybrid Systems Verification
%Hybrid Systems: Computation and Control
%Lecture Notes in Computer Science Volume 2993, 2004, pp 326-341

function model = createRoomHeating(r,h,A,B,C)
	% Check input parameters
	if nargin < 5
		if r ~= 3
			error('Supply, A, B, C parameters.');
		% Set default parameters for r = 3
		else
            A = [-0.9  0.5  0.0
				  0.5 -1.3  0.5
				  0.0  0.5 -0.9];
			B = [0.4;0.3;0.4];
			C = diag([6,7,8]);
		end
	end

	u =		4; % outside temperature
	offT =	21	*ones(1,r);
	onT =	20	*ones(1,r);
	getT =	18	*ones(1,r);
	difT =	1	*ones(1,r);	
	x0L =	20  *ones(1,r);
	x0H =	20.5*ones(1,r);
	

	nroomcomb = nchoosek(r,h);
	nheatercomb = 2^h;
	model.nmodes = nroomcomb*nheatercomb;
	model.nstates = r;
	
	model.mode0 = 1;
	model.x0 = [x0L;x0H];
	
	roomcombs = nchoosek(1:r,h);
	hc = zeros(nroomcomb,r);
	for i=1:nroomcomb
		hc(i,roomcombs(i,:)) = 1;
	end
	
	for i=1:nroomcomb
		for j=1:nheatercomb
			m = (i-1)*nheatercomb + j;
			
			% Set dynamics of mode
			hc1  = hc(i,:);
			heaterson1 = dec2binvec(j-1,h);
			hc1(roomcombs(i,~heaterson1)) = 0;
			model.modes(m).ode = @(t,x) (A*x+B*u+C*hc1');
			
			% Set transitions of mode: heater on/off switch
			kk = 1; % Guard counter
			for k=1:nheatercomb
				mm = (i-1)*nheatercomb + k;
				hc2  = hc(i,:);
				heaterson2 = dec2binvec(k-1,h);
				hc2(roomcombs(i,~heaterson2)) = 0;
				
				hcdiff = hc1 - hc2;
				
				% Exactly 1 on/off difference
				if sum(hcdiff~=0)==1
					room = find(hcdiff~=0);
					% If on in current mode
					if hc1(hcdiff~=0)==1
						model.modes(m).guards(kk).formula = @(x) (x(:,room) >= offT(room));
						model.modes(m).guards(kk).target = mm;
					% Iff off in current mode
					else
						model.modes(m).guards(kk).formula = @(x) (x(:,room) <= onT(room));
						model.modes(m).guards(kk).target = mm;
					end
					kk = kk + 1; % Added a new guard
				end
			end
			
			
			% Set transitions of mode: heater placement
			for k=1:nroomcomb
				% Rooms that are different
				rooms = find(hc(i,:)~=hc(k,:));
				mm = (k-1)*nheatercomb + j;
				% There is difference in 1 heater position
				if length(rooms)==2
					% Heater goes from room 1->2
					if hc(i,rooms(1))== 1
						model.modes(m).guards(kk).formula = @(x) ((x(:,rooms(2)) < getT(rooms(2))) & (x(:,rooms(1))-x(:,rooms(2)) >= difT(rooms(2))));
						model.modes(m).guards(kk).target = mm;
					% Heater goes from room 2->1
					else
						model.modes(m).guards(kk).formula = @(x) ((x(:,rooms(1)) < getT(rooms(1))) & (x(:,rooms(2))-x(:,rooms(1)) >= difT(rooms(1))));
						model.modes(m).guards(kk).target = mm;
					end
					kk = kk + 1; % Added a new guard
				end
			end
			
			% Add this field just for debugging purposes
			% It shows -1 if heater in room is off
			% It shows 1 if heater in room is on
			% It shows 0 if there is no heater in room
			model.modes(m).heaterconf = zeros(1,r);
			for z=1:h
				if heaterson1(z)==1
					model.modes(m).heaterconf(roomcombs(i,z)) = 1;
				else
					model.modes(m).heaterconf(roomcombs(i,z)) = -1;
				end
			end
		end
	end
      
    % ODE solver to simulate the model
	model.solver = @ode45;
end