%C2: F<=500([M4] && F<=500(G<=100([M1])))
function truth = quantpropertyC1(q,x)
	dt = 0.5;
	nq1 = 100/dt;
    u = x(:,1);
	q4 = find(x(:,1)>=1.2,1,'first');
	if ~isempty(q4)
		ctr = 0;
		for i = (q4+1):length(q)
			if u(i)<=0.006
				ctr = ctr + 1;
			else
				ctr = 0;
			end
			if ctr == nq1
				truth = true;
				return;
			end
		end
	end
	truth = false;
end
