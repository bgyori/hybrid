%C3: F<=500(G<=1([u>=1.4]) && F<=500([0.8<=u<=1.1] && F<=500(G<=100(u>=1.1)))))

function truth = quantpropertyC3(q,x)
	dt = 0.5;
	nq0 = 1/dt;
	nq1 = 100/dt;
	truth = false;

	u = x(:,1);

	%u

        ctr = 0;
	for i = 1:length(q)		
		if u(i) >= 1.4
			ctr = ctr + 1;
		else
			ctr = 0;
		end
		if ctr == nq0
			break;
		end
	end
	
	%fprintf('i=%d\n q.length=%d',i,length(q));

	if i < length(q)
		u2 = u(i:length(u)); 
		q0 = find(((u2>=0.8)+(u2<=1.1))>1,1,'first');
		if ~isempty(q0)
			ctr = 0;
			for j = (q0+1):length(u2)	
				if u2(j)>=1.1
					ctr = ctr + 1;
				else
					ctr = 0;
				end
				%fprintf('j=%d ctr=%d u2=%d\n',j,ctr,u2(j));
				if ctr == 50
					truth = true;
					return;
				end
			end
		end
	end
end
