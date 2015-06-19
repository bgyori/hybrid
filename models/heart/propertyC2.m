%$\mathbf{F}^{\leq 500}([q==4]) \wedge \mathbf{F}^{\leq 500}(\mathbf{G}^{\leq 100}([q==1]))$.
function truth = propertyC2(q,x)
	dt = 0.5;
	nq1 = 100/dt;
	q4 = find(q==4,1,'first');
	if ~isempty(q4)
		ctr = 0;
		for i = (q4+1):length(q)
			if (q(i)==1)||(q(i)==8)
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
