%C2: F([x10<=1]/\[x1>=2]/\[x3>=2] /\ F([x10>=1.5]/\[x1<=1]/\[x3<=1] /\ 
% F([x10<=1]/\[x1>=2]/\[x3>=2] /\ F([x10>=1.5]/\[x1<=1]/\[x3<=1]))))
function truth = quantpropertyC2(q,x)
	dt = 0.1;
	truth = false;
	bmal0 = x(:,10);
        per0 = x(:,1);
        cry0 = x(:,3);
	q0 = find(((bmal0<=1)+(per0>=2)+(cry0>=2))>2,1,'first');
	if ~isempty(q0)
		bmal1 = bmal0(q0:length(bmal0));
		per1 = per0(q0:length(per0));
		cry1 = cry0(q0:length(cry0));
		q1 = find(((bmal1>=1.5)+(per1<=1)+(cry1<=1))>2,1,'first');
			if ~isempty(q1)
				bmal2 = bmal1(q1:length(bmal1));
				per2 = per1(q1:length(per1));
				cry2 = cry1(q1:length(cry1));
				q2 = find(((bmal2<=1)+(per2>=2)+(cry2>=2))>2,1,'first');
					if ~isempty(q2)
						bmal3 = bmal2(q2:length(bmal2));
						per3 = per2(q2:length(per2));
						cry3 = cry2(q2:length(cry2));
						q3 = find(((bmal3>=1.5)+(per3<=1)+(cry3<=1))>2,1,'first');
							if ~isempty(q3)
								truth = true;
								return;
							end
					end
			end
	end
end

