%C1: F<=500([x10>=1.5] && F<=500([x10<=1] && F<=500([X10>=1.5] && F<=500([X10<=1] && F<=500([x10>=1.5])))))
function truth = quantpropertyC1(q,x)
	dt = 0.1;
	truth = false;
	bmal = x(:,10);
	q0 = find(bmal>=1.5,1,'first');

	if ~isempty(q0)
		bmal2 = bmal(q0:length(bmal));
		q1 = find(bmal2<=1.0,1,'first');
			if ~isempty(q1)
				bmal3 = bmal2(q1:length(bmal2));
				q2 = find(bmal3>=1.5,1,'first');
					if ~isempty(q2)
						bmal4 = bmal3(q2:length(bmal3));
						q3 = find(bmal4<=1.0,1,'first');
							if ~isempty(q3)
								bmal5 = bmal4(q3:length(bmal4));
								q4 = find(bmal5>=1.5,1,'first');
									if ~isempty(q4)
										truth = true;
										return;
									end
							end
					end
			end
	end
end

