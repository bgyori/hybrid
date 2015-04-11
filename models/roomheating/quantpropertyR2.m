function truth = quantpropertyR1(m,x)
    truth = true;
	until = 0;
	
    for t=1:40
        if ((x(t,1) >= 18) && (x(t,2) >= 18) && (x(t,3) >= 18) && (x(t,1) <= 22) && (x(t,2) <= 22) && (x(t,3) <= 22))
            until = t;
            break;
        end
    end
	if  (t == 40) && (until == 0)
		truth = false;
        return;
	end
    for t=until:length(x)
        if((x(t,1) < 18) || (x(t,2) < 18) || (x(t,3) < 18) || (x(t,1) > 22) || (x(t,2) > 22) || (x(t,3) > 22))
            truth=false;
            break;
        end
    end
end


