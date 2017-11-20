function [infectious, t] = SEKI_simulate(transmission_rate, exposure_rate, N, e1, Tmax, nclasses)

% Nclasses can be 1 to 3

switch nclasses
    case 1
        i=1;
        t(1)=0;
        infectious(1)=0;
        
        s=N-e1;
        
        while (t(i)<=Tmax && (N-e1-s)<N)
            
            % Rate of leaving current state
            q_exposure = (transmission_rate * (N-e1-s)* s /(N-1));
            q_e1_to_inf = exposure_rate*e1;
            % Total rate of leaving current state
            q = (q_exposure + q_e1_to_inf);
            % Time next event occurs
%             t(i+1)=t(i)+exprnd(1/q);
            t(i+1)=t(i)-log(rand)/q;
            
            u = rand;
            if (u < (q_exposure/q) && s>0)
                % Exposure event
                s=s-1;
                e1=e1+1;
            else
                % Movement between e1 and infectious
                e1=e1-1;
            end
            infectious(i+1)=N-e1-s;
            i=i+1;
        end
    
        
    case 2
        i=1;
        t(1)=0;
        infectious(1)=0;
        
        e2=0;
        s=N-e1;
        
        while (t(i)<=Tmax && (N-e1-e2-s)<N)
            
            % Rate of leaving current state
            q_exposure = (transmission_rate * (N-e1-e2-s)* s /(N-1));
            q_e1_to_e2 = exposure_rate*e1;
            q_e2_to_inf = exposure_rate*e2;
            % Total rate of leaving current state
            q = (q_exposure + q_e1_to_e2 + q_e2_to_inf);
            % Time next event occurs
            % t(i+1)=t(i)+exprnd(1/q);
            t(i+1)=t(i)-log(rand)/q;
            
            u = rand;
            if (u < (q_exposure/q) && s>0)
                % Exposure event
                s=s-1;
                e1=e1+1;
            elseif (u<((q_exposure+q_e1_to_e2)/q) && e1>0)
                % Movement between e1 and e2
                e1=e1-1;
                e2=e2+1;
            else
                % Movement from e2 to infectious
                e2=e2-1;
            end
            infectious(i+1)=N-e1-e2-s;
            i=i+1;
        end




    case 3

        i=1;
        t(1)=0;
        infectious(1)=0;
        
        e3=0;
        e2=0;
        s=N-e1;
        
        while (t(i)<=Tmax && (N-e1-e2-s)<N)
            
            % Rate of leaving current state
            q_exposure = (transmission_rate * (N-e1-e2-e3-s)* s /(N-1));
            q_e1_to_e2 = exposure_rate*e1;
            q_e2_to_e3 = exposure_rate*e2;
            q_e3_to_inf = exposure_rate*e3;
            % Total rate of leaving current state
            q = (q_exposure + q_e1_to_e2 + q_e2_to_e3 + q_e3_to_inf);
            % Time next event occurs
%             t(i+1)=t(i)+exprnd(1/q);
            t(i+1)=t(i)-log(rand)/q;
            
            u = rand;
            if (u < (q_exposure/q) && s>0)
                % Exposure event
                s=s-1;
                e1=e1+1;
            elseif (u<((q_exposure+q_e1_to_e2)/q) && e1>0)
                % Movement between e1 and e2
                e1=e1-1;
                e2=e2+1;
            elseif (u<((q_exposure + q_e1_to_e2 + q_e2_to_e3)/q) && e2>0)
                % Movement between e2 and e3
                e2=e2-1;
                e3=e3+1;
            else
                % Movement from e3 to infectious
                e3=e3-1;
            end
            infectious(i+1)=N-e1-e2-e3-s;
            i=i+1;
        end
        
        
        
        
end
