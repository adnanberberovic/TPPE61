function [dG] = gradientOIS(n, m, T_k, T_s, f_tilde)


% input example
%T_k    = [0 1/12 3/12 6/12 9/12 1 2 3 4 5 6 7 8 9 10]; m = length(T_k)-1;
%T_s    = [0 1/12 3/12 6/12 9/12 1       5         10]; n = length(T_s)-1;
%f_tilde = 0.01*ones(1, 4*nSplines); 





delta = diff(T_k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Täljare  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = zeros(m,1);

for c = 1:m
    
    T(c) = 1 - exp(-spotRate(T_s, f_tilde, T_k(c+1)));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Nämnare  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = zeros(m,1);

for c = 1:m
    
    for i = 1:c
        
        N(c) = N(c) + delta(i)*exp(-spotRate(T_s, f_tilde, T_k(i+1))*T_k(i+1));
        
    end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  dTäljare  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dT = zeros(4*n, m);

for c = 1:m
    
    maturity = T_k(c+1);
    
    dot = maturity*exp(-spotRate(T_s, f_tilde, maturity)*maturity);
    
    d = 1;
    for s = 1:n
        
        if maturity > T_s(s+1)
            dT(d,c) = dot*(1/4)*(T_s(s+1)-T_s(s))^4;
            dT(d+1,c) = dot*(1/3)*(T_s(s+1)-T_s(s))^3; 
            dT(d+2,c) = dot*(1/2)*(T_s(s+1)-T_s(s))^2;
            dT(d+3,c) = dot*(T_s(s+1)-T_s(s));
            d = d + 4;
        else 
            dT(d,c) = dot*(1/4)*(maturity-T_s(s))^4;
            dT(d+1,c) = dot*(1/3)*(maturity-T_s(s))^3; 
            dT(d+2,c) = dot*(1/2)*(maturity-T_s(s))^2;
            dT(d+3,c) = dot*(maturity-T_s(s));
            
            dT(d+4:end,c) = 0;
            break;
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  dNämnare  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dN = zeros(4*n, m);

for c = 1:m
    
    maturity = T_k(c+1);
    
    d = 1;
    for s = 1:n
        
        if maturity > T_s(s+1)
               
            for i = 1:length(T_k)-1
                if T_k(i) >= T_s(s) && T_k(i) < T_s(s+1)
                    
                    r = spotRate(T_s, f_tilde, T_k(i+1));
                    deltaT = delta(i);
                    
                    dN(d,c)   = dN(d,c)   + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*(1/4)*(T_s(s+1)-T_s(s))^4;
                    dN(d+1,c) = dN(d+1,c) + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*(1/3)*(T_s(s+1)-T_s(s))^3;
                    dN(d+2,c) = dN(d+2,c) + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*(1/2)*(T_s(s+1)-T_s(s))^2;
                    dN(d+3,c) = dN(d+3,c) + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*      (T_s(s+1)-T_s(s));
                
                end
            end
            
            d = d + 4;
        else 
            
            for i = 1:length(T_k)-1
                if T_k(i) >= T_s(s) && T_k(i) < T_s(s+1)
                    
                    r = spotRate(T_s, f_tilde, T_k(i+1));
                    deltaT = delta(i);
                    
                    dN(d,c)   = dN(d,c)   + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*(1/4)*(T_s(s+1)-T_s(s))^4;
                    dN(d+1,c) = dN(d+1,c) + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*(1/3)*(T_s(s+1)-T_s(s))^3;
                    dN(d+2,c) = dN(d+2,c) + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*(1/2)*(T_s(s+1)-T_s(s))^2;
                    dN(d+3,c) = dN(d+3,c) + deltaT*(-T_k(i+1))*exp(-r*T_k(i+1))*      (T_s(s+1)-T_s(s));
                
                end
            end
            
            dN(d+4:end,c) = 0;
            break;
        end

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Gradient  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dG = (dT.*repmat(N', 4*n, 1) - dN.*repmat(T', 4*n, 1)) ./ (dN.^2);
dG(isnan(dG)) = 0;

end



