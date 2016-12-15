function [dG] = gradientOIS(tMat, tSpline, coeffVec)

nContracts = length(tMat)-1;
nSplines = length(tSpline)-1;
delta = diff(tMat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Täljare  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dT = zeros(nContracts, 4*nSplines);

for c = 1:nContracts
    
    maturity = tMat(c+1);
    
    dot = maturity*exp(-spotRate(tSpline, coeffVec, maturity)*maturity);
    
    d = 1;
    for s = 1:nSplines
        
        if maturity > tSpline(s+1)
            dT(c,d) = dot*(1/4)*(tSpline(s+1)-tSpline(s))^4;
            dT(c,d+1) = dot*(1/3)*(tSpline(s+1)-tSpline(s))^3; 
            dT(c,d+2) = dot*(1/2)*(tSpline(s+1)-tSpline(s))^2;
            dT(c,d+3) = dot*(tSpline(s+1)-tSpline(s));
            d = d + 4;
        else 
            dT(c,d) = dot*(1/4)*(maturity-tSpline(s))^4;
            dT(c,d+1) = dot*(1/3)*(maturity-tSpline(s))^3; 
            dT(c,d+2) = dot*(1/2)*(maturity-tSpline(s))^2;
            dT(c,d+3) = dot*(maturity-tSpline(s));
            
            dT(c,d+4:end) = 0;
            break;
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Nämnare  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dN = zeros(nContracts, 4*nSplines);

for c = 1:nContracts
    
    maturity = tMat(c+1);
    
    d = 1;
    for s = 1:nSplines
        
        if maturity > tSpline(s+1)
               
            for i = 1:length(tMat)-1
                if tMat(i) >= tSpline(s) && tMat(i) < tSpline(s+1)
                    
                    r = spotRate(tSpline, coeffVec, tMat(i+1));
                    deltaT = delta(i);
                    
                    dN(c,d)   = dN(c,d)   + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*(1/4)*(tSpline(s+1)-tSpline(s))^4;
                    dN(c,d+1) = dN(c,d+1) + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*(1/3)*(tSpline(s+1)-tSpline(s))^3;
                    dN(c,d+2) = dN(c,d+2) + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*(1/2)*(tSpline(s+1)-tSpline(s))^2;
                    dN(c,d+3) = dN(c,d+3) + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*      (tSpline(s+1)-tSpline(s));
                
                end
            end
            
            d = d + 4;
        else 
            
            for i = 1:length(tMat)-1
                if tMat(i) >= tSpline(s) && tMat(i) < tSpline(s+1)
                    
                    r = spotRate(tSpline, coeffVec, tMat(i+1));
                    deltaT = delta(i);
                    
                    dN(c,d)   = dN(c,d)   + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*(1/4)*(tSpline(s+1)-tSpline(s))^4;
                    dN(c,d+1) = dN(c,d+1) + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*(1/3)*(tSpline(s+1)-tSpline(s))^3;
                    dN(c,d+2) = dN(c,d+2) + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*(1/2)*(tSpline(s+1)-tSpline(s))^2;
                    dN(c,d+3) = dN(c,d+3) + deltaT*(-tMat(i+1))*exp(-r*tMat(i+1))*      (tSpline(s+1)-tSpline(s));
                
                end
            end
            
            dN(c,d+4:end) = 0;
            break;
        end

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Gradient  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dG = dT ./ dN;
dG(isnan(dG)) = 0;



end



