function [r] = spotRate(tSpline, coeffVec, tSpot)

r = 0;


for s = 1:length(tSpline)-1
    
    
    if tSpot > tSpline(s+1)
        r = r + (1/4)*coeffVec(s)*(tSpline(s+1) - tSpline(s))^4 + (1/3)*coeffVec(s+1)*(tSpline(s+1) - tSpline(s))^3 + (1/2)*coeffVec(s+2)*(tSpline(s+1) - tSpline(s))^2 + coeffVec(s+3)*(tSpline(s+1) - tSpline(s));
    else 
        r = r + (1/4)*coeffVec(s)*(tSpot - tSpline(s))^4 + (1/3)*coeffVec(s+1)*(tSpot - tSpline(s))^3 + (1/2)*coeffVec(s+2)*(tSpot - tSpline(s))^2 + coeffVec(s+3)*(tSpot - tSpline(s));
        break;
    end
    
    
    
end


end



% För att istället integrera över varje segment i tMat:
% for i = 1:find(tMat == tSpot);
%     
%     if any(tSpline == tMat(i)) == 1 % if tMat(i) exists in tSpline use that one
%         s = find(tSpline == min(tSpline(tSpline == tMat(i)))); %if this exists in tSpline set equality 
%     else % else use the smallest one of the above
%         s = find(tSpline == min(tSpline(tSpline > tMat(i)))); %if it does not exists in tSpline set inequality 
%     end
%     
%     r = r + (1/4)*coeffVec(s)*(tMat(i+1) - tMat(i))^4 + (1/3)*coeffVec(s+1)*(tMat(i+1) - tMat(i))^3 + (1/2)*coeffVec(s+2)*(tMat(i+1) - tMat(i))^2 + coeffVec(s+3)*(tMat(i+1) - tMat(i));
%     
% end