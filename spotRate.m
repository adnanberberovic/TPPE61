function [r] = spotRate(t, f, n, T_s)

r = 0;

for s = 1:n
    tau = min(T_s(s+1), t);
    r = r + (1/4)*f((s-1)*4 + 0)*(tau - T_s(s))^4 + ...
            (1/3)*f((s-1)*4 + 1)*(tau - T_s(s))^3 + ...
            (1/2)*f((s-1)*4 + 2)*(tau - T_s(s))^2 + ...
                  f((s-1)*4 + 3)*(tau - T_s(s));

    if tau == t
        break;
    end
end

r = r / t;

end



% F�r att ist�llet integrera �ver varje segment i tMat:
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