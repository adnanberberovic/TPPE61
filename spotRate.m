function [r] = spotRate(t, f, n, T_s)

r = 0;

for s = 1:n
    tau = min(T_s(s+1), t);
    r = r + (1/4)*f((s-1)*4 + 1)*(tau - T_s(s))^4 + ...
            (1/3)*f((s-1)*4 + 2)*(tau - T_s(s))^3 + ...
            (1/2)*f((s-1)*4 + 3)*(tau - T_s(s))^2 + ...
                  f((s-1)*4 + 4)*(tau - T_s(s));

    if tau == t
        break;
    end
end

r = r / t;

end
