function [ P ] = getPermutationMatrix( n )
% Generate the permutation matrix
    
    P = zeros(4*n);
    xn_length = 3+n;
    xb_length = 4*n-xn_length;
    
    P(xb_length+1:xb_length+4,1:4) = eye(4);
    
    for i = 5:4*n
        if mod(i-1,4)==0
           P(xb_length+i,i) = 1;
        else
            P(i-5,i) = 1;
        end
    end
    
end
