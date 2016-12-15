n = 5; % Number of spline
m = 10; % Number of assets
P = zeros(4*n); % Permutation matrix

regul = zeros(4*n);% generateRegularistation(); % Computes regularisation matrix

p = 1; % Penalty for pricing error
E = p*eye(m); % Penalty matrix, zEz 
F = eye(m);
b_e = zeros(m,1); % Actual prices

f = zeros(4*n,1); % Parameters for each spline
f_bar = zeros(size(f)); % Permutated initial guess (f_tilde)
x = zeros(size(f)); % Permutated perturbations (delta_f)
x_B = zeros(3*(n-1),1); % Basic variables
x_N = zeros(n+3,1); % Non-basic variables
f_bar_B = f_bar(1:length(x_B));
f_bar_N = f_bar(length(x_B)+1:end);

% Re-compute every iteration
g = zeros(m,1); % Takes forward interest to ois
grad_g = zeros(4*n,m); % Gradient of function g 
grad_g_bar = P*grad_g; % Gradient with permutated rows
grad_g_bar_B = grad_g_bar(1:length(x_B),:);
grad_g_bar_N = grad_g_bar(length(x_B)+1:end,:);

A_B = -F\grad_g_bar_B';
A_N = -F\grad_g_bar_N';
a = -F\(g-b_e);

B_B = -eye(size(x_B,1));
B_N = zeros(size(x_B,1),size(x_N,1));
b = zeros(size(x_B)); % R.h.s for spline conditions

H = P*regul*P';
H_BB = H(1:length(x_B),1:length(x_B)); 
H_BN = H(1:length(x_B),length(x_B)+1:end);
H_NN = H(length(x_B)+1:end,length(x_B)+1:end);
H_NB = H(length(x_B)+1:end,1:length(x_B));

H_tilde = (B_B\B_N)'*H_BB*(B_B\B_N) - 2*(B_B\B_N)'*H_BN + H_NN + ...
          (A_N - A_B*(B_B\B_N))'*E*(A_N - A_B*(B_B\B_N));
% Assume b = 0
h_tilde = ((-a)'*E*(A_N-A_B*(B_B\B_N)) + f_bar_B'*(H_BN - H_BB*(B_B\B_N)) + ...
          f_bar_N'*(H_NN - H_NB*(B_B\B_N)))';
      
x_N = -H_tilde\h_tilde;
x_B = -(B_B\B_N)*x_N;

f = P'*([f_bar_B; f_bar_N] + [x_B; x_N]);










