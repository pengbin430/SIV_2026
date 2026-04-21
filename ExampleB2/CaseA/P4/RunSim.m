function RunSim(n, no_sim)
% This function implements Example B2
%
% Inputs:       
% n            The number of observations
% no_sim       The number of replications

% Define a few variables
bet = 1;          % beta_0
k   = (3:8)';     % Truncation parameters
L_k = length(k);

sim_b_LS   = zeros(no_sim, 1);
sim_b_1SLS = zeros(no_sim, 3);
sim_k      = zeros(no_sim, 1);
sim_L_nk   = zeros(no_sim, 1);

for i_sim = 1:no_sim
    % DGP
    x = pi*rand(n, 1);
    e = randn(n, 1);
    y = x*bet + 0.75*sqrt(ceil(3*log(log(n))))*n^(-1/2)*m_x(x) + e;      % Endogeneity
    
    % OLS regression
    b_LS = (x'*x)\(x'*y);
    sim_b_LS(i_sim, 1) = b_LS;
    
    % 1SLS regression with different k's
    IC = zeros(L_k, 1);           % Record the GCV value
    b_1SLS_temp = zeros(L_k, 1);  % Record the estimated b
    for j = 1:L_k
        V = Gen_V(x, k(j));
        M = eye(n) - V/(V'*V)*(V');

        b_1SLS_temp(j, 1) = (x'*M*x)\(x'*M*y);
        
        gam = (V'*V)\(V'*(y - x*b_1SLS_temp(j, 1)));
        Den = (1 - (k(j)+1)/n)^2;
        IC(j, 1) = mean((y - x*b_1SLS_temp(j, 1) - V*gam).^2, 1) / Den;
    end
    [~, ID] = min(IC);
    k0 = k(ID(1)); % Generate k0

    sim_k(i_sim, 1) = k0;
    V = Gen_V(x, k0);
    M = eye(n) - V/(V'*V)*(V');
    sim_b_1SLS(i_sim, :) = (x'*M*x)\(x'*M*y);
    
    % Run endogeneity test
    [No_Endo, ~] = Test_Endo(y, x, k0);
    sim_L_nk(i_sim, 1) = No_Endo;
end

b_LS   = mean(sim_b_LS, 1) - bet;
b_1SLS = mean(sim_b_1SLS, 1) - bet;
k_1S   = mean(sim_k, 1);

sd_LS   = std(sim_b_LS, 0, 1);
sd_1SLS = std(sim_b_1SLS, 0, 1);
sd_k_1S = std(sim_k, 0, 1);

sz = mean(sim_L_nk, 1); 

output = [b_LS, b_1SLS, k_1S; sd_LS, sd_1SLS, sd_k_1S]';

% Store results
filename = sprintf('RE_n%d.txt', n);
save(filename, 'output', '-ascii')

% Store results
filename = sprintf('TestRate_n%d.txt', n);
save(filename, 'sz', '-ascii')

end