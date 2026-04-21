function RunSim(n, no_sim)
% This function implements Example B.1
%
% Inputs:       
% n            The number of observations
% no_sim       The number of replications

% Define a few variables
bet = 1;          % beta_0
k   = 3*floor(n^(1/3));     % Truncation parameters

sim_b_LS   = zeros(no_sim, 1);
sim_b_1SLS = zeros(no_sim, 1);
sim_b_LSO1 = zeros(no_sim, 1);

sim_Sel = zeros(no_sim, 1);
sim_Sel_All = zeros(no_sim, 20);

for i_sim = 1:no_sim
    % DGP
    x = pi*rand(n, 1);
    e = randn(n, 1);
    y = 1 + x*bet + m_x(x, 1000) + e;
    
    % OLS regression
    x_all = [ones(n, 1), x];
    b_LS = (x_all'*x_all)\(x_all'*y);
    sim_b_LS(i_sim, 1) = b_LS(2);
    
    % 1SLS regression with different k's
    V  = Gen_V(x, k);
    M  = eye(n) - V/(V'*V)*(V');
 
    sim_b_1SLS(i_sim, 1) = (x'*M*x)\(x'*M*y); % Record the estimated b

    % LASSO
    V_LSO = Gen_V(x, k);
    M_x   = eye(n) - x/(x'*x)*(x');
    gam_LSO = lasso(M_x*V_LSO, M_x*y, 'Lambda', (log(n)/n)^(1/2));

    % Check purity
    temp = zeros(length(gam_LSO), 1);
    temp(1, 1) = 1;
    i = 0;
    while i <= k-1
        i = i + 1;
        if rem(i-1, 4) == 3
            temp(i, 1) = 1;
        end
    end
    sim_Sel(i_sim, 1) = mean(abs(temp - (gam_LSO ~= 0)), 1);
    sim_Sel_All(i_sim, 1:length(gam_LSO)) = (gam_LSO ~= 0)';
    
    % Post LASSO
    V_LSO = V_LSO(:, gam_LSO ~= 0);
    temp = 0; % Check whether ones(n,1) is included
    for i = 1:size(V_LSO, 2)
        if V_LSO(:, i) == ones(n, 1)
            temp = temp + 1;
        end
    end
    if temp == 0
        V_LSO = [ones(n, 1), V_LSO]; %#ok<AGROW>
    end
    M = eye(n) - V_LSO/(V_LSO'*V_LSO)*(V_LSO');
    sim_b_LSO1(i_sim, 1) = (x'*M*x)\(x'*M*y);
end

b_LS    = mean(sim_b_LS, 1) - bet;
b_1SLS  = mean(sim_b_1SLS, 1) - bet;
b_LSO1  = mean(sim_b_LSO1, 1) - bet;
LSO_Sel = mean(sim_Sel, 1);

LSO_All = mean(sim_Sel_All, 1)';

sd_LS      = std(sim_b_LS, 0, 1);
sd_1SLS    = std(sim_b_1SLS, 0, 1);
sd_LSO1    = std(sim_b_LSO1, 0, 1);
sd_LSO_Sel = std(sim_Sel, 0, 1);

output = [b_LS, b_1SLS, b_LSO1, LSO_Sel; ...
    sd_LS, sd_1SLS, sd_LSO1, sd_LSO_Sel]';

% Store results
filename = sprintf('RE_n%d.txt', n);
save(filename, 'output', '-ascii')

% Store results
filename = sprintf('Selection_n%d.txt', n);
save(filename, 'LSO_All', '-ascii')

end