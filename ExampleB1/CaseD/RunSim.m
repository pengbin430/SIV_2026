function RunSim(n, no_sim)
% This function implements Example B.1
%
% Inputs:       
% n            The number of observations
% no_sim       The number of replications

% Define a few variables
bet = 1;          % beta_0
bet_w = 1;
k   = 2*floor(n^(1/3));     % Truncation parameters

sim_b_LS   = zeros(no_sim, 2);
sim_b_1SLS = zeros(no_sim, 2);
sim_b_LSO1 = zeros(no_sim, 2);

sim_Sel = zeros(no_sim, 1);

sim_Sel_All = zeros(no_sim, 20);

for i_sim = 1:no_sim
    % DGP
    x = pi*rand(n, 1);
    w = randn(n, 1);
    e = randn(n, 1);
    y = 1 + x*bet + w*bet_w + m_x(x) + e;
    
    % OLS regression
    x_all = [ones(n, 1), x, w];
    b_LS = (x_all'*x_all)\(x_all'*y);
    sim_b_LS(i_sim, :) = b_LS(2:3, 1)';

    V  = Gen_V(x, k);
    M  = eye(n) - V/(V'*V)*(V');
 
    sim_b_1SLS(i_sim, 2) = ( (x'*M*x)\(x'*M*y) )'; % Record the estimated b

    % LASSO
    xw = [x, w];
    V_LSO   = Gen_V(x, k);
    M_xw    = eye(n) - xw/(xw'*xw)*(xw');
    gam_LSO = lasso(M_xw*V_LSO, M_xw*y, 'Lambda', (log(n)/n)^(1/2));

    % Need to make sure length(gam_LSO) >= 6
    temp = zeros(length(gam_LSO), 1);
    temp(1) = 1;
    if length(gam_LSO) >= 4
        temp(4) = 1;
    end
    if length(gam_LSO) >= 6
        temp(6) = 1;
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
    sim_b_LSO1(i_sim, :) = ( (xw'*M*xw)\(xw'*M*y) )';
end

b_LS    = mean(sim_b_LS, 1) - [bet, bet_w];
b_1SLS  = mean(sim_b_1SLS, 1) - [bet, bet_w];
b_LSO1  = mean(sim_b_LSO1, 1) - [bet, bet_w];
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