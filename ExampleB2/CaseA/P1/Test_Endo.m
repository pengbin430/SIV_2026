function [No_Endo, L_nk_stat] = Test_Endo(y, x, k0)
% This file implements the test procedure
%
% Inputs:
% y            The response, n-by-1
% x            The regressor, n-by-1
% k0           The number of basis functions involved
%
% Outputs:
% No_Endo      1 means no endogeneity
% L_nk_stat    The value of test statistic

% Produce the variables for test
n = length(y);
V = Gen_V(x, k0);
M = eye(n) - V/(V'*V)*(V');
b_hat   = (x'*M*x)\(x'*M*y);
ehat    = y - x*b_hat;
sigma_e = sqrt(mean(ehat.^2, 1));

k_test = ceil(2*log(log(n))):(3*log(log(n)));
V_k = zeros(n, n, length(k_test));
for ik = 1:length(k_test)
    V = Gen_V(x, k_test(ik));
    V_k(:, :, ik) = V*(V');
end

% Generate test statistic
L_nk = 0;
sigma_n2 = 0;
for i = 1:n
    for j = 1:n
        if j ~= i
            temp = 0;
            for ik = 1:length(k_test)
                L_nk = L_nk + ehat(i) * ehat(j) * V_k(i, j, ik);
                temp = temp + V_k(i, j, ik);
            end
            sigma_n2 = sigma_n2 + temp^2;
        end
    end
end
sigma_n2  = sqrt(sigma_n2 * 2 * (sigma_e^4));
L_nk_stat = L_nk/sigma_n2;

% Implement boostrap procedure to obtain critical values
no_BST = 400;

L_nk_BST = zeros(no_BST, 1);
for i_BST = 1:no_BST
    y_BST = x*b_hat + ehat.*randn(n, 1);
    btemp = (x'*M*x)\(x'*M*y_BST);
    e_BST = y_BST - x*btemp;

    L_temp = 0;
    for i = 1:n
        for j = 1:n
            if j ~= i
                for ik = 1:length(k_test)
                    L_temp = L_temp + e_BST(i) * e_BST(j) * V_k(i, j, ik);
                end
            end
        end
    end
    L_nk_BST(i_BST, 1) = L_temp/sigma_n2;
end

% 1 stands for no endogeneity, 0 stands for otherwise
if L_nk_stat <= quantile(L_nk_BST, 0.975) && L_nk_stat >= quantile(L_nk_BST, 0.025)
    No_Endo = 1;
else
    No_Endo = 0;
end

end