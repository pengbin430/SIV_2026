function output = Gen_V(x, k)
% GEN_V generates the V matrix which is n-by-k
% 
% Inputs:
% x       The regressor, n-by-1
% k       The number of basis functions involved
%
% Output: The V matrix which is n-by-k

n = size(x, 1);
output = zeros(n, k);
for j = 1:k
    if j == 1
        output(:, j) = ones(n, 1);
    else
        output(:, j) = sqrt(2)*cos((j-1)*x);
    end
end

end