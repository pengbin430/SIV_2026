function output = m_x(x, k0)
% This function implements the m(x) function of Example B.1
% 
% Inputs:
% x    The input of m(x), n-by-1
% k0   The number of basis function

output = zeros(size(x, 1), 1);
for i = 1:k0
    j = (i-1)*4+3;
    output = output + ( 0.9^j )*4*cos(j*x);
end

end