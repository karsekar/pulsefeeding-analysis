function [y] = normdist(b,x)
%Gompertz Summary of this function goes here
%   b(1) is the scaling factor, b(2) is the mu, b(3) is the std dev
    y = b(1) * exp( -(x-b(2)).^2/ (2*b(3)^2));
    
end

