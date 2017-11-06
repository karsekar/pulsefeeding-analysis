function [y] = thres_decay(b,x)
%THRESHOLD Summary of this function goes here
%   Detailed explanation goes here
    y(x > b(1)) = b(3)*exp(-b(2)*(x(x > b(1))-b(1)));
    y(x <= b(1)) = b(3);
    y(find(y == Inf))=5000; %prevents from exploding
    y = y';
end

