function [y] = threshold_set(b,x)
%THRESHOLD Summary of this function goes here
%   Detailed explanation goes here
    b(1) = 14.179;
    y(x > b(1)) = (x(x > b(1))-b(1))*b(2)+b(3);
    y(x <= b(1)) = b(3);
    if b(1) < 0 %keep linear
        b(1) = 0;
    end
    y = y';
end

