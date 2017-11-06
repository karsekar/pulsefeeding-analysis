function [ a ] = stderr( inputarray )
%STDERR - Standard Error
%   calcualte the std error
    a = std(inputarray)/sqrt(length(inputarray));

end

