function [data] = ztransform(inputdata,baseline)
%ZTRANS Performs a z transform on the data based on the average/stdev of a
%baseline
%   User provides the time course data and defines the baseline to
%   transform to
    num_of_ions = length(inputdata(1,:));
    num_of_tp = length(inputdata(:,1));
    avg_vector = mean(inputdata(baseline,:));
    std_dev = std(inputdata(baseline,:));
    replicating = repmat(1:num_of_ions,num_of_tp,1);

    data = (inputdata-avg_vector(replicating))  ./std_dev(replicating);

end

