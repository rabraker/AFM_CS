% [ trace_ind ] = get_trace_indeces(nperiods, samps_per_period)

% Returns a set of indeces corresponding only to the trace data, so that we
% can discard the retrace.
%
% use like 
% [ trace_ind ] = get_trace_indeces(nperiods, samps_per_period)
% xdat_trace = xdat(trace_ind)
function [ trace_ind ] = get_trace_indeces(nperiods, samps_per_period)


samps_per_line = samps_per_period/2;
trace_ind = [];
for i=0:nperiods-1
    ind = [i*samps_per_period+1:(i*samps_per_period+samps_per_line)]';
    trace_ind = [trace_ind; ind];
end

end

