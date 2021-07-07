function [y] = convolution_spikes_variable_kernel(time, spike_times, lags, stdv)

y = zeros(size(time));
if length(spike_times) ~= 0
    latency_spikes = interp1(time,lags,spike_times);
    
    for sp = 1:length(spike_times)
        y = y + gaussian_filter_shifted(time,spike_times(sp),latency_spikes(sp),stdv);
        
    end
end
        

    