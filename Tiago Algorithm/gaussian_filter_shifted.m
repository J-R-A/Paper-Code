function [kernel] = gaussian_filter_shifted(time,spike,lag,stdv)

time = time - spike;
n_const = 1/(stdv*(2*pi)^0.5);
a_exp =-1/2*((time-lag)/stdv).^2;
kernel = n_const*exp(a_exp);