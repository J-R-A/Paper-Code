function [kernel,diffKernel] = GenerateKernelNew(sampling,stdv)
clearvars -except sampling stdv


kernRange = (0-3*stdv-sampling):sampling:(0+3*stdv+sampling);


%%% gaussian Kernel

n_const = 1/(stdv*(2*pi)^0.5);
a_exp =-1/2*((kernRange)/stdv).^2;
kernel = n_const*exp(a_exp);

%%% gaussian derivative
n_const = -(kernRange)/(stdv^3*(2*pi)^0.5);
diffKernel = n_const.*exp(a_exp);

end
