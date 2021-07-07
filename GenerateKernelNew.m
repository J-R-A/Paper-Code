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


% l = 0.1;
% 
% n_const = 1/(stdv*(2*pi)^0.5);
% a_exp =-1/2*((kernRange - l)/stdv).^2;
% kernel_r = n_const*exp(a_exp);
% 
% n_const = 1/(stdv*(2*pi)^0.5);
% a_exp =-1/2*((kernRange)/stdv).^2;
% kernel = n_const*exp(a_exp);
% 
% 
% n_const = -(kernRange)/(stdv^3*(2*pi)^0.5);
% a_exp =-1/2*((kernRange)/stdv).^2;
% kernel_d = n_const.*exp(a_exp);
% 
% kernel_a = kernel - kernel_d*l;
% 
% 
% 
% 
% 
% 
% 
% tms = [0:length(dummy)-1]*sampling;
% 
% 
% 
% 
% plot(tms,conv(dummy,kernel,'same'))
% hold on
% plot(tms,conv(dummy,kernel_r,'same'))
% hold on
% plot(tms,conv(dummy,kernel_a,'same'))
% hold on
% plot(tms,conv(dummy,kernel,'same') - conv(dummy,fliplr(kernel_d),'same')*l)
% 
%     
% 
% 
% plot(kernRange,kernel)
% hold on
% plot(kernRange,kernel_r)
% hold on
% plot(kernRange,kernel_a)
% 
% 
% 
% 
% 
