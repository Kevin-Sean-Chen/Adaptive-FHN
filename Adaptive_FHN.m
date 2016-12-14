%Adaptive_FHN
 
% %clear
% %AFHN_Brownian_Interval
% clearvars -except T
% t = 10000; %in milliseconds
% t = 0:.1:t; step = t(2)-t(1); in = 1/step;

%Stimuli 
n = length(t); I = []; m = max(abs(T)); T = T/m;
b = 20; d1 = 20; d2 = 20; ratio = 1.; %bright and dark durations
pulse = 200; Ab = 1.; Ad = -0;  %p for pulses and A for amplitude
for i = 1:pulse
    I = [I Ad*ones(1,floor((1+T(i))*d1*in)) ratio*Ab*ones(1,floor((1+T(i))*b*in)) ];  %making stimuli
    I = [I Ad*ones(1,floor((1+T(i))*d2*in)) Ab*ones(1,floor((1+T(i))*b*in)) ]; %%two period exp.
    %I = [I Ad*ones(1,floor(d1*in)) ratio*Ab*ones(1,floor(b*in)) ];  %making stimuli
    %I = [I Ad*ones(1,floor(d2*in)) Ab*ones(1,floor(b*in)) ]; %%two period exp.
end

%H = 400;
%I = [I 1*ones(1,H*in)];
%I = [I Ad*ones(1,floor(d1*in)) ratio*Ab*ones(1,floor(b*in)) ];  %making stimuli
%I = [I Ad*ones(1,floor(d2*in)) Ab*ones(1,floor(b*in)) ]; %%two period exp.

k = n-length(I); I = [I -0*ones(1,k)];   %%last state exp.
if length(I)>n
    I(n+1:end) = [];  %confirm same size
end
%A_FHN
tau1 = 10; taua = 200; noise = 0; %%%p = 1.3; a0 = 1.1;
dx = zeros(3,n);
dx(:,1) = [-1,-2/3,1.3];
for i = 2:n
    dx(1,i) = dx(1,i-1) + ((dx(1,i-1)-((dx(1,i-1)^3)/3) - dx(2,i-1) + I(1,i-1)))*step+noise*randn*step^0.5;
    dx(2,i) = dx(2,i-1) + ((1/tau1)*(dx(1,i-1)+dx(3,i-1)))*step;
    %%% dynamics of 'a' %%%
    dx(3,i) = dx(3,i-1) + (((((1-p)*a0+p*((a0^3)/3))-p*dx(2,i-1)) - dx(3,i-1))/taua)*step;
    %dx(3,i) = a0;
    
end

% subplot(2,1,2);plot(t,dx(1,:),t,dx(2,:),t,dx(3,:))
% legend('v','w','a')
% subplot(2,1,1);plot(t,I)
% legend('stimuli')
% figure
% plot(t,dx(3,:))
