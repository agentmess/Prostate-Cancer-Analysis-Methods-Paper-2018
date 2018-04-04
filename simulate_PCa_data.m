% simulated data curves

% std_noise seems off, SNR for in vivo is overestimated
exp.R1P = 1/30;  exp.R1L =1/25;  exp.kPL = 0.02; exp.std_noise = 0.004; exp.Tarrival = 4; exp.Tbolus = 12;

N = 21;
load clinical_flips_20160928_42s.mat
flips_all = flips;
[Sscale, Mzscale] = flips_scaling_factors(flips_all, N);
flips = acos(Mzscale);
acq.flips = flips;
TR = 2;


t = [0:N-1]*TR;
Tbolus = exp.Tbolus;
Tarrival = exp.Tarrival;
input_function = gampdf(t+TR - Tarrival,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1

%%

[Mxy Mz] = simulate_2site_model([0 0], [exp.R1P exp.R1L], [exp.kPL 0], acq.flips, TR, input_function);
Sn = Mxy/exp.std_noise + randn(size(Mxy)); % SNR units

figure
plot(t, Sn)

%%
clear Mxy Mz
Nsim = 100;
A = [-exp.R1P-exp.kPL 0
    +exp.kPL -exp.R1L];

Mxy(1:2,1) = [0,0];
Mz(1:2,1) = [0,0];
for n = 2:N
    %        Mz_m = Mz(:,n-1);
    
    % more accurate to spread out input over a number of samples to
    % avoid unrealistically large signal jumps
    for ni = 1:Nsim-1
        Mz(:,(n-2)*Nsim + ni+1) = expm(A*TR/Nsim) * (Mz(:,(n-2)*Nsim + ni) + [input_function(n-1)/Nsim;0]);
    end
    
    Mxy(:,n) =  Mz(:,(n-1)*Nsim) .* sin(flips(:,n));
    Mz(:,(n-1)*Nsim+1) = Mz(:,(n-1)*Nsim) .* cos(flips(:,n));
end

ti = [0:Nsim*(N-1)]/Nsim * TR + TR;
subplot(211)
plot(ti,Mz, t, input_function)
xlabel('time (s)'), ylabel('M_Z')
legend('Pyruvate', 'Lactate', 'Input, u(t)')
subplot(212)
plot(t,Mxy, 'x-')
xlabel('time (s)'), ylabel('M_{XY}')

%%
figure
plot([1:168]*.25, flips_all*180/pi)
xlabel('time (s)'), ylabel('\theta (degrees)')
legend('Pyruvate', 'Lactate')

% kPL, SNR, Tbolus, Tarrival
%%
figure
for kPL_test = [0.00 .02 .04]
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [kPL_test 0], acq.flips, TR, input_function);
    Sn = Mxy/exp.std_noise + randn(size(Mxy)); % SNR units
    plot(t, Sn)
    hold on
end
hold off
title('kPL')


figure
for std_test = [0.00 .01 .02]
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [exp.kPL 0], acq.flips, TR, input_function);
    Sn = Mxy + randn(size(Mxy))*std_test; % MZ units
    plot(t, Sn)
    hold on
end
hold off
title('SNR')


%%
figure
for Tarrival_test = [-5 0 5]+Tarrival
    input_function_test = gampdf(t - Tarrival_test,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
    input_function_test = input_function_test/sum(input_function_test); % normalize so total input magnetization = 1
    
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [exp.kPL 0], acq.flips, TR, input_function_test);
    Sn = Mxy/exp.std_noise;% + randn(size(Mxy)); % SNR units
    plot(t, Sn)
    hold on
end
hold off
title('Tarrival')

%%
figure
for Tbolus_test = [8 12 16]
    input_function_test = gampdf(t - Tarrival,4,Tbolus_test/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
    input_function_test = input_function_test/sum(input_function_test); % normalize so total input magnetization = 1
    
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [exp.kPL 0], acq.flips, TR, input_function_test);
    Sn = Mxy/exp.std_noise + randn(size(Mxy)); % SNR units
    plot(t, Sn)
    hold on
end
hold off
title('Tbolus')



%% match pc9154

Tarrival = 6; Tbolus = 8;
input_function = gampdf(t+TR - Tarrival,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
std_noise = .004;

for kPL_test = [.015]% .015]
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [kPL_test 0], acq.flips, TR, input_function);
    Sn = Mxy/std_noise + randn(size(Mxy)); % SNR units
    %figure
    plot(t, Sn)
    title(['match pc9154, kPL = ' num2str(kPL_test)])
end
%% match pc9375

Tarrival = 12; Tbolus = 6;
input_function = gampdf(t+TR - Tarrival,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
std_noise = .005;

for kPL_test = [.017]
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [kPL_test 0], acq.flips, TR, input_function);
    Sn = Mxy/std_noise + randn(size(Mxy)); % SNR units
    %figure
    plot(t, Sn)
    title(['match pc9375, kPL = ' num2str(kPL_test)])
end
%% match pc9654

Tarrival = 0; Tbolus = 8;
input_function = gampdf(t+TR - Tarrival,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
std_noise = .008;

for kPL_test = [.015]% .033]  %.005
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [kPL_test 0], acq.flips, TR, input_function);
    Sn = Mxy/std_noise + randn(size(Mxy)); % SNR units
    %figure
    plot(t, Sn)
    title(['match pc9654, kPL = ' num2str(kPL_test)])
end

%% match pc10174

Tarrival = 0; Tbolus = 8;
input_function = gampdf(t - Tarrival,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
std_noise = .006;

for kPL_test = .036%[0.003 .041]
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [kPL_test 0], acq.flips, TR, input_function);
    Sn = Mxy/std_noise + randn(size(Mxy)); % SNR units
    %figure
    plot(t, Sn)
    title(['match pc10174, kPL = ' num2str(kPL_test)])
end


%% VFA TRAMP

% 3D dynamic, initial
N = 18;
load tramp_vfa_n144_flips.mat
flips_all = flips;
[Sscale, Mzscale] = flips_scaling_factors(flips_all, N);
flips = acos(Mzscale);
TR = 2;

acq.flips = flips;
t = [1:N]*TR;


%% match pc8190

Tarrival = 8; Tbolus = 14;
input_function = gampdf(t - Tarrival,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
std_noise = .00275;

for kPL_test = [.017]
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [kPL_test 0], acq.flips, TR, input_function);
    Sn = Mxy/std_noise + randn(size(Mxy)); % SNR units
    figure
    plot(t, Sn)
    title(['match pc8190, kPL = ' num2str(kPL_test)])
end


%% match pc8484

Tarrival = 4; Tbolus = 10;
input_function = gampdf(t - Tarrival,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
std_noise = .0048;

for kPL_test = [.018]
    Mxy = simulate_2site_model([0 0], [exp.R1P exp.R1L], [kPL_test 0], acq.flips, TR, input_function);
    Sn = Mxy/std_noise + randn(size(Mxy)); % SNR units
    figure
    plot(t, Sn)
    title(['match pc8484, kPL = ' num2str(kPL_test)])
end