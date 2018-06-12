clear all

quick_test = 1;

% default experiment values
exp.R1P = 1/30;  exp.R1L =1/25;  exp.kPL = 0.02; exp.std_noise = 0.004; exp.Tarrival = 4; exp.Tbolus = 8;
%exp.std_noise = 0; %noiseless test
R1P_est = 1/30; R1L_est = 1/25; kPL_est = .02;
R1P_lb = 1/40;  R1P_ub = 1/15;
R1L_lb = 1/35;  R1L_ub = 1/15; % +- 10 s

% good for gammainput:
Rinj_est = 1; Tarrival_est = 4;  A_est = 4; Tbolus_est = 8; B_est = Tbolus_est/4;   % A - shape, B - scale (~ Tbolus/4 for A=4...)
Tarrival_lb = 0; Tarrival_ub = 12;
B_lb = 6 /4; B_ub = 10 /4;

% % for boxcar input
% Tarrival_est = 1.25; Tbolus_est = 16; Rinj_est = 0.06; % boxcar, flips = 2
% Tarrival_est = 2; Tbolus_est = 17; Rinj_est = 0.06; % boxcar, flips = 3

params_all = {'kPL', 'R1L', 'R1P', 'Rinj', 'Tarrival', 'A','B'};

% Options: Fix R1L, Fix Tarrive, Tbolus
% realistically, need to fit kPL and Rinj, but for simulation purposes
% maybe too much
fit_params = {{'kPL','R1L', 'Rinj', 'Tarrival', 'B'}, ...  % 1 - fit everything
    {'kPL', 'Rinj', 'Tarrival','B'}, ...  % 2 - fix R1L, fit bolus
    {'kPL', 'Rinj', 'B'}, ...  % 3- fix R1L, fit bolus duration - hard
    {'kPL', 'Rinj', 'Tarrival'}, ...  % 4- fix R1L, fit bolus arrival - looks easier
    {'kPL','R1L', 'Rinj'}, ... % 5- fix entire bolus, but fit amplitude (Rinj)
    {'kPL', 'Rinj'}, ... % 6 -fix bolus except rate, fix R1L
        {'kPL'}}; % 7- fix R1L and ALL bolus characteristics

fit_desc = {'fitall', ...
    'fixR1L',...
    'fixR1L+Tarrival', ...
    'fixR1L+Tbolus', ...
    'fixbolus', ...
    'fixR1L+bolus', ...
    'fixR1L+bolus+Rinj'};

% % fix Rinj since fitting is just modulating overall SNR
% Fitting Rinj looks to work a little better
% fit_params = {{'kPL','R1L', 'Tarrival', 'Tbolus', 'B'}, ...  % fit everything
%     {'kPL', 'Tarrival', 'Tbolus','B'}, ...  % fix R1L, fit bolus
%     {'kPL', 'Tbolus', 'B'}, ...  % fix R1L, fit bolus duration
%     {'kPL', 'Tarrival'}, ...  % fix R1L, fit bolus arrival
%     {'kPL','R1L'}, ... % fix entire bolus, but fit amplitude (Rinj), maybe add Tarrival for timing?
%     {'kPL'}}; % fix R1L and bolus characteristics


% less apparent variance when fitting Rinj, some more bias though, 
% fixed Rinj more T1 variations

% Models: inputless, boxcar input, gamma variate input [not sure how to put in Tbolus though B = Tbolus/2 A=1 ??  FWHM of gamma function?]
fitting_model = {@fit_kPL, @fit_kPL_withinput, @fit_kPL_withgammainput};
% these correspond to {input-less fitting, fitting with boxcar input, fitting with gamma-variate input shape}

% test multiple fitting options
% Ifitting_mode
for Ifit_params = 2%[2,6] %1:length(fit_params)  % choose which parameter combinations (fit_params) to fit
    for Ifitting_model = [1]  %1:length(fitting_model)  % choose which fitting models (fitting_model) to apply
        
        if Ifit_params>2 && Ifitting_model ==1
            % input-less fitting (fit_kPL) doesn't include bolus parameters
            continue
        end

        % we'll assume best case of a good initial guess, then perturb in further
        % simulations
        
        clear params_est params_fixed acq fitting
        
        for Iparams = 1:length(params_all)
            if find(strcmp(fit_params{Ifit_params}, params_all{Iparams}))
                eval(['params_est.(params_all{Iparams}) =' params_all{Iparams} '_est;'])
                if exist([params_all{Iparams} '_lb'])
                    eval(['params_est.([params_all{Iparams} ''_lb'']) =' params_all{Iparams} '_lb;'])
                end
                if exist([params_all{Iparams} '_ub'])
                    eval(['params_est.([params_all{Iparams} ''_ub'']) =' params_all{Iparams} '_ub;'])
                end
                
            else
                eval(['params_fixed.(params_all{Iparams}) =' params_all{Iparams} '_est;'])
            end
        end
        
        for flip_scheme = [1,4]% [1:4] % choose which flip angle scheme (flip_scheme) to apply
            
            clear flips_all flips
            
            switch flip_scheme
                case 1
                    % 2D dynamic 10/10 flips
                    flip_desc = 'const10';
                    Tacq = 90; TR = 5; N = Tacq/TR;
                    Npe = 8; Nall = N * Npe;
                    flips_all(1:2,1:Nall) = repmat([10*pi/180; 10*pi/180], [1 Nall]);
                    % for same Mz usage
                    flips(1:2,1:N) = repmat(acos(cos([10*pi/180; 10*pi/180]).^Npe), [1 N]);
                    
                case 2
                    % 2D dynamic 10/20 flips
                    flip_desc = 'mband10-20';
                    Tacq = 90; TR = 5; N = Tacq/TR;
                    Npe = 8; Nall = N * Npe;
                    flips_all(1:2,1:Nall) = repmat([10*pi/180; 20*pi/180], [1 Nall]);
                    % for same Mz usage
                    flips(1:2,1:N) = repmat(acos(cos([10*pi/180; 20*pi/180]).^Npe), [1 N]);
                    
                case 3
                    % 3D dynamic, initial
                    flip_desc = 'vfa-tramp';
                    N = 18;
                    load tramp_vfa_n144_flips.mat
                    flips_all = flips;
                    [Sscale, Mzscale] = flips_scaling_factors(flips_all, N);
                    flips = acos(Mzscale);
                    TR = 2;
                    %k12 = 0.05;
                    %flips(1:2,1:N) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ...
                    %    vfa_const_amp(N, pi/2, exp(-TR * ( - k12)))];
                case 4
                    flip_desc = 'vfa-clinical';
                    N = 21;
                    load clinical_flips_20160928_42s.mat
                    flips_all = flips;
                    [Sscale, Mzscale] = flips_scaling_factors(flips_all, N);
                    flips = acos(Mzscale);
                    TR = 2;
                    
            end
                        
            acq.TR = TR; acq.N = N; acq.flips = flips;
            fitting.fit_fcn = fitting_model{Ifitting_model};
            fitting.params_est = params_est; fitting.params_fixed = params_fixed;
            
            if quick_test
                fitting.NMC = 40; %fast testing option
            end
            
            exp_desc = [func2str(fitting_model{Ifitting_model}) '_' fit_desc{Ifit_params} '_' flip_desc];
            
            tic

%            [results(Ifit_params, Ifitting_model, flip_scheme) hdata hsim] = HP_montecarlo_evaluation( acq, fitting, exp );

            [results(Ifit_params, Ifitting_model, flip_scheme) hdata hsim] = HP_montecarlo_evaluation2( acq, fitting, exp );
            toc
            
            set(hsim,'Name', exp_desc);
            if ~quick_test
                figure(hsim)
                print([ exp_desc], '-dpdf')
            end
            
        end
    end
end


return
