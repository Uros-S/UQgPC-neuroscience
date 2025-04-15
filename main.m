% Uros Sutulovic, 04/2025

% Script initialisation
clear; close all; clc;
addpath([pwd,'/util/'])
figID = 0;

% Model setting
model = 1;              % 1 = Hindmarsh-Rose model
                        % 2 = Jansen-Rit model
                        % 3 = Epileptor model

recompute = 0;          % 0 = use precomputed data of the paper
                        % 1 = recompute everything from scratch

% Monte Carlo settings
mc_samples = 5000;

switch model 
    %% Hindmarsh-Rose model
    case 1
        % Initial conditions setting
        x_0 = [0,0,0];
        
        % Nominal model parameters settings
        a = 1.0;
        c = 1.0;
        d = 5.0;
        r = 0.01;
        s = 4.0;
        x_R = -8/5;
        
        parameters_vec_A = [a, 3.1, c, d, 3.5, r, s, x_R];
        % uncertainty = {parameter number of interest (as in parameters_vec),
        % length of uncertainty interval,
        % position of uncertainty w.r.t. value in parameters_vec (1=value is the left extreme of interval)}
        uncertainty_A = {2,3.3-3.1,1};
        parameters_vec_B = [a, 3.0, c, d, 2.4, r, s, x_R];
        uncertainty_B = {2,3.15-3.0,1};
        parameters_vec_C = [a, 2.6, c, d, 2.6, r, s, x_R];
        uncertainty_C = {2,2.8-2.6,1};
        parameters_vec_D = [a, 2.5, c, d, 3.8, r, s, x_R];
        uncertainty_D = {5,4.2-3.8,1};

        % Construction of structs for PoCET toolbox
        [states_A,parameters_A,inputs_A] = HR_struct_create(parameters_vec_A,x_0,uncertainty_A);
        [states_B,parameters_B,inputs_B] = HR_struct_create(parameters_vec_B,x_0,uncertainty_B);
        [states_C,parameters_C,inputs_C] = HR_struct_create(parameters_vec_C,x_0,uncertainty_C);
        [states_D,parameters_D,inputs_D] = HR_struct_create(parameters_vec_D,x_0,uncertainty_D);

        max_expansion_order_gal = 35;   % maximum order expansion for gPC Galerkin

    %% Jansen-Rit model
    case 2
        % Initial conditions setting
        x_0 = [0,0,0,0,0,0];    
        
        % Nominal model parameters settings
        A = 3.25;
        B = 22;
        a = 100;
        b = 50;
        V_0 = 6;
        nu_max = 5;
        r = 0.56;
        p = 120;

        parameters_vec_1 = [68, A, B, a, b, V_0, nu_max, r, p];
        parameters_vec_2 = [95, A, B, a, b, V_0, nu_max, r, p];
        parameters_vec_3 = [128, A, B, a, b, V_0, nu_max, r, p];
        parameters_vec_4 = [135, A, B, a, b, V_0, nu_max, r, p];
        parameters_vec_5 = [155, A, B, a, b, V_0, nu_max, r, p];
        parameters_vec_6 = [173, A, B, a, b, V_0, nu_max, r, p];
        % uncertainty = {parameter number of interest (as in parameters_vec),
        % length of uncertainty interval,
        % position of uncertainty w.r.t. value in parameters_vec (1=value is the left extreme of interval)}
        uncertainty = {1,320-120,1};
        
        % Construction of structs for PoCET toolbox
        [states_1,parameters_1,inputs_1] = JR_struct_create(parameters_vec_1,x_0,uncertainty);
        [states_2,parameters_2,inputs_2] = JR_struct_create(parameters_vec_2,x_0,uncertainty);
        [states_3,parameters_3,inputs_3] = JR_struct_create(parameters_vec_3,x_0,uncertainty);
        [states_4,parameters_4,inputs_4] = JR_struct_create(parameters_vec_4,x_0,uncertainty);
        [states_5,parameters_5,inputs_5] = JR_struct_create(parameters_vec_5,x_0,uncertainty);
        [states_6,parameters_6,inputs_6] = JR_struct_create(parameters_vec_6,x_0,uncertainty);

    %% Epileptor model
    case 3
        % Initial conditions setting
        x_0 = [0,-5,3,0,0,0];

        % Nominal model parameters settings
        r_1 = -1.6;
        r_2 = 1;
        I_1 = 3.1;
        I_2 = 0.42;
        tau_0 = 2500;
        m = 0;
        tau_2 = 10;
        gamma = 0.01;
        
        % uncertainty = {parameter number of interest (as in parameters_vec),
        % length of uncertainty interval,
        % position of uncertainty w.r.t. value in parameters_vec (1=value is the left extreme of interval)}
        parameters_vec_1 = [r_1, r_2, 3.1, 0.2, tau_0, m, tau_2, gamma];
        uncertainty_1 = {3,3.6-3.1,1};
        parameters_vec_2 = [r_1, r_2, 3.1, I_2, tau_0, m, tau_2, gamma];
        uncertainty_2 = {3,3.6-3.1,1};
        parameters_vec_3 = [r_1, r_2, I_1, 0.2, tau_0, m, tau_2, gamma];
        uncertainty_3 = {4,0.35-0.2,1};
        parameters_vec_4 = [r_1, r_2, I_1, I_2, 2500, m, tau_2, gamma];
        uncertainty_4 = {5,2750-2500,1};
        
        % Construction of structs for PoCET toolbox
        [states_1,parameters_1,inputs_1] = Epileptor_struct_create(parameters_vec_1,x_0,uncertainty_1);
        [states_2,parameters_2,inputs_2] = Epileptor_struct_create(parameters_vec_2,x_0,uncertainty_2);
        [states_3,parameters_3,inputs_3] = Epileptor_struct_create(parameters_vec_3,x_0,uncertainty_3);
        [states_4,parameters_4,inputs_4] = Epileptor_struct_create(parameters_vec_4,x_0,uncertainty_4);

    otherwise
        disp('Model selected not valid!');
        return;
end

%% Simulation settings
t_initial = 0;
switch model
    case 1
        t_final = 1200; 
        delta_t = 0.01;
        t_trans = 600;
    case 2
        t_final = 2.5;
        delta_t = 0.0001;
        t_trans = 1.75;
    case 3
        t_final = 4500;
        delta_t = 0.01;
end

integrator = 'ode45';

% Overall simulation options
simulation_opts = {t_initial,t_final,delta_t,integrator,mc_samples,0,0};

%% Computation of the desired performance
%% Load ground truths
switch model
    case 1
        tmp = load('A) b unif [3.1, 3.3], I=3.5.mat');
        ground_truth_A = tmp.mean_potential;
        tmp = load('B) b unif [3.0, 3.15], I=2.4.mat');
        ground_truth_B = tmp.mean_potential;
        tmp = load('C) b unif [2.6, 2.8], I=2.6.mat');
        ground_truth_C = tmp.mean_potential;
        tmp = load('D) b=2.5, I unif [3.8, 4.2].mat');
        ground_truth_D = tmp.mean_potential;
    case 2
        tmp = load('C=68.mat');
        ground_truth_1 = tmp.results;
        tmp = load('C=95.mat');
        ground_truth_2 = tmp.results;
        tmp = load('C=128.mat');
        ground_truth_3 = tmp.results;
        tmp = load('C=135.mat');
        ground_truth_4 = tmp.results;
        tmp = load('C=155.mat');
        ground_truth_5 = tmp.results;
        tmp = load('C=173.mat');
        ground_truth_6 = tmp.results;
    case 3
        tmp = load('I_1 unif [3.1, 3.6], I_2=0.2.mat');
        ground_truth_1 = tmp.results;
        tmp = load('I_1 unif [3.1, 3.6], I_2=0.42.mat');
        ground_truth_2 = tmp.results;
        tmp = load('I_2 unif [0.2, 0.35].mat');
        ground_truth_3 = tmp.results;
        tmp = load('tau_0 unif [2500, 2750].mat');
        ground_truth_4 = tmp.results;
end

%% Reference Monte Carlo with 5000 samples
switch model
    case 1
        switch recompute
            case 0
                execution_time_MC = 1205;
                rmse_mean_MC = 0.0062;
                rmse_var_MC = 0.012;
            case 1
                [MC_results_A,simulation_time_MC_A] = surrogate_model(states_A,parameters_A,inputs_A,2,3,simulation_opts,model);
                [MC_results_B,simulation_time_MC_B] = surrogate_model(states_B,parameters_B,inputs_B,2,3,simulation_opts,model);
                [MC_results_C,simulation_time_MC_C] = surrogate_model(states_C,parameters_C,inputs_C,2,3,simulation_opts,model);
                [MC_results_D,simulation_time_MC_D] = surrogate_model(states_D,parameters_D,inputs_D,2,3,simulation_opts,model);

                % Select best performance
                execution_time_MC = min([simulation_time_MC_A,simulation_time_MC_B, ...
                                         simulation_time_MC_C,simulation_time_MC_D]);

                cut_off_indx = find(ground_truth_A.time > 600,1);   % same for all cases

                rmse_mean_MC_A = rms(MC_results_A.x.moments(1,cut_off_indx:end)-ground_truth_A.x.moments(1,cut_off_indx:end));
                rmse_var_MC_A = rms(MC_results_A.x.moments(2,cut_off_indx:end)-ground_truth_A.x.moments(2,cut_off_indx:end));
                rmse_mean_MC_B = rms(MC_results_B.x.moments(1,cut_off_indx:end)-ground_truth_B.x.moments(1,cut_off_indx:end));
                rmse_var_MC_B = rms(MC_results_B.x.moments(2,cut_off_indx:end)-ground_truth_B.x.moments(2,cut_off_indx:end));
                rmse_mean_MC_C = rms(MC_results_C.x.moments(1,cut_off_indx:end)-ground_truth_C.x.moments(1,cut_off_indx:end));
                rmse_var_MC_C = rms(MC_results_C.x.moments(2,cut_off_indx:end)-ground_truth_C.x.moments(2,cut_off_indx:end));
                rmse_mean_MC_D = rms(MC_results_D.x.moments(1,cut_off_indx:end)-ground_truth_D.x.moments(1,cut_off_indx:end));
                rmse_var_MC_D = rms(MC_results_D.x.moments(2,cut_off_indx:end)-ground_truth_D.x.moments(2,cut_off_indx:end));

                rmse_mean_MC = min([rmse_mean_MC_A,rmse_mean_MC_B,rmse_mean_MC_C,rmse_mean_MC_D]);
                rmse_var_MC = min([rmse_var_MC_A,rmse_var_MC_B,rmse_var_MC_C,rmse_var_MC_D]);

            otherwise
                disp('Recomputation option selected not valid!');
                return;
        end
    case 2
        switch recompute
            case 0
                execution_time_MC = 112;
                rmse_mean_MC = 0.017;
                rmse_var_MC = 0.055;
            case 1
                [MC_results_1,simulation_time_MC_1] = surrogate_model(states_1,parameters_1,inputs_1,2,3,simulation_opts,model);
                [MC_results_2,simulation_time_MC_2] = surrogate_model(states_2,parameters_2,inputs_2,2,3,simulation_opts,model);
                [MC_results_3,simulation_time_MC_3] = surrogate_model(states_3,parameters_3,inputs_3,2,3,simulation_opts,model);
                [MC_results_4,simulation_time_MC_4] = surrogate_model(states_4,parameters_4,inputs_4,2,3,simulation_opts,model);
                [MC_results_5,simulation_time_MC_5] = surrogate_model(states_5,parameters_5,inputs_5,2,3,simulation_opts,model);
                [MC_results_6,simulation_time_MC_6] = surrogate_model(states_6,parameters_6,inputs_6,2,3,simulation_opts,model);
        
                % Select best performance
                execution_time_MC = min([simulation_time_MC_1,simulation_time_MC_2, ...
                                         simulation_time_MC_3,simulation_time_MC_4, ...
                                         simulation_time_MC_5,simulation_time_MC_6]);
        
                cut_off_indx = 1;       % transient not discarded
        
                rmse_mean_MC_1 = rms(MC_results_1.y_2.moments(1,cut_off_indx:end)-ground_truth_1.y_2.moments(1,cut_off_indx:end));
                rmse_var_MC_1 = rms(MC_results_1.y_2.moments(2,cut_off_indx:end)-ground_truth_1.y_2.moments(2,cut_off_indx:end));
                rmse_mean_MC_2 = rms(MC_results_2.y_2.moments(1,cut_off_indx:end)-ground_truth_2.y_2.moments(1,cut_off_indx:end));
                rmse_var_MC_2 = rms(MC_results_2.y_2.moments(2,cut_off_indx:end)-ground_truth_2.y_2.moments(2,cut_off_indx:end));
                rmse_mean_MC_3 = rms(MC_results_3.y_2.moments(1,cut_off_indx:end)-ground_truth_3.y_2.moments(1,cut_off_indx:end));
                rmse_var_MC_3 = rms(MC_results_3.y_2.moments(2,cut_off_indx:end)-ground_truth_3.y_2.moments(2,cut_off_indx:end));
                rmse_mean_MC_4 = rms(MC_results_4.y_2.moments(1,cut_off_indx:end)-ground_truth_4.y_2.moments(1,cut_off_indx:end));
                rmse_var_MC_4 = rms(MC_results_4.y_2.moments(2,cut_off_indx:end)-ground_truth_4.y_2.moments(2,cut_off_indx:end));
                rmse_mean_MC_5 = rms(MC_results_5.y_2.moments(1,cut_off_indx:end)-ground_truth_5.y_2.moments(1,cut_off_indx:end));
                rmse_var_MC_5 = rms(MC_results_5.y_2.moments(2,cut_off_indx:end)-ground_truth_5.y_2.moments(2,cut_off_indx:end));
                rmse_mean_MC_6 = rms(MC_results_6.y_2.moments(1,cut_off_indx:end)-ground_truth_6.y_2.moments(1,cut_off_indx:end));
                rmse_var_MC_6 = rms(MC_results_6.y_2.moments(2,cut_off_indx:end)-ground_truth_6.y_2.moments(2,cut_off_indx:end));
                
                rmse_mean_MC = min([rmse_mean_MC_1,rmse_mean_MC_2,rmse_mean_MC_3, ...
                                    rmse_mean_MC_4,rmse_mean_MC_5,rmse_mean_MC_6]);
                rmse_var_MC = min([rmse_var_MC_1,rmse_var_MC_2,rmse_var_MC_3, ...
                                   rmse_var_MC_4,rmse_var_MC_5,rmse_var_MC_6]);
            otherwise
                disp('Recomputation option selected not valid!');
                return;
        end
    case 3
        switch recompute
            case 0
                execution_time_MC = 5228;
                rmse_mean_MC = 0.0034;
                rmse_var_MC = 0.002;
            case 1
                [MC_results_1,simulation_time_MC_1] = surrogate_model(states_1,parameters_1,inputs_1,2,3,simulation_opts,model);
                [MC_results_2,simulation_time_MC_2] = surrogate_model(states_2,parameters_2,inputs_2,2,3,simulation_opts,model);
                [MC_results_3,simulation_time_MC_3] = surrogate_model(states_3,parameters_3,inputs_3,2,3,simulation_opts,model);
                [MC_results_4,simulation_time_MC_4] = surrogate_model(states_4,parameters_4,inputs_4,2,3,simulation_opts,model);
        
                % Select best performance
                execution_time_MC = min([simulation_time_MC_1,simulation_time_MC_2, ...
                                         simulation_time_MC_3,simulation_time_MC_4]);
        
                cut_off_indx = 1;       % transient not discarded
        
                rmse_mean_MC_1 = rms(MC_results_1.x_1.moments(1,cut_off_indx:end)-ground_truth_1.x_1.moments(1,cut_off_indx:end));
                rmse_var_MC_1 = rms(MC_results_1.x_1.moments(2,cut_off_indx:end)-ground_truth_1.x_1.moments(2,cut_off_indx:end));
                rmse_mean_MC_2 = rms(MC_results_2.x_1.moments(1,cut_off_indx:end)-ground_truth_2.x_1.moments(1,cut_off_indx:end));
                rmse_var_MC_2 = rms(MC_results_2.x_1.moments(2,cut_off_indx:end)-ground_truth_2.x_1.moments(2,cut_off_indx:end));
                rmse_mean_MC_3 = rms(MC_results_3.x_1.moments(1,cut_off_indx:end)-ground_truth_3.x_1.moments(1,cut_off_indx:end));
                rmse_var_MC_3 = rms(MC_results_3.x_1.moments(2,cut_off_indx:end)-ground_truth_3.x_1.moments(2,cut_off_indx:end));
                rmse_mean_MC_4 = rms(MC_results_4.x_1.moments(1,cut_off_indx:end)-ground_truth_4.x_1.moments(1,cut_off_indx:end));
                rmse_var_MC_4 = rms(MC_results_4.x_1.moments(2,cut_off_indx:end)-ground_truth_4.x_1.moments(2,cut_off_indx:end));
                
                rmse_mean_MC = min([rmse_mean_MC_1,rmse_mean_MC_2,rmse_mean_MC_3,rmse_mean_MC_4]);
                rmse_var_MC = min([rmse_var_MC_1,rmse_var_MC_2,rmse_var_MC_3,rmse_var_MC_4]);
            otherwise
                disp('Recomputation option selected not valid!');
                return;
        end
end

MC_data = [execution_time_MC,rmse_mean_MC,rmse_var_MC];

%% gPC Galerkin (only for Hindmarsh-Rose model)
switch model
    case 1
        switch recompute
            case 0
                load('HR A Galerkin.mat');
                load('HR B Galerkin.mat');
                load('HR C Galerkin.mat');
                load('HR D Galerkin.mat');
            case 1
                accuracy_results_A = gPC_accuracy(states_A,parameters_A,inputs_A,ground_truth_A,simulation_opts,1,model);
                accuracy_results_B = gPC_accuracy(states_B,parameters_B,inputs_B,ground_truth_B,simulation_opts,1,model);
                accuracy_results_C = gPC_accuracy(states_C,parameters_C,inputs_C,ground_truth_C,simulation_opts,1,model);
                accuracy_results_D = gPC_accuracy(states_D,parameters_D,inputs_D,ground_truth_D,simulation_opts,1,model);
        end

        galerkin_data = {accuracy_results_A; accuracy_results_B; accuracy_results_C; accuracy_results_D};

        % Plot gPC Galerkin results
        figID = plot_galerkin_results(galerkin_data,MC_data,figID);
end


%% gPC collocation
switch model
    case 1
        switch recompute
            case 0
                load('HR A collocation.mat');
                load('HR B collocation.mat');
                load('HR C collocation.mat');
                load('HR D collocation.mat');
            case 1
                accuracy_results_A = gPC_accuracy(states_A,parameters_A,inputs_A,ground_truth_A,simulation_opts,2,model);
                accuracy_results_B = gPC_accuracy(states_B,parameters_B,inputs_B,ground_truth_B,simulation_opts,2,model);
                accuracy_results_C = gPC_accuracy(states_C,parameters_C,inputs_C,ground_truth_C,simulation_opts,2,model);
                accuracy_results_D = gPC_accuracy(states_D,parameters_D,inputs_D,ground_truth_D,simulation_opts,2,model);
        end
 
        collocation_data = {accuracy_results_A; accuracy_results_B; accuracy_results_C; accuracy_results_D};
    case 2
        switch recompute
            case 0
                load('JR 1 collocation.mat');
                load('JR 2 collocation.mat');
                load('JR 3 collocation.mat');
                load('JR 4 collocation.mat');
                load('JR 5 collocation.mat');
                load('JR 6 collocation.mat');
            case 1
                accuracy_results_1 = gPC_accuracy(states_1,parameters_1,inputs_1,ground_truth_1,simulation_opts,2,model);
                accuracy_results_2 = gPC_accuracy(states_2,parameters_2,inputs_2,ground_truth_2,simulation_opts,2,model);
                accuracy_results_3 = gPC_accuracy(states_3,parameters_3,inputs_3,ground_truth_3,simulation_opts,2,model);
                accuracy_results_4 = gPC_accuracy(states_4,parameters_4,inputs_4,ground_truth_4,simulation_opts,2,model);
                accuracy_results_5 = gPC_accuracy(states_5,parameters_5,inputs_5,ground_truth_5,simulation_opts,2,model);
                accuracy_results_6 = gPC_accuracy(states_6,parameters_6,inputs_6,ground_truth_6,simulation_opts,2,model);
        end
        
        collocation_data = {accuracy_results_1; accuracy_results_2; accuracy_results_3; ...
                            accuracy_results_4; accuracy_results_5; accuracy_results_6};
    case 3
        switch recompute
            case 0
                load('Epileptor 1 collocation.mat');
                load('Epileptor 2 collocation.mat');
                load('Epileptor 3 collocation.mat');
                load('Epileptor 4 collocation.mat');
            case 1
                accuracy_results_1 = gPC_accuracy(states_1,parameters_1,inputs_1,ground_truth_1,simulation_opts,2,model);
                accuracy_results_2 = gPC_accuracy(states_2,parameters_2,inputs_2,ground_truth_2,simulation_opts,2,model);
                accuracy_results_3 = gPC_accuracy(states_3,parameters_3,inputs_3,ground_truth_3,simulation_opts,2,model);
                accuracy_results_4 = gPC_accuracy(states_4,parameters_4,inputs_4,ground_truth_4,simulation_opts,2,model);
        end
        
        collocation_data = {accuracy_results_1; accuracy_results_2; accuracy_results_3; accuracy_results_4};
end

figID = plot_collocation_results(collocation_data,MC_data,model,figID);

