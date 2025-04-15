function [results,tEnd] = surrogate_model(states,parameters,inputs,moment_order,surr_model,simulation_opts,model)
% This function computes the first moment_order statistics on the uncertainty
% structure specified by the structs "states" and "parameters" according to
% the surrogate model specified by "surr_model" and the simulation options
% specified by "simulation_opts". The time taken to compute the results is
% returned as well and stored in tEnd.

    % Setting of simulation parameters
    simoptions.tspan = [simulation_opts{1,1}, simulation_opts{1,2}];
    simoptions.dt = simulation_opts{1,3};
    simoptions.setup = odeset;
    simoptions.solver = simulation_opts{1,4};
    
    %% Galerkin approach
    if surr_model == 1

        pce_order = simulation_opts{1,7};

        tStart = tic;

        % Compose the PCE system and write files
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
        PoCETwriteFiles(sys,'PoCET_deterministic_expanded_system',[],'PoCET_nominal_system');

        % Run galerkin-PCE simulation
        results = PoCETsimGalerkin(sys,'PoCET_deterministic_expanded_system',[],simoptions);  
        
        % Compute moments from simulation results
        sys.MomMats = PoCETmomentMatrices(sys,moment_order);
        results.x.moments = PoCETcalcMoments(sys,sys.MomMats,results.x.pcvals);
        
        % Computation time
        tEnd = toc(tStart);
    
    %% Collocation approach
    elseif surr_model == 2
        
        pce_order = simulation_opts{1,7};
        colloc_samples = simulation_opts{1,6};

        tStart = tic;

        % Compose the PCE system and write files
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
        if model ~= 3
            PoCETwriteFiles(sys,'PoCET_deterministic_expanded_system',[],'PoCET_nominal_system');
        end

        % Run collocation-PCE simulation
        basis = PoCETsample(sys,'basis',colloc_samples);    % 'basis' samples from the stochastic basis \xi and not the actual random variable
        if model == 3
            results = PoCETsimCollocation(sys,'PoCET_Epileptor_nominal_system',[],basis,simoptions);  % separate file needed for automatic generation problems with piecewise functions for f_1 and f_2
        else
            results = PoCETsimCollocation(sys,'PoCET_nominal_system',[],basis,simoptions);
        end  
        
        % Compute moments from simulation results for quantities of interest
        sys.MomMats = PoCETmomentMatrices(sys,moment_order);
        if model == 1
            results.x.moments = PoCETcalcMoments(sys,sys.MomMats,results.x.pcvals);
        elseif model == 2
            results.y_2.moments = PoCETcalcMoments(sys,sys.MomMats,results.y_2.pcvals);
        elseif model == 3
            results.x_1.moments = PoCETcalcMoments(sys,sys.MomMats,results.x_1.pcvals);
        end

        % Computation time
        tEnd = toc(tStart);
    
    %% Monte Carlo
    elseif surr_model == 3

        mc_samples = simulation_opts{1,5};
        pce_order = simulation_opts{1,7};

        if pce_order <= 0
            pce_order = 1;  % not important for Monte Carlo
        end

        tStart = tic;

        % Compose the PCE system and write files
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);
        if model ~= 3
            PoCETwriteFiles(sys,'PoCET_deterministic_expanded_system',[],'PoCET_nominal_system');
        end

        % Run Monte-Carlo simulations
        if model == 3
            results = PoCETsimMonteCarlo(sys,'PoCET_Epileptor_nominal_system',[],mc_samples,simoptions,'method','moments');   % separate file needed for automatic generation problems with piecewise functions for f_1 and f_2
        else
            results = PoCETsimMonteCarlo(sys,'PoCET_nominal_system',[],mc_samples,simoptions,'method','moments');
        end

        % Computation time
        tEnd = toc(tStart);
    end

end