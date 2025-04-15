function output = gPC_accuracy(states,parameters,inputs,ground_truth,simulation_opts,modality,model)
% This function gives computation times and RMSEs w.r.t. the results in 
% ground_truth for increasing order of expansion (until max_expansion_order)
% and, if model == 2 or model == 3, collocation samples for the specified model. 
% If modality == 1 Gaelerkin is used and if modality == 2 collocation is used.
% All the results are then collected in the struct "output".

    if model == 1
        cut_off_indx = find(ground_truth.time > 600,1);
    else
        cut_off_indx = 1;   % transient not discarded
    end
    
    %% Galerkin vs ground truth
    if modality == 1                              
        if model == 2 || model == 3
            disp('Galerkin approach not available for the selected model');
            return;
        end
        
        max_expansion_order = 35;
        
        % Variables to be saved after completing simulation
        execution_time_galerkin = zeros(1,max_expansion_order);
    
        signal_mean_gt = ground_truth.x.moments(1,cut_off_indx:end);
        signal_var_gt = ground_truth.x.moments(2,cut_off_indx:end);
    
        rmse_mean_x = zeros(1,max_expansion_order);
        rmse_variance_x = zeros(1,max_expansion_order);
        
        
        for i=1:max_expansion_order
            simulation_opts{1,7} = i;
        
            [galerkin_results,simulation_time_G] = surrogate_model(states,parameters,inputs,2,1,simulation_opts,model);
        
            execution_time_galerkin(i) = simulation_time_G;
        
            signal_mean_galerkin = galerkin_results.x.moments(1,cut_off_indx:end);
            signal_var_galerkin = galerkin_results.x.moments(2,cut_off_indx:end);
    
            rmse_mean_x(i) = rms(signal_mean_galerkin-signal_mean_gt);
            rmse_variance_x(i) = rms(signal_var_galerkin-signal_var_gt);
            disp(newline);
            disp(['Expansion order reached: ',num2str(i),'/',num2str(max_expansion_order)]);
            disp(newline);
        end
        
        % Collect all the results
        output.execution_time_galerkin = execution_time_galerkin;
        output.rmse_mean = rmse_mean_x;
        output.rmse_variance = rmse_variance_x;
    
    %% collocation vs ground truth
    else   
        if model == 1
            max_expansion_order = 15;
            max_colloc_samples = 500;
            delta_samples = 50;

            signal_mean_gt = ground_truth.x.moments(1,cut_off_indx:end);
            signal_var_gt = ground_truth.x.moments(2,cut_off_indx:end);

        elseif model == 2
            max_expansion_order = 30;
            max_colloc_samples = 1000;
            delta_samples = 50;
            
            signal_mean_gt = ground_truth.x_1.moments(1,cut_off_indx:end);
            signal_var_gt = ground_truth.x_1.moments(2,cut_off_indx:end);

        else
            max_expansion_order = 15;
            max_colloc_samples = 500;
            delta_samples = 50;

            signal_mean_gt = ground_truth.y_2.moments(1,cut_off_indx:end);
            signal_var_gt = ground_truth.y_2.moments(2,cut_off_indx:end);

        end
        
        % Variables to be saved
        execution_time_collocation = zeros(max_expansion_order,int32(max_colloc_samples/delta_samples));
        rmse_mean = zeros(max_expansion_order,int32(max_colloc_samples/delta_samples));
        rmse_variance = zeros(max_expansion_order,int32(max_colloc_samples/delta_samples));
        
        for i=1:max_expansion_order
        
            for j=1:int32(max_colloc_samples/delta_samples)
                simulation_opts{1,7} = i;
                simulation_opts{1,6} = j*delta_samples;
                [colloc_results,simulation_time_C] = surrogate_model(states,parameters,inputs,2,2,simulation_opts,model);
                
                % Expansion order increases from the last row to the first row, while number of samples increases along the column number
                execution_time_collocation(end-i+1,j) = simulation_time_C; 
        
                if model == 1
                    signal_mean_colloc = colloc_results.x.moments(1,cut_off_indx:end);
                    signal_var_colloc = colloc_results.x.moments(2,cut_off_indx:end);     
                elseif model == 2
                    signal_mean_colloc = colloc_results.x_1.moments(1,cut_off_indx:end);
                    signal_var_colloc = colloc_results.x_1.moments(2,cut_off_indx:end);
                elseif model == 3
                    signal_mean_colloc = colloc_results.y_2.moments(1,cut_off_indx:end);
                    signal_var_colloc = colloc_results.y_2.moments(2,cut_off_indx:end);
                end

                rmse_mean(end-i+1,j) = rms(signal_mean_colloc-signal_mean_gt);
                rmse_variance(end-i+1,j) = rms(signal_var_colloc-signal_var_gt);

                disp(newline);
                disp(['Expansion order reached: ',num2str(i),'/',num2str(max_expansion_order)]);
                disp(['Number of collocation samples reached: ',num2str(j*delta_samples),'/',num2str(max_colloc_samples)]);
                disp(newline);
            end
        end

        output.execution_time_collocation = execution_time_collocation;
        output.rmse_mean = rmse_mean;
        output.rmse_variance = rmse_variance;
    end
end