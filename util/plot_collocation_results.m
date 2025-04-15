function figID = plot_collocation_results(collocation_data,MC_data,model,figID)
% This function plots the accuracy results for gPC collocation.

    % Plot data initialisation
    execution_time_collocation_max = zeros(size(collocation_data{1,1}.execution_time_collocation,1),size(collocation_data{1,1}.execution_time_collocation,2));
    execution_time_collocation_min = Inf(size(collocation_data{1,1}.execution_time_collocation,1),size(collocation_data{1,1}.execution_time_collocation,2));
    rmse_mean_max = zeros(size(collocation_data{2,1}.rmse_mean,1),size(collocation_data{2,1}.rmse_mean,2));
    rmse_mean_min = Inf(size(collocation_data{2,1}.rmse_mean,1),size(collocation_data{2,1}.rmse_mean,2));
    rmse_variance_max = zeros(size(collocation_data{3,1}.rmse_variance,1),size(collocation_data{3,1}.rmse_variance,2));
    rmse_variance_min = Inf(size(collocation_data{3,1}.rmse_variance,1),size(collocation_data{3,1}.rmse_variance,2));
    
    % Extract collocation data
    for i=1:size(collocation_data,1)
        execution_time_collocation_max = max(execution_time_collocation_max,collocation_data{i,1}.execution_time_collocation);
        execution_time_collocation_min = min(execution_time_collocation_min,collocation_data{i,1}.execution_time_collocation);
        rmse_mean_max = max(rmse_mean_max,collocation_data{i,1}.rmse_mean);
        rmse_mean_min = min(rmse_mean_min,collocation_data{i,1}.rmse_mean);
        rmse_variance_max = max(rmse_variance_max,collocation_data{i,1}.rmse_variance);
        rmse_variance_min = min(rmse_variance_min,collocation_data{i,1}.rmse_variance);
    end

    % Discard data with 50 and 100 samples and flip data for plotting
    execution_time_collocation_max = flip(execution_time_collocation_max);
    execution_time_collocation_max = execution_time_collocation_max(:,3:end);
    execution_time_collocation_min = flip(execution_time_collocation_min);
    execution_time_collocation_min = execution_time_collocation_min(:,3:end);
    rmse_mean_max = rmse_mean_max(:,3:end);
    rmse_mean_max = flip(rmse_mean_max);
    rmse_mean_min = rmse_mean_min(:,3:end);
    rmse_mean_min = flip(rmse_mean_min);
    rmse_variance_max = rmse_variance_max(:,3:end);
    rmse_variance_max = flip(rmse_variance_max);
    rmse_variance_min = rmse_variance_min(:,3:end);
    rmse_variance_min = flip(rmse_variance_min);

    % Extract reference Monte Carlo with 5000 samples
    execution_time_monte_carlo_value = MC_data(1);
    RMSE_mean_monte_carlo_value = MC_data(2);
    RMSE_var_monte_carlo_value = MC_data(3);

    % Plot settings for each model
    mask = ones(size(execution_time_collocation_max,1),size(execution_time_collocation_max,2)); % for excluding the undersampling region
    switch model
        case 1
            delta_z_mean = 0.025;
            delta_z_var = 0.25;
            max_colloc_samples = 500;  
            delta_samples = 50;
            max_expansion_order = 15;
            mask(1:5,1) = NaN(5,1);
            mask(1:4,2) = NaN(4,1);
            mask(1:3,3) = NaN(3,1);
            mask(1:2,4) = NaN(2,1);
            mask(1,5) = NaN;
        case 2
            delta_z_mean = 0.25;
            delta_z_var = 20;
            max_colloc_samples = 1000;  
            delta_samples = 50;
            max_expansion_order = 30; 
            mask(1:17,1) = NaN(17,1);
            mask(1:16,2) = NaN(16,1);
            mask(1:15,3) = NaN(15,1);
            mask(1:14,4) = NaN(14,1);
            mask(1:13,5) = NaN(13,1);
            mask(1:12,6) = NaN(12,1);
            mask(1:11,7) = NaN(11,1);
            mask(1:10,8) = NaN(10,1);
            mask(1:9,9) = NaN(9,1);
            mask(1:8,10) = NaN(8,1);
            mask(1:7,11) = NaN(7,1);
            mask(1:6,12) = NaN(6,1);
            mask(1:5,13) = NaN(5,1);
            mask(1:4,14) = NaN(4,1);
            mask(1:3,15) = NaN(3,1);
            mask(1:2,16) = NaN(2,1);
            mask(1,17) = NaN;
        case 3
            delta_z_mean = 0.01;
            delta_z_var = 0.05;
            max_colloc_samples = 500;
            delta_samples = 50;
            max_expansion_order = 15;
            mask(1:5,1) = NaN(5,1);
            mask(1:4,2) = NaN(4,1);
            mask(1:3,3) = NaN(3,1);
            mask(1:2,4) = NaN(2,1);
            mask(1,5) = NaN;
    end
    mask = flip(mask);
    
    % Grid for surface plot
    x_values = 50:delta_samples:max_colloc_samples;
    x_values = x_values(3:end);                      % discard columns with 50 and 100 samples
    y_values = 1:max_expansion_order;
    
    [X,Y] = meshgrid(x_values,y_values);
    
    epsilon = 0.5;                                 % relative height of Monte Carlo surface for execution time
    
    %% Plot execution time
    % Flat surface of Monte Carlo execution time
    execution_time_monte_carlo = (1+epsilon)*max(execution_time_collocation_max,[],'all')*ones(max_expansion_order,int32(max_colloc_samples/delta_samples)-2);
    % Plot execution times collocation
    figID = figID + 1;
    figure(figID)
    surf(X,Y,execution_time_collocation_max.*mask,'FaceAlpha',0.35);
    hold on;
    zt = round(sort([0:250:max(execution_time_collocation_max,[],'all'),max(execution_time_collocation_max,[],'all'),(1+epsilon)*max(execution_time_collocation_max,[],'all')]),2,'significant');
    zticks(zt);
    zticklabels([zt(1:end-1),execution_time_monte_carlo_value]);
    surf(X,Y,execution_time_collocation_min.*mask,'EdgeColor','none');
    surf(X,Y,execution_time_monte_carlo.*mask,'FaceAlpha',0.7,'EdgeColor','none','LineStyle','none','FaceColor',[0.7 0.7 0.7]);
    xlabel('Collocation samples'); ylabel('Expansion order'); zlabel('Execution time \tau_C [s]');
    title('gPC collocation approach');
    clim([min(execution_time_collocation_min.*mask,[],'all'),max(execution_time_collocation_max.*mask,[],'all')]);
    ax = gca;
    ax.FontSize = 25;
    
    % Insert z-axis break
    axes('Position',[.105 .68 .05 .05]);
    px=[1 5];
    py1=[1 2];
    height=1;
    py2=py1+height;
    plot(px,py1,'k','LineWidth',2);hold all;
    plot(px,py2,'k','LineWidth',2);hold all;
    fill([px flip(px)],[py1 flip(py2)],'w','EdgeColor','none');
    box off;
    axis off;
    hold off;
    
    %% Plot RMSE mean
    % Flat surface of Monte Carlo RMSE for mean
    RMSE_mean_monte_carlo = RMSE_mean_monte_carlo_value*ones(max_expansion_order,int32(max_colloc_samples/delta_samples)-2);
    % Plot RMSE of mean with collocation
    figID = figID + 1;
    figure(figID)
    surf(X,Y,rmse_mean_max.*mask,'FaceAlpha',0.35);
    hold on;
    zt = sort([0:delta_z_mean:max(rmse_mean_max.*mask,[],'all'),RMSE_mean_monte_carlo_value]);
    zt = zt(2:end);     % eliminate zero tick
    zticks(round(zt,2,'significant'));
    zticklabels(round(zt,2,'significant'));
    surf(X,Y,rmse_mean_min.*mask,'EdgeColor','none');
    if model == 2
        surf(X,Y,RMSE_mean_monte_carlo.*mask,'FaceAlpha',0.5,'EdgeColor','none','LineStyle','none','FaceColor',[0.7 0.7 0.7]);
    else
        surf(X,Y,RMSE_mean_monte_carlo.*mask,'FaceAlpha',1,'EdgeColor','none','LineStyle','none','FaceColor',[0.7 0.7 0.7]);
    end
    xlabel('Collocation samples'); ylabel('Expansion order'); 
    if model == 1
        zlabel('RMSE mean (x)');
    elseif model == 2
        zlabel('RMSE mean (y_2)');
    else
        zlabel('RMSE mean (x_1)');
    end
    title('gPC collocation approach');
    clim([min(rmse_mean_min.*mask,[],'all'),max(rmse_mean_max.*mask,[],'all')]);
    set(gca,'Xdir','reverse','Ydir','reverse');
    ax = gca;
    ax.FontSize = 25;
    
    %% Plot RMSE variance
    % Flat surface of Monte Carlo RMSE for variance
    RMSE_var_monte_carlo = RMSE_var_monte_carlo_value*ones(max_expansion_order,int32(max_colloc_samples/delta_samples)-2);
    % Plot RMSE of variance with collocation
    figID = figID + 1;
    figure(figID)
    surf(X,Y,rmse_variance_max.*mask,'FaceAlpha',0.35);
    hold on;
    zt = sort([0:delta_z_var:max(rmse_variance_max.*mask,[],'all'),RMSE_var_monte_carlo_value,min(rmse_variance_min,[],'all')]);
    zt = zt(2:end);     % eliminate zero tick
    zticks(round(zt,2,'significant'));
    zticklabels(round(zt,2,'significant'));
    surf(X,Y,rmse_variance_min.*mask,'EdgeColor','none');
    if model == 2
        surf(X,Y,RMSE_var_monte_carlo.*mask,'FaceAlpha',0.5,'EdgeColor','none','LineStyle','none','FaceColor',[0.7 0.7 0.7]);
    else
        surf(X,Y,RMSE_var_monte_carlo.*mask,'FaceAlpha',1,'EdgeColor','none','LineStyle','none','FaceColor',[0.7 0.7 0.7]);
    end
    xlabel('Collocation samples'); ylabel('Expansion order'); 
    if model == 1
        zlabel('RMSE mean (x)');
    elseif model == 2
        zlabel('RMSE mean (y_2)');
    else
        zlabel('RMSE mean (x_1)');
    end
    title('gPC collocation approach');
    set(gca,'Xdir','reverse','Ydir','reverse');
    clim([min(rmse_variance_min.*mask,[],'all'),max(rmse_variance_max.*mask,[],'all')]);
    ax = gca;
    ax.FontSize = 25;

end
