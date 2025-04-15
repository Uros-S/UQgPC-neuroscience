function figID = plot_galerkin_results(galerkin_data,MC_data,figID)
% This function plots the accuracy results for gPC Galerkin (Hindmarsh-Rose
% model only).

    % Extract gPC collocation data
    accuracy_results_A = galerkin_data{1,1};
    accuracy_results_B = galerkin_data{2,1};
    accuracy_results_C = galerkin_data{3,1};
    accuracy_results_D = galerkin_data{4,1};

    % Extract reference Monte Carlo with 5000 samples
    execution_time_MC = MC_data(1);
    rmse_mean_MC = MC_data(2);
    rmse_var_MC = MC_data(3);

    % Execution time
    figID = figID+1;
    figure(figID);
    plot(accuracy_results_A.execution_time_galerkin,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_B.execution_time_galerkin,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_C.execution_time_galerkin,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_D.execution_time_galerkin,'LineWidth',1.5);
    hold on;
    yline(execution_time_MC,'--','LineWidth',1.5,'Color',[0.1,0.1,0.1])
    grid on;
    set(gca, 'YScale','log');
    title('gPC Galerkin approach');
    xlabel('Expansion order');
    ylabel('Execution time [s]');
    legend('Regime (A)','Regime (B)','Regime (C)','Regime (D)','MC 5000 samples');
    ylim([0,1.2*execution_time_MC])
    ax = gca;
    ax.FontSize = 25;

    % RMSE mean
    figID = figID+1;
    figure(figID);
    plot(accuracy_results_A.rmse_mean,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_B.rmse_mean,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_C.rmse_mean,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_D.rmse_mean,'LineWidth',1.5);
    hold on;
    yline(rmse_mean_MC,'--','LineWidth',1.5,'Color',[0.1,0.1,0.1])
    title('gPC Galerkin approach');
    xlabel('Expansion order');
    ylabel('RMSE mean');
    legend('Regime (A)','Regime (B)','Regime (C)','Regime (D)','MC 5000 samples');
    ax = gca;
    ax.FontSize = 25;
    yl = yticklabels;
    yl{1,1} = num2str(rmse_mean_MC);
    yticks(sort(str2double(yl)));

    % RMSE variance
    figID = figID+1;
    figure(figID);
    plot(accuracy_results_A.rmse_variance,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_B.rmse_variance,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_C.rmse_variance,'LineWidth',1.5);
    hold on;
    plot(accuracy_results_D.rmse_variance,'LineWidth',1.5);
    hold on;
    yline(rmse_var_MC,'--','LineWidth',1.5,'Color',[0.1,0.1,0.1])
    title('gPC Galerkin approach');
    xlabel('Expansion order');
    ylabel('RMSE variance');
    legend('Regime (A)','Regime (B)','Regime (C)','Regime (D)','MC 5000 samples');
    ax = gca;
    ax.FontSize = 25;
    yl = yticklabels;
    yl{1,1} = num2str(rmse_var_MC);
    yticks(sort(str2double(yl)));

end
