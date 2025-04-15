function [states,parameters,inputs] = Epileptor_struct_create(nom_parameters,x_0,uncertainty)
    % This function creates the structs needed for the PoCET toolbox for
    % gPC expansions. First the nominal system is created and then uniform
    % uncertainty is added. 

    parameters_of_interest = uncertainty{1,1};
    uncertainty_data = uncertainty{1,2};
    
    % Creation of structs for nominal system
    states(1).name = 'x_1';
    states(1).dist = 'none';
    states(1).data = x_0(1);
    states(1).rhs = 'x_2 - f_1 - x_3 + I_1';
            
    states(2).name = 'x_2';
    states(2).dist = 'none';
    states(2).data = x_0(2);
    states(2).rhs = 'r_2 - 5*x_1^2 - x_2';
            
    states(3).name = 'x_3';
    states(3).dist = 'none';
    states(3).data = x_0(3);
    states(3).rhs = '1/tau_0*(4*(x_1 - r_1) - x_3)';

    states(4).name = 'x_4';
    states(4).dist = 'none';
    states(4).data = x_0(4);
    states(4).rhs = '-x_5 + x_4 -x_4^3 + I_2 + 2*u - 0.3*(x_3 - 3.5)';

    states(5).name = 'x_5';
    states(5).dist = 'none';
    states(5).data = x_0(5);
    states(5).rhs = '1/tau_2*(-x_5 + f_2)';

    states(6).name = 'u';
    states(6).dist = 'none';
    states(6).data = x_0(6);
    states(6).rhs = '-gamma*(u-0.1*x_1)';
            

    parameters(1).name = 'r_1';
    parameters(1).dist = 'none';
    parameters(1).data = nom_parameters(1);
            
    parameters(2).name = 'r_2';
    parameters(2).dist = 'none';
    parameters(2).data = nom_parameters(2);
            
    parameters(3).name = 'I_1';
    parameters(3).dist = 'none';
    parameters(3).data = nom_parameters(3);
            
    parameters(4).name = 'I_2';
    parameters(4).dist = 'none';
    parameters(4).data = nom_parameters(4);
            
    parameters(5).name = 'tau_0';
    parameters(5).dist = 'none';
    parameters(5).data = nom_parameters(5);
            
    parameters(6).name = 'm';
    parameters(6).dist = 'none';
    parameters(6).data = nom_parameters(6);
            
    parameters(7).name = 'tau_2';
    parameters(7).dist = 'none';
    parameters(7).data = nom_parameters(7);
            
    parameters(8).name = 'gamma';
    parameters(8).dist = 'none';
    parameters(8).data = nom_parameters(8);
    

    inputs(1).name = 'f_1';
    inputs(1).rhs = 'piecewise(x_1<0, x_1^3-3*x_1^2, x_1>=0, (x_4-0.6*(x_3-4)^2-m)*x_1)';
    
    inputs(2).name = 'f_2';
    inputs(2).rhs = 'piecewise(x_4<-0.25, 0, x_4>=-0.25, 6*(x_4+0.25))';
    
    % Addition of uncertainty on parameter of interest
    for i=1:length(parameters_of_interest)
        if uncertainty{1,3} == 1
                parameters(parameters_of_interest(i)).dist = 'uniform';
                parameters(parameters_of_interest(i)).data = [nom_parameters(parameters_of_interest(i)), ...
                                                  nom_parameters(parameters_of_interest(i))+uncertainty_data(i)];
        elseif uncertainty{1,3} == 2
                parameters(parameters_of_interest(i)).dist = 'uniform';
                parameters(parameters_of_interest(i)).data = [nom_parameters(parameters_of_interest(i))-uncertainty_data(i)/2, ...
                                                  nom_parameters(parameters_of_interest(i))+uncertainty_data(i)/2];
        elseif uncertainty{1,3} == 3
            parameters(parameters_of_interest(i)).dist = 'uniform';
            parameters(parameters_of_interest(i)).data = [nom_parameters(parameters_of_interest(i))-uncertainty_data(i), ...
                                                  nom_parameters(parameters_of_interest(i))];
        elseif uncertainty{1,3} == 4
            parameters(parameters_of_interest(i)).dist = 'normal';
            parameters(parameters_of_interest(i)).data = [nom_parameters(parameters_of_interest(i)),uncertainty_data(i)];
        end
    end

end