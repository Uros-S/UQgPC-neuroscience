function [states,parameters,inputs] = JR_struct_create(nom_parameters,x_0,uncertainty)
    % This function creates the structs needed for the PoCET toolbox for
    % gPC expansions. First the nominal system is created and then uniform
    % uncertainty is added. 

    parameters_of_interest = uncertainty{1,1};
    uncertainty_data = uncertainty{1,2};
    
    % Creation of structs for nominal system
    states(1).name = 'y_1';
    states(1).dist = 'none';
    states(1).data = x_0(1);
    states(1).rhs = 'y_4';
            
    states(2).name = 'y_2';
    states(2).dist = 'none';
    states(2).data = x_0(2);
    states(2).rhs = 'y_5';
            
    states(3).name = 'y_3';
    states(3).dist = 'none';
    states(3).data = x_0(3);
    states(3).rhs = 'y_6';

    states(4).name = 'y_4';
    states(4).dist = 'none';
    states(4).data = x_0(4);
    states(4).rhs = 'A*a*(nu_max/(1+exp(r*(V_0-(y_2-y_3))))) - 2*a*y_4 - a^2*y_1';

    states(5).name = 'y_5';
    states(5).dist = 'none';
    states(5).data = x_0(5);
    states(5).rhs = 'A*a*(p + 0.8*C*(nu_max/(1+exp(r*(V_0-(C*y_1)))))) - 2*a*y_5 - a^2*y_2';

    states(6).name = 'y_6';
    states(6).dist = 'none';
    states(6).data = x_0(6);
    states(6).rhs = 'B*b*0.25*C*(nu_max/(1+exp(r*(V_0-(0.25*C*y_1))))) - 2*b*y_6 - b^2*y_3';


            
    parameters(1).name = 'C';
    parameters(1).dist = 'none';
    parameters(1).data = nom_parameters(1);
            
    parameters(2).name = 'A';
    parameters(2).dist = 'none';
    parameters(2).data = nom_parameters(2);
            
    parameters(3).name = 'B';
    parameters(3).dist = 'none';
    parameters(3).data = nom_parameters(3);
            
    parameters(4).name = 'a';
    parameters(4).dist = 'none';
    parameters(4).data = nom_parameters(4);
            
    parameters(5).name = 'b';
    parameters(5).dist = 'none';
    parameters(5).data = nom_parameters(5);
            
    parameters(6).name = 'V_0';
    parameters(6).dist = 'none';
    parameters(6).data = nom_parameters(6);
            
    parameters(7).name = 'nu_max';
    parameters(7).dist = 'none';
    parameters(7).data = nom_parameters(7);
            
    parameters(8).name = 'r';
    parameters(8).dist = 'none';
    parameters(8).data = nom_parameters(8);

    parameters(9).name = 'p';
    parameters(9).dist = 'none';
    parameters(9).data = nom_parameters(9);
    
    
    inputs = [];

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