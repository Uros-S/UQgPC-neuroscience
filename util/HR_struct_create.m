function [states,parameters,inputs] = HR_struct_create(nom_parameters,x_0,uncertainty)
    % This function creates the structs needed for the PoCET toolbox for
    % gPC expansions. First the nominal system is created and then uniform
    % uncertainty is added. 

    parameters_of_interest = uncertainty{1,1};
    uncertainty_data = uncertainty{1,2};
    
    % Creation of structs for nominal system
    states(1).name = 'x';
    states(1).dist = 'none';
    states(1).data = x_0(1);
    states(1).rhs = '-a*x^3+b*x^2+y-z+I';
            
    states(2).name = 'y';
    states(2).dist = 'none';
    states(2).data = x_0(2);
    states(2).rhs = 'c-d*x^2-y';
            
    states(3).name = 'z';
    states(3).dist = 'none';
    states(3).data = x_0(3);
    states(3).rhs = 'r*(s*(x-x_R)-z)';


            
    parameters(1).name = 'a';
    parameters(1).dist = 'none';
    parameters(1).data = nom_parameters(1);
            
    parameters(2).name = 'b';
    parameters(2).dist = 'none';
    parameters(2).data = nom_parameters(2);
            
    parameters(3).name = 'c';
    parameters(3).dist = 'none';
    parameters(3).data = nom_parameters(3);
            
    parameters(4).name = 'd';
    parameters(4).dist = 'none';
    parameters(4).data = nom_parameters(4);
            
    parameters(5).name = 'I';
    parameters(5).dist = 'none';
    parameters(5).data = nom_parameters(5);
            
    parameters(6).name = 'r';
    parameters(6).dist = 'none';
    parameters(6).data = nom_parameters(6);
            
    parameters(7).name = 's';
    parameters(7).dist = 'none';
    parameters(7).data = nom_parameters(7);
            
    parameters(8).name = 'x_R';
    parameters(8).dist = 'none';
    parameters(8).data = nom_parameters(8);
    
    
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