function [pressure,rho] = isotherm(N,L,T,b)
%% create an isotherm for a phase diagram diagram (P Vs. rho) with a
%% monte carlo simulation, 2D, periodic boundary conditions, Lennad-Jonse
%% interaction.
%       input:
%       N - number of particles
%       L - a row\column vectors with all box sizes
%       T - temprature in reduced units: 4*T[kelvin*bolzman factor]/a
%       b - prameter for the Lennard-Jones potantial: 
%               U(r) = a*[(b/r)^12 -(b/r)^6]

%       output:
%       pressure, rho (density) on the isotherm T. both are reduced:
%       Preduced = 4*p*b^2/a, rhoReduced = rho*b^2

%       this function uses the function "my_num2str(num)" that turns every
%       dot to '_'. ths function is as follows:
        % function str = my_num2str(num)
        % % the number will be converted into a string.
        % % if there is a dot in the number, it will be turned into _
        % 
        %     str = num2str(num);
        %     dotPlace = find('.' == str);
        %     if ~isempty(dotPlace)
        %         str(dotPlace) = '_';
        %     end
        %     
        % end


       len = length(L);
       pressure = zeros(1,len);
       rho = N./(L.^2);
       rho = rho*b^2; %make the density in reduced units
       
       for i = 1:len
           results = MonteCarlo2D('T',T,'Nsteps',10^6,'savMAT'...
               ,[],'N',N,'L',L(i),'pressure',true);
           fileName = ['results_T' my_num2str(T) '_N' num2str(N)...
               '_rho' my_num2str(rho(i))];
           save(fileName);
           pressure(i) = mean(results.pressure);
       end
       
       save(['isotherm_T' my_num2str(T) '_N' num2str(N)]);
end