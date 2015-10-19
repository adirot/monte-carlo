function [bins,histo] = RDF(allDist,L,NumOfBins)
%% calculate 2D radial distribution function, with pediodic boundary
%% conditions.
% this function was built for Monte Carlo calculations with the function 
% 'MonteCarlo2D', but it can be used without it.

%         input:
%         ~~~~~~
%         allDist - distance matrix of all particles. this can be created
%         with the function 'MonteCarlo2D'. 

%         L - size of the board.
%
%         NumOfBins - number of bins in the histogram
%
%         output:
%         ~~~~~~~
%         bins - x axis of the RDF histogram
%         histo - each colmun is the RDF axis of a diffrent step in the
%         monte carlo simulation. (so to plot the RDF for the 10'th step we
%         need to write: plot(bins,histo(:,10));

%       method:
%       ~~~~~~~
%       the RDF is calculated by binnig all pair partical distances into 
%       a histogram, and normalizing each bin with it's Ideal gas number of
%       particals. 
%       when binning the pair distances, we take into acount Periodic
%       Boudary Conditions (PBC)
%       finaly, to ansure that when r->infinity : RDF ->1 , we
%       normalize by multiplying RDF*2/(N-1) where N is the number of
%       particals. 
%       for more information
%   http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
%       page 48 - "Radial distribution function"
           
          [N,~,steps] = size(allDist); 
          histo = zeros(NumOfBins,steps);      
          
          for step = 1:steps
                
                dist = allDist(:,:,step);
                d = reshape(dist,1,[]);
                d = nonzeros(d);
                d = d(d < L*0.5);
                bins = linspace(0,L*0.5,NumOfBins);
                histo(:,step) = hist(d,bins);
                increment = bins(2) - bins(1);
                rho = N/(L^2);
                
                % each bin should be normalized according to its volume
                for bin = 1:NumOfBins
                    
                        % rVal is the number of particles in some layer of area 
                        % da(r)=2pi*r*dr, a distance r from the central cell
                        rVal = bins(bin);
                        next_rVal = increment + rVal;
                        
                        % Calculate the area of the bin (a ring of radii r,
                        % r+dr)
                        ereaBin = pi*next_rVal^2 - pi*rVal^2;
                      
                        % Calculate the number of particles expected in this bin in
                        % the ideal case
                        nIdeal = ereaBin*rho;
                        
                        % Normalize the bin
                        histo(bin,step) = histo(bin,step) / nIdeal;
                        
                end
                    histo(:,step) = 2*histo(:,step)/(N-1);
          end

        end