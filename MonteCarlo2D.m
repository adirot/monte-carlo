function results = MonteCarlo2D(varargin)
        
        %% 3 %% in this version I will:
        


        %% Monte Carlo in 2D for Lennard-Jones like potential with hard 
        %% discs, NVT ansamble, PBC (piriodic boundary conditions)
        
        % this code and all the functions it uses are in the git
        % repository: https://github.com/adirot/monte-carlo.git
        % some general functions this code uses are in: 
        % https://github.com/adirot/general-matlab-functions.git
        
        % This function runs a Monte Carlo simulation, according to the
        % Metropolis algorithm. for more information about Monte Carlo and
        % Metropolis see: Allen, Tildesley computer simulation of liquids,
        % section 4.4
        
        % we sett a potential between pairs of particals:
        % u = 4*epsilon*[(sigma/r)^n - (sigma/r)^m]
        % when n = 12, m = 6 This is the Lennard-Jonse potential.
        % for more information: https://en.wikipedia.org/wiki/Lennard-Jones_potential
        
        % In reduced units:
        %       U' = 4*[(1/r')^n - (1/r')^m]
        %       T' = kB*T/epsilon
        %       r' = r/sigma
        %       U' = U/epsilon
        %       P' = P*sigma^2/epsilon
        
        % The particals in this simulation have a radius sett by the user,
        % and we assume hard core repultion (the particals cannot overlap).
        
        % This code can also compute the Radial distribution function. for
        % more information about RDF, see
        % http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
        % page 48 - "Radial distribution function"
        
        
        %% Particles system parameters
        % N - number of particles (default: 256)
        % L - box length in reduced units, user cannot set this
        % rho - density of the particles in reduced units. (default: 0.3)
        % n - the pair potential power law in the distance 
        %     (pair potential is u = 4*epsilon*[(sigma/r)^n - (sigma/r)^m] 
        %     where r is the pair distance and a is some constant)
        %     (default: 12)
        % m - the 'm' constant in the pair potantial (default: 6)
        % T - reduced Temperature (default: 0.6)
        % r - partical radius in reduced units: allways (2^(-1/6))/2, 
        %       user cannot sett this
        
        % important note: if you change the potantial(n,m parameters) you
        % might get a long range potantial. in this case you should
        % consider changing the truncation of the energy calculation. 
        
        %% Simulation parameters
        % Nsteps - number of Monte carlo steps (default: 500000)
        % dr - max displacement of a particle in each step, this is only an
        %      initial value that will change during the run.
        %       reduced units.
        %      (default: 0.5)*
        % rCutoff - largest r for which we calculate the energy. reduced
        %           units.
        %           (recomended: 2.5)
        % sampleFreq - the frequency of sampling (for ensemlse average)
        %               (recommended - 5*N to 10*N, set to 5*N)
        % optimize_dr - true of false. should dr be optimized during run or
        %               kept at it's initial value (see the function
        %               'optimizedr' below for details on how optimization 
        %               is done. notice that for some combinations of
        %               tempratures and densities optimization is not 
        %               recommended, the function should stop and alart you
        %               to set optimize_dr to false in those cases.
        % initialConfig - options: 'random' , 'hex' or coordinate list.
        %                 'random': random initial configuration with no
        %                 ovelaps. not recomanded for high densities!
        %                 'hex' - hexagonal packing
        %                 2 by N vector - initial coordinates to use.
        %
        % * if we choose a dr that is too high or too low, our simulation
        %   will converge slowly. it is customary to make sure that dr is
        %   chosen so that half of the steps are kept. for this end we
        %   check the acceptance rate every 100 steps, and change dr 
        %   accordingly. for more information: Allen, Tildesley computer
        %   simulation of liquids, section 1.5.2
        
        
        %% RDF parameters
        % NumOfBins - nubers of bins in RDF histogram (default: 7).
        %           if you sett NumOfBins = [], RDF will not be calculated.
        
        %% Pressure calcuation
        % pressure - true or false, calculate the pressure or not.(defaut:
        %           false). reduced units
        
        %% Display save and plot options
        % savMAT - saved file name (mat file). will be followed by a number
        %           noteing the step number saved in this file.
        %           if empty, the result will not be saved
        % maxVarSize - save to mat file every time a variable reaches 'maxVarSize'
        %           (default: 1000)
        
        %% The algorithm:
        % 1. choose starting possitions of the particals. 
        % 
        % 2. Monte Carlo, Metropolis:
        %       a. calculate the energy in reduced units.
        %       b. try to move a partical: 
        %          choose a random particle. randomly
        %          choose a displacement within a box of size dr^2 aroud
        %          the chosen particle.
        %       c. calculate dU - the change of energy in reduced units. 
        %          (during this calculation we check for ovarlaps, if we
        %          find an overlap the step is rejected)
        %       d. if dU < 0 - keep step (save the new energy U + dU).
        %       e. if dU > 0 - keep step with probability exp(-dU/T). 
        %           both dU and T are in reduced units.
        %           importante: if we rejecte a step, we save the old 
        %           energy.
        %       
        % 3. calculate Radial distribution function (RDF):
        %       using the distances between particals already found,
        %       calculate the RDF.
        %
        % 4. calculate the pressure in reduced units: 
        %       using the distances between particles already found,
        %       calculate the pressure. according to the virial theoram:
        %       P' = rho'T' - (2*rho'/N)*sum(6(1/r'_ij)-12(1/r'_ij)), 
        %       where r'_ij is the distance between particles i,j in 
        %       reduced units. the sum is on all particle 
        %       pairs, each pair is counted once. (see the note about 
        %       reduced units in this git repository for more information:
        %       https://github.com/adirot/monte-carlo.git
        %
        % 5. save simulation parameters every 5*N steps.
        %
        % 6. accumulate averages: energy, RDF, pressure.
        %
        % 7. calculate the energy of the final configuration, just to make
        %    sure we did nothing wrong.
        
        %% The Verelet nieghbour list
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Check input arguments %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
            p = inputParser();

            %% system parameters
            addOptional(p, 'N', 256);
            addOptional(p, 'rho', 0.3);
            addOptional(p, 'Un', 12);
            addOptional(p, 'Um', 6); 
            addOptional(p, 'T', 0.6);
            r = (2^(-1/6))/2; %particle radius

            %% simulation parameters
            addOptional(p, 'Nsteps', 500000);
            addOptional(p, 'dr', 0.5);
            addOptional(p, 'rCutoff', 2.5);
            addOptional(p, 'sampleFreq', 5);
            addOptional(p, 'optimize_dr', true);
            addOptional(p, 'initialConfig', 'random');
            addOptional(p, 'startFromOldRun_fileName', []);
            
            
            %% RDF parameters
            addOptional(p, 'NumOfBins', 7);
            
            %% pressure calculation
            addOptional(p, 'pressure', false);
            
            %% display save and plot options
            addOptional(p, 'savMAT', 'outputMC-'); % saved file name
                        % if empty, the result will not be saved
            addOptional(p, 'maxVarSize', 1000); % this option will save the 
            % variables when they have reached the spacified size, and
            % clear the workspace for better preformance.

            parse(p, varargin{:});
            Results = p.Results;

            N = Results.N;
            rho = Results.rho;
            L = sqrt(N/rho);
            n = Results.Un;
            m = Results.Um;
            T = Results.T; 
            Nsteps = Results.Nsteps;
            dr = Results.dr;
            init_dr = dr;
            sampleFreq = Results.sampleFreq;
            NumOfBins = Results.NumOfBins;
            optimize_dr = Results.optimize_dr;
            initialConfig = Results.initialConfig;
            startFromOldRun_fileName = Results.startFromOldRun_fileName;
            pressure = Results.pressure;
            savMAT = Results.savMAT;
            maxVarSize = Results.maxVarSize;
            rCutoff = Results.rCutoff;
            
            results.N = Results.N;
            results.rho = Results.rho;
            results.L = L;
            results.n = Results.Un;
            results.m = Results.Um;
            results.T = Results.T;
            results.Nsteps = Results.Nsteps;
            results.dr = Results.dr;
            results.initialConfig = Results.initialConfig;
            results.sampleFreq = Results.sampleFreq;
            results.NumOfBins = Results.NumOfBins;
            results.pressure = Results.pressure;
            results.savMAT = Results.savMAT;
            results.maxVarSize = Results.maxVarSize;
            results.rCutoff = Results.rCutoff;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Set start possitions of N particals %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % also calculate all pair distances, in PBC
        
        if ischar(initialConfig)
            if strcmp(initialConfig,'random')
                % create a random initial configuration
                particlesPosition = randomStart;
                if isempty(particlesPosition)
                     results = [];
                     return
                end

            else
                if strcmp(initialConfig,'hex')
                    % create Hexagonal close packed initial configuration
                    particlesPosition = hcp;
                     if isempty(particlesPosition)
                         results = [];
                         return
                     end
                else
                    if ~isempty(startFromOldRun_fileName)
                        
                        % get old run data
                        load(startFromOldRun_fileName);
                        N = results.N;
                        rho = results.rho;
                        L = sqrt(N/rho);
                        n = results.Un;
                        m = results.Um;
                        T = results.T; 
                        dr = results.dr;
                        init_dr = dr;
                        sampleFreq = results.sampleFreq;
                        NumOfBins = results.NumOfBins;
                        optimize_dr = results.optimize_dr;
                        initialConfig = results.initialConfig;
                        pressure = results.pressure;
                        savMAT = results.savMAT;
                        maxVarSize = results.maxVarSize;
                        rCutoff = results.rCutoff;
                        
                        % get last particles configuration
                        particlesPosition = particlesPosition(:,:,end);
                        startStep = results.Nsteps + 1;
                        Nsteps = results.Nsteps + Nsteps;
                    end
                        
                end
            end
        else
            [xsize,ysize] = size(initialConfig);
            if xsize == 2 && ysize == N
                particlesPosition = initialConfig;
            else
                if xsize == N && ysize == 2
                    particlesPosition = initialConfig';
                else
                    results = [];
                    disp('initialConfig must be "hex" , "random" or a coordiantes list');
                    return
                end
            end
        end
                  

        
        
        % calculate all pair distances
        dist = zeros(N);
        for par = 1:N
            
            x = particlesPosition(1,par);
            y = particlesPosition(2,par);
            dist((par+1):N,par) = ...
                distPBC(x,y,particlesPosition(:,(par+1):N));
        end
        
            % you can plot the particles if you want, to see 
            % the initial configuration:
            
            % plotParticles(particlesPosition,L,r)
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Run Monte Carlo simulation %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
            %% calculate the initial energy
            d = reshape(dist,1,[]);
            d = nonzeros(d);
            U = pairU(d);

            %% initialize parameters
            if ~isempty(startFromOldRun_fileName)
                moveCount = 0; % counts eccepted steps
                fileNameList = {};
            end
            allU = zeros(1,maxVarSize); %the energy in every sample
            allDist = zeros(N,N,maxVarSize);%the distances in every sample
            allCoords = zeros(2,N,maxVarSize);%the coordinates in every sample
            t = 0; % allU index
            
            
            ec1 = 0 ; ec2 = 0; ec3 = 0; 
            displace = [];
            displacedInd = [];
            
            %% start Monte Carlo
            for step = startStep:Nsteps
                    step
                    %% accumulate averages, save distances and coordinates
                    if  mod(step,sampleFreq*N) == 0
                        
                                t = t + 1;
                                allU(1,t) = U;
                                allDist(:,:,t) = dist;
                                allCoords(:,:,t) = particlesPosition; 
                                %step
                                
                    end
       
                    %% move a particale (don't allow overlaps,set PBC)
                                
                            % chose particle to move
                            movedParticle = randi([1 N]);
                            
                            % choose displacement:
                            displacex = dr*rand - (dr/2);
                            displacey = dr*rand - (dr/2);

                            % move particle
                            newParticlesPosition = movePBC;

                            % calculate new distances
                            newDist = reCalcDist;
                            
                    % make sure there are no overlaps
                    overlap = checkOverlap;
                    
                    if ~overlap
                        %% calculate the change in energy after movement
                        dU = Uchange;

                        %% if dU < 0 eccept move
                        if dU < 0  
                                U = U + dU;
                                dist = newDist;
                                particlesPosition = newParticlesPosition;
                                moveCount = moveCount + 1;
                                %step
                                ec1 = ec1 + 1; 
                                displace = [displace...
                                    sqrt(displacex^2 + displacey^2)];
                                displacedInd = [displacedInd movedParticle];
                        else
                                %% otherwise,
                                % keep the new state with a 
                                % probability corresponding to the Boltzmann
                                % factor. if the new state is rejected,
                                % recount the old configuration.
                                if rand < exp(-(1/T)*dU)
                                    U = U + dU;
                                    dist = newDist;
                                    particlesPosition = newParticlesPosition;
                                    moveCount = moveCount + 1;
                                 %   step
                                    
                                    displace = [displace...
                                    sqrt(displacex^2 + displacey^2)];
                                    displacedInd = [displacedInd movedParticle];
                                    ec2 = ec2 + 1;
                                else
                                    ec3 = ec3 + 1;
                                end
                        end
                    end
                    
%                     ec1
%                     ec2
%                     ec3
                    %% find the optimal dr (check it every 100 steps)
                    % this find the optimal dr, for a faster run. sometimes
                    % this is not recomended, the function should alart you
                    % if thats the case.
                    
                    %moveCount/step
                    if (mod(step,100) == 0)&&((moveCount/step)~=0.5...
                            &&optimize_dr)
                            
                            dr = optimizedr;
                            if dr > L/2 || dr < L*0.0001
                                % Display an error message
                                disp(dr);
                                results.dr = dr;
                                disp('you should consider turning off the optimize_dr option - dr became too large or too small. the function will tun off optimize_dr now');
                                
                                % automaticly turns off
                                % dr_optimitation:
                                dr = init_dr;
                                optimize_dr = false;
                                results.commet = 'optimize_dr have been turned off'
                                
                            end
                                
                    end
                        
                    %% save to mat file, clear veriables for better 
                    %% preformance, Calculate ensamble avrages
                    if t == maxVarSize
                            
                            t = 0;
                            
                            %% Calculate ensamble avrages %%

                                if ~isempty(NumOfBins)
                                    % calculate RDF: 
                                    % bins is the x axis of the RDF. histo is the y axis, for all
                                    % samples. each  column in "histo" is a sample (so if you want
                                    % to plot sample 10, you use: plot(bins,histo(:,10))

                                    [bins,histo] = RDF(allDist,L,NumOfBins);
                                    results.bins = bins;
                                    results.histo = histo;
                                end

                                if pressure
                                    allPressure = ...
                                        calcPressure(allDist,rho,n,m,T); 
                                    results.allPressure = allPressure;
                                end
                                
                                if ~isempty(savMAT)
                                    saveToMat(step,maxVarSize,savMAT);
                                end
                                
                    end
                    
            end
            
            
            %% calculate the final energy, just for debugging
             d = reshape(dist,1,[]);
             d = nonzeros(d);
             U = pairU(d);
             
            results.allU = allU;
            results.finalU = U;
            results.allDist = allDist;
            results.allCoords = allCoords;
            
            results.displace = displace;
            results.displacedInd = displacedInd;
            results.fileNameList = fileNameList;
            results.lastParticlesPosition = particlesPosition(:,:,end);
            results.moveCount = moveCount;
                
        %% Calculate ensamble avrages %%
        
        len = length(fileNameList);
        
        allU = zeros(1,len);
        for ii = 0:(len-1)
              mat = matfile(fileNameList{1,ii+1});
              allU(1,(1:maxVarSize)+ii*maxVarSize) = mat.allU;
        end
        results.meanU = my_mean(allU);
        
        if ~isempty(NumOfBins)
            % calculate RDF: 
            % bins is the x axis of the RDF. histo is the y axis, for all
            % samples. each  column in "histo" is a sample (so if you want
            % to plot sample 10, you use: plot(bins,histo(:,10))
               
            bins = mat.bins;
            histo = zeros(length(bins),len);
            for ii = 0:(len-1)
                mat = matfile(fileNameList{1,ii+1});
                histo(:,(1:maxVarSize)+ii*maxVarSize) = mat.histo;
            end
            
            % mean RDF
            meanHisto = mean(histo,2);
            
            results.bins = bins;
            results.histo = histo;
            results.meanHisto = meanHisto;
            clear histo;
            
        end
        
        if pressure
            
            allPressure = zeros(1,len);
            for ii = 0:(len-1)
                mat = matfile(fileNameList{1,ii+1});
                allPressure(:,(1:maxVarSize)+ii*maxVarSize)...
                    = mat.allPressure;
            end
            
            results.allPressure = allPressure;
            results.meanPressure = my_mean(allPressure);
            
        end
        
        
        if ~isempty(savMAT)
            save(savMAT);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Functions used in this code %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function particlesPosition = hcp
    
                % Make sure N is a perfect square
                intRoot = floor(sqrt(N));
                if (sqrt(N) - intRoot) > 1e-7
                    % Display an error message
                    disp('Number of particles should be a perfect square');
                    particlesPosition = [];
                    return
                end

                % Calculate the seperation length between particles centers
                sepDist = L/sqrt(N);
                
                % Make sure the density is not too high
                if sepDist < 2*r
                    % Display an error message
                    disp('density is too high');
                    particlesPosition = [];
                    return
                end

                % Find the box size
                Lx = sepDist * sqrt(N);

                % Create a vector of linearly spaced points along the
                % x-direction
                xPos = linspace(sepDist/2, Lx-sepDist/2, sqrt(N));
                % And find the corresponsing y-direction increments
                yPos = (sqrt(3)/2)*xPos;

                % Create a matrix with all combinations of x and y
                [X,Y] = meshgrid(xPos,yPos);
                % Shift coordinates to the be at the center of each
                % particle
                X(1:2:end,:) = X(1:2:end,:) + sepDist/2;

                % Reshape the matrix to be 1D in X and 1D in Y
                % (numel returns the number of elements in a given array)
                particlesPosition =...
                    [reshape(X,1,numel(X));reshape(Y,1,numel(Y))];
                
                % make the board in: [-L/2 L/2]x[-L/2 L/2]
                particlesPosition = particlesPosition - L/2;

        end

        
        function particlesPosition = randomStart
                % randomize first particle possition in the box 
                % [-L/2,L/2] x [-L/2,L/2]
                particlesPosition(1,1) = L*rand - (L/2);
                particlesPosition(2,1) = L*rand - (L/2);
                dist = zeros(N);

                for j = 2:N
                    
                      % choose random possition
                      particlesPosition(1,j) = L*rand - (L/2);
                      particlesPosition(2,j) = L*rand - (L/2);

                      % calculate PBC distances
                      xj = particlesPosition(1,j);
                      yj = particlesPosition(2,j);
                      dist(j,1:j) = distPBC(xj,yj,particlesPosition);

                      % check for piriodic boundary condition overlaps,
                      % randomize new possition if overlaps are found.
                      overlapPBC = sum(dist(j,1:(j-1)) < 2*r) > 0;
                      while overlapPBC
                              particlesPosition(1,j) = L*rand - (L/2);
                              particlesPosition(2,j) = L*rand - (L/2);

                              % calculate PBC distances
                              xj = particlesPosition(1,j);
                              yj = particlesPosition(2,j);
                              dist(j,1:j) = distPBC(xj,yj,particlesPosition);

                              overlapPBC = sum(dist(j,1:(j-1)) < 2*r) > 0;
                      end

                end

       end
        
        function dist = distPBC(x,y,allPossition)
            % calculate the PBC distances of the point (x,y) from all
            % particle possitions
            
                distx = abs(x - allPossition(1,:));
                disty = abs(y - allPossition(2,:));
                bigDistx = distx > (L/2);
                bigDisty = disty > (L/2);
                distx = distx - L.*bigDistx;
                disty = disty - L.*bigDisty;
                dist = sqrt(distx.^2 + disty.^2);
        end
    
            
        function newparticlesPosition = movePBC
            % if the particle gets out of the board, apply PBC
                newparticlesPosition = particlesPosition;
                x = particlesPosition(1,movedParticle)+ displacex;
                y = particlesPosition(2,movedParticle)+ displacey;
                
                if x > L/2
                    x = x - L;
                end
                
                if x < (-L/2) 
                    x =  x + L;
                end
                
                if y > L/2 
                    y = y - L;
                end
                
                if y < (-L/2) 
                    y = y + L;
                end
                
                newparticlesPosition(1,movedParticle) = x;
                newparticlesPosition(2,movedParticle) = y;
        end
        
        function newdist = reCalcDist
            
            % recalculates pair distances after moving a particle
            
                xi = newParticlesPosition(1,movedParticle);
                yi = newParticlesPosition(2,movedParticle);
                newdist = dist;
                
                % recalculate the relevent row elements in dist matrix
                if movedParticle > 1
                    newdist(movedParticle,1:(movedParticle-1)) =...
                        distPBC(xi,yi,...
                            particlesPosition(:,1:(movedParticle-1)));
                end
                
                % recalculate the relevent column elements in dist matrix
                if movedParticle < N
                    newdist((movedParticle + 1):N,movedParticle) =...
                        distPBC(xi,yi,...
                            particlesPosition(:,(movedParticle+1):N));
                end
        end
    
        function overlap = checkOverlap
            % checking for overlaps in the newDist matrix
                
                overlap = false;
                
                % check relevent rows
                if movedParticle > 1
                    row = ...
                        sum(newDist(movedParticle,1:(movedParticle-1))<...
                            (2*r)) > 0;
                else
                    row = false;
                end
                
                if row
                   overlap = true;
                else
                   % check relevent columns
                   if movedParticle < N
                     col = sum(newDist((movedParticle+1):N,movedParticle)...
                       < (2*r)) > 0;
                   else
                       col = false;
                   end
                   
                   if col
                       overlap = true;
                   end
                   
                end
        end
    
        function U = pairU(dist)
            % input is a row vector of all pair distances
            % calculates the reduced energy according to the pair potantial.
            % U is the total energy, u are the energies of each pair
            
                dist1 = dist(dist < rCutoff);
                u = 4*(((1./dist1).^n)-((1./dist1).^m));
                U = sum(u);
        end
    
        function dU = Uchange
            % only calculates the change in energy after a particle has
            % moved
                
                % calculate oldU for the relevant row
                if movedParticle > 1
                    oldUrow = ...
                        pairU(dist(movedParticle,1:(movedParticle - 1)));
                else 
                    oldUrow = 0;
                end
                
                % calculate oldU for the relevant column
                if movedParticle < N
                    oldUcol = ...
                        pairU(dist((movedParticle + 1):N,movedParticle));
                else 
                    oldUcol = 0;
                end
                
                % calculate oldU
                oldU = oldUrow + oldUcol;
                
                % calculate newU for the relevant row
                if movedParticle > 1
                    newUrow = ...
                        pairU(newDist(movedParticle,1:(movedParticle - 1)));
                else 
                    newUrow = 0;
                end
                
                % calculate oldU for the relevant column
                if movedParticle < N
                    newUcol = ...
                        pairU(newDist((movedParticle + 1):N,movedParticle));
                else 
                    newUcol = 0;
                end
                
                newU = newUrow + newUcol;
                
                dU = newU - oldU;
        end
    
        function newdr = optimizedr
            %   it is customary to make sure that dr is
            %   chosen so that half of the steps are kept. we fix dr
            %   accordindly:
                   
                    if (moveCount/step) > 0.5
                           newdr = dr*1.05;
                    end

                    if(moveCount/step) < 0.5
                           newdr = 0.95*dr;
                    end

        end

        function saveToMat(step,maxVarSize,savMAT)
            
		    fileInd = (step/(sampleFreq*N))/maxVarSize;		
                    fileName = strcat(savMAT,num2str((fileInd-1)*maxVarSize),'-',...
                                    num2str(fileInd*maxVarSize));
                    mat = matfile(fileName,'Writable',true);
                    fileNameList{1,fileInd} = fileName;
                    mat.results = results;
                    mat.allU = allU;
                    mat.allDist = allDist;
                    mat.allCoords = allCoords;
                    
                    if ~isempty(NumOfBins)
                        mat.bins = results.bins;
                        mat.histo = results.histo;
                        clear histo;
                    end
                    
                    if pressure
                        mat.allPressure = allPressure;
                        clear allPressure;
                    end
        end
        
    function meanProp = my_mean(prop)
        % input: some property of the simulation in all steps.
        % output: meanProp(i) is the mean of all the values of property on
        % steps 1 to i.
            
            lenProp = length(prop);
            meanProp = zeros(1,lenProp); 
            meanProp(1) = prop(1);
            for i = 2:lenProp
                meanProp(i) = meanProp(i-1) + prop(i);
            end
            one2len = 1:lenProp;
            meanProp = meanProp./one2len;
    end
end
