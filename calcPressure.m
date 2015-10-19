function pressure = calcPressure(allDist,L,n,m,a,b,T)

%% calculate the pressure of a system of 2D particles with pair potantial
%       using the distances between particles 'allDist',
%       calculate the pressure. according to the virial theoram:

%        ==================================
%       | PV = NT - (1/2)sum(r_ij (dU/dr)) |
%        ==================================

%       where r_ij is the distance between particles i,j. the sum is on all
%       particle pairs, each pair is counted once. U is the pair energy:
%       U = a*[(b/r_ij)^n - (b/r_ij)^m]
%       so that: (dU/dr) = a*[-n(b/r_ij)^n/r_ij + m(b/r_ij)^m/r_ij] 
%       finaly:
%%                P = NT/V - (1/2V)sum (-n(b/r_ij)^n + m(b/r_ij)^m) 
%       for more information: http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
%
%       input: 
%       allDist - a matrix containing all pair distances for each step in a
%                 monte carlo simulation, created with 'MonteCarlo2D'
%       L - box length (default: 40)*% 
%       n - the pair potential power law in the distance 
%           (pair potential is u = a*[(b/r)^n - (b/r)^m] 
%           where r is the pair distance and a is some constant)
%       m - the 'm' constant in the pair potantial
%       a - the 'a' constant in the pair potatial
%       b - the 'b' constant in the pair potantial 
%       T - Temperature 

        V = L^2;
        
       [N,~,steps] = size(allDist); 
       P1 = allDist;
       P2 = allDist;
       
       a1 = 'calc P1';
       b1 = 'calc P2';
       
       for i = 1:steps
            i
            % P1 = (1/r_ij)^(-n)
            a1
            P1(:,:,i) = tril(allDist(:,:,i).^(-n),-1);
            % P2 = (1/r_ij)^(-m)
            b1
            P2(:,:,i) = tril(allDist(:,:,i).^(-m),-1);
       end
       
       % P = NT/V - (1/2V)sum (-n(b/r_ij)^n + m(b/r_ij)^m)   
       % P = NT/V - (a/2V)sum (-n b^n P1 + m b^m P2)
       c = 'calc pr';
       c
       pr = (N*T/V) - (a/(2*V))*(sum(sum(-n*(b^n)*P1 + m*(b^m)*P2)));
       c = 'done calc pr';
       c
       % make pr a row vector
       pressure(1,:) = pr(1,1,:);

end