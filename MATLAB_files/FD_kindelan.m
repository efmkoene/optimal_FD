function [kindelan,k_extrema] = FD_kindelan( L, grerr )
% FD_KINDELAN Compute the optimized finite-difference coefficients cf. Kindelan (1990) 
%   FD_kindelan( L, grerr ) finds finite-difference coefficients that
%   approximate first-order derivatives using a least-squares minimization
%   of the group velocity error, within an error bound.
%
%   This function is based on Kindelan, Kamel and Sguazzero (1990): https://library.seg.org/doi/pdf/10.1190/1.1442763
%
%   For staggered derivatives of type 'f(x+1/2)-f(x-1/2)'.
%
% Written by Erik Koene, erik.koene@erdw.ethz.ch.
% ETH Zurich, 15-May-2018. 
% Copyright goes to the respective authors of the papers and their publishers. 

% SET INITIAL CONDITIONS
X0 = [FD_taylor(L)';linspace(0,sqrt(L),L)'];

% SET PROBLEM
problem.objective = @(X) F(L,grerr,X);
problem.x0 = X0;
problem.solver = 'lsqnonlin';
problem.options = optimoptions('lsqnonlin', ...
                               'Algorithm','trust-region-reflective',...
                               'SpecifyObjectiveGradient',true,...
                               'ScaleProblem','jacobian',...
                               'MaxFunEvals',1000,...
                               'MaxIterations',1000,...
                               'FunctionTolerance',eps*10, ...
                               'StepTolerance',eps*10, ...
                               'OptimalityTolerance',eps*10,...
                               'Display','off');
% Set upper and lower boundary for the values
problem.lb = [-2*ones( L,1);
                0;
                0.1*ones(L-1,1)];
problem.ub = [ 2*ones( L,1);
                eps;
               3*ones(L-1,1)];

% SOLVE PROBLEM FIRST
[x,~,flag] = lsqnonlin( problem );

% Solve again now that we are closer to the final solution, making sure k 
% is sorted upwards (0,...,k_max) and no (near) duplicates appear
for i=1:3
problem.x0 = [FD_taylor(L)';linspace(0,max(x(L+1:end)),L)'];
[x,~,flag] = lsqnonlin( problem );
end

kindelan =x(1:L)';
k_extrema=x(L+1:end)';
end

function [f,g] = F(L,grerr,X)
% Get coefficients (d) and wavenumbers (k)
d = X(1:L    )';
k = X(L+1:end)';

% Error-bounds
E = (1+(-1).^((1:L)+mod(L,2))*grerr)';

% LHS functions
groupvel =  d .* 2.*([1:L+1/2]-1/2)     * cos( ([1:L+1/2]'-1/2)*k ); % Objective function 1
groupveld= -d .* 2.*([1:L+1/2]-1/2).^2  * sin( ([1:L+1/2]'-1/2)*k ); % Objective function 2

% JACOBIAN FOR F1
groupveltod = (0*d+1) .* 2.*([1:L+1/2]-1/2)     .* cos( k'*([1:L+1/2]-1/2) ) ; % Grad_D
groupveltok =  sin( k'*([1:L+1/2]-1/2) ) *    (-d .* 2.*([1:L+1/2]-1/2).^2)' ; % Grad_K

% JACOBIAN FOR F2
groupveldtod = -(0*d+1) .* 2.*([1:L+1/2]-1/2).^2     .* sin( k'*([1:L+1/2]-1/2) ); % Grad_D
groupveldtok = cos( k'*([1:L+1/2]-1/2) ) *    (-d .* 2.*([1:L+1/2]-1/2).^3)'    ;  % Grad_K

% Build homogeneous objective function of type F(X)=E -> F(X)-E=0.
f = [groupvel' - E;
     groupveld'];

if nargout == 2
    % --- GENERATE JACOBIAN OF TYPE J(i,k)=d(F_i)/d(x_k), i.e.:
    %   [dF_1/da_1 dF_1/da_2 dF_1/dk_1 dF_1/dk_2 ...;
    %    dF_2/da_1 dF_2/da_2 dF_2/dk_1 dF_2/dk_2 ...;
    %                        ... ]
    g = [groupveltod     diag(groupveltok) ;
         groupveldtod    diag(groupveldtok)];
end

end
