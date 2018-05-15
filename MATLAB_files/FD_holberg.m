function [holberg,k] = FD_holberg( L, grerr )
% FD_HOLBERG  Compute the optimized finite-difference coefficients cf. Holberg (1987)
%   [holberg] = FD_holberg(L,grerr,taylor) returns the finite-difference
%   coefficients that 'optimally' approximate the first-order derivative
%   for an operator of half-order L, with a maximum error in group velocity
%   of 'grerr'.
%
%   [holberg,k] = FD_holberg(L,grerr,taylor) returns the wavenumber vector
%   'k' from - until k_critical.
%
%   For staggered derivatives of type 'f(x+1/2)-f(x-1/2)'.
%
% Written by Erik Koene, erik.koene@erdw.ethz.ch.
% ETH Zurich, 07-May-2018. 
% Copyright goes to the respective authors of the papers and their publishers. 


% --- We look for compact operators, which for half-offset is -> (f(x+1/2dx) - f(x-1/2dx))/dx even number of elements desired
L=L-1/2;

% Equation (12, Holberg) integrand: 
Dxj = @(d,k) real(sum( d' .* ((-L:L)') .* exp(1i*((-L:L)')*k) ,1 ))-1;

% Equation (12, Holberg) integrated (through simple summation)
Dxj_integrated =@(d,k,dk) ( sum(Dxj(d,k).^6)  );


% =========================================================================
% Optimize FORNBERG for increasing wavenumbers
% =========================================================================
% --- Set first (k=0) bound
problem.options = optimoptions(@fmincon,'Algorithm','sqp','Display','off','FinDiffType','central');
problem.solver = 'fmincon';
taylor = FD_taylor(L+1/2);
problem.x0 = [-taylor(end:-1:1) taylor];
problem.Aeq = ones(1,2*L+1); % / all coefficients must sum to 0
problem.beq = 0;             % \

% --- Set cut-off point
lim = eps*10;

% --- COARSE SEARCHING KMAX
dk = 0.1;
for kmax = 0:dk:pi
    % --- Create wavenumber range
    k=0:dk:kmax;
    
    % --- MATLAB optimization for largest k (kmax)
    problem.objective = @(d) Dxj_integrated(d,k,dk);    
    [problem.Aineq,problem.bineq] = create_mats( L, k, grerr );    
    
    % --- Find minimum...global cannot be guaranteed...
    [x,~] = fmincon(problem);
    
    % --- Perform check if we must leave the loop (i.e., any error > grerr)
    c = abs(((-L:L)) .* exp( k'*1i*((-L:L)) ) * x' - 1) - grerr;
    % >>> if any c > 0 -> terminate the loop
    if any(c(1:length(k)) >= lim)
        break
    end    
end

% --- FINE SEARCHING OF KMAX using 5 consecutive halvings of dk (identical code block as above) 
for j=1:5
dk = dk*0.5;
for kmax = kmax-2*dk:dk:kmax
    % --- Create wavenumber range
    k=0:dk:kmax;
    
    % --- MATLAB optimization for largest k (kmax)
    problem.objective = @(d) Dxj_integrated(d,k,dk);    
    [problem.Aineq,problem.bineq] = create_mats( L, k, grerr );    
    
    % --- Find minimum...global cannot be guaranteed...
    [x,~] = fmincon(problem);
    
    % --- Perform check if we must leave the loop (i.e., any error > grerr)
    c = abs(((-L:L)) .* exp( k'*1i*((-L:L)) ) * x' - 1) - grerr;
    % >>> if any c > 0 -> terminate the loop
    if any(c(1:length(k)) >= lim)
        break
    end    
end
end

% Do the final optimization given the found parameters
% --- Create wavenumber range
k = k(1:end-1);
% --- MATLAB optimization for largest k (kmax)
problem.objective = @(d) Dxj_integrated(d,k,dk);    
[problem.Aineq,problem.bineq] = create_mats( L, k, grerr );    
% --- Find minimum...global cannot be guaranteed...
[x] = fmincon(problem);

k=[0:dk:k(end)+50*dk];
holberg = x(end/2+1:end);

end

function [A,b] = create_mats( L, k, grerr )
    % --- Create matrix A. L=row vector [-L,-L+1,...,L], k=column vector [0;dk;...;kmax]
    Af = ((-L:L)) .* exp( k'*1i*((-L:L)) );
    A = [real(Af);
         imag(Af);
        -real(Af)];
    % Set upper and lower bounds
    b = [ 1+grerr * ones(2*length(k),1);
         -1+grerr * ones(  length(k),1)];
end
