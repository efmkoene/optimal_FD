function [mittet] = FD_mittet( L, grerr )
% FD_MITTET uses least-squares optimization cf. Mittet (2000)
%   FD_mittet( L, grerr ) finds finite-difference coefficients that
%   approximate first-order derivatives using a least-squares minimization
%   of the group velocity error, within an error bound.
%
%   This function is based on Mittet (2000): http://www.ipt.ntnu.no/~barn/Myarticles/Mittet2000b.pdf
%
%   For staggered derivatives of type 'f(x+1/2)-f(x-1/2)'.
%
%   Note: doesn't work from ~ L=25 onwards...
%
% Written by Erik Koene, erik.koene@erdw.ethz.ch.
% ETH Zurich, 07-May-2018. 
% Copyright goes to the respective authors of the papers and their publishers. 
warning('off','optimlib:lsqncommon:SwitchToLineSearch');

% Wavenumber error function (to be minimized)
Dxj = @(d,k) (sum( d .* 2.*([1:L]-1/2)    * cos( ([1:L]'-1/2)*k ),1) -1);

% =========================================================================
% Optimize MITTET least squares for increasing wavenumbers
% =========================================================================
% --- Set cut-off point
lim = 0;

% --- Set first (k=0) bound
problem.options = optimoptions(@lsqnonlin,'OptimalityTolerance',eps^3,'FunctionTolerance',eps^3,'display','off');
problem.solver = 'lsqnonlin';
problem.x0 = FD_taylor(L);

% --- COARSE SEARCHING FUNCTION
dk = 0.1;
for kmax = dk:dk:pi
    % --- Create wavenumber range
    k=dk:dk:kmax;
    
    % --- MATLAB optimization for largest k (kmax)
    problem.objective = @(d) Dxj(d,k).^4;    
    x = lsqnonlin(problem);
    
    % --- Perform check if we must leave the loop (i.e., any error > grerr)
    if any(abs(Dxj(x,k)) - grerr >= lim)
        break                             % Takes away with it the kmax that caused the break!!!
    end    
end

% --- FINE SEARCHING FUNCTION (successive halvings)
for j=1:9
dk = dk*0.5;
for kmax = kmax-2*dk:dk:kmax
    % --- Create wavenumber range
    k=dk:dk:kmax;
    
    % --- MATLAB optimization for largest k (kmax)
    problem.objective = @(d) Dxj(d,k).^4;    
    x = lsqnonlin(problem);
    
    % --- Perform check if we must leave the loop (i.e., any error > grerr)
    if any(abs(Dxj(x,k)) - grerr >= lim)
        break                             % Takes away with it the kmax that caused the break!!!
    end    
end
end


mittet = x;
end