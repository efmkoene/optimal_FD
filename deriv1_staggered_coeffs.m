function [taylor,holberg,kindelan,mittet,liu] = deriv1_staggered_coeffs(L,grerr)
% DERIV1_STAGGERED_COEFFS  Compute staggered-grid finite-difference coefficients
%   [taylor,holberg,mittet,liu] = DERIV1_STAGGERED_COEFFS(L,grerr) returns
%   finite-difference coefficients that approximate the first-order
%   derivative using 2*L points.
%    taylor  = Taylor-series coefficients, no error for low wavenumbers but
%              covers small wavenumber range.
%    holberg = Coefficients that trade errors over the entire wavenumber
%              range (including k=0) for coverage of a larger range.
%    kindelan= Alternative (faster) way to derive Holberg's coefficients 
%              using a constrained system of equations.
%    mittet  = Least-squares approximation to Holberg's coefficients. These
%              coefficients are different from Holberg.
%    liu     = Mixes the properties of the two approaches, covering a
%              large wavenumber range.
%
%   A figure is additionally produced to show the group-velocity error that
%   we aim to minimize, and a table provides information about the number
%   of grid-points required per minimum wavelength 
%
%   A staggered finite-difference is of type :
%             f(x+1/2dx) - f(x-1/2dx)
%    f'(x)=~= -----------------------,
%                        dx
%   or generally for higher-order schemes with 'half-order' operator length L:
%               L      f(x+i/2 dx) - f(x-i/2 dx)
%    f'(x)=~= sum  a_i -------------------------, with a_i the finite-difference coefficients.
%              i=1               i * dx
%
%   The stencil thus works for data at equally spaced points ONLY!
%
% -------------------------------------------------------------------------
%   # EXAMPLES
%   ## Example run 1:
%           taylor = deriv1_coeffs_holberg(2,0)
%   ## Returns: taylor = [ 1.1250   -0.0417 ].
%     which is the general finite-difference scheme [9/8, -1/24] (Fornberg, 1988; Table 2)
%
%   ## Example run 2:
%           [taylor,holberg,liu] = deriv_coeffs_holberg(2,0.00235)
%   ## Returns  taylor = [ 1.1250   -0.0417 ],
%              holberg = [ 1.1449   -0.0491 ],
%             kindelan = [ 1.1449   -0.0491 ],
%               mittet = [ 1.1422   -0.0480 ],
%                  liu = [ 1.1405   -0.0468 ].
%     which contains additionally the 'optimized' finite-difference 
%     coefficients that allow a wavenumber-error in the group-velocity of 
%     maximally 0.235%.
%
% -------------------------------------------------------------------------
%   # BACKGROUND INFORMATION
%   ## Taylor-series information: https://amath.colorado.edu/faculty/fornberg/Docs/MathComp_88_FD_formulas.pdf
%   I derive the Taylor-series coefficients straight from polynomial
%   interpolation coefficients, but they can be derived in many ways. Their
%   'problem' is that they have perfect accuracy at low wavenumbers but
%   poor accuracy at high wavenumbers.
%   The optimized methods cover a larger wavenumber-range at the cost of
%   introducing errors in the lower wavenumbers.
%   ## Holberg (1987): https://doi.org/10.1111/j.1365-2478.1987.tb00841.x
%   The optimization function is Equation 12. Note that an error here is
%   allowed even at k=0! This is why it, in the end, can cover the largest
%   wavenumber range. It is based on minimizing the global group velocity
%   error.
%   ## Kindelan, Kamel and Sguazzero (1990): https://library.seg.org/doi/pdf/10.1190/1.1442763
%   The scheme here recognizes that there are L extrema in the error plot
%   of Holberg at 1±grerr, creating a constraint system of equations that
%   finds the coefficients and locations of extrema. This method is much
%   faster than Holberg (1987) at finding the coefficients, though required
%   some work to truly constrain the system on my side.
%   ## Mittet (2000): http://www.ipt.ntnu.no/~barn/Myarticles/Mittet2000b.pdf
%   This scheme uses a least-squares procedure to obtain the coefficients,
%   more in line with Holberg's intentions perhaps, and much faster to
%   approximate. However, the coefficients are less perfect in the sense of
%   filling the wavenumber-error space. See e.g., Table 1 (p. F37) of this
%   publication https://library.seg.org/doi/pdf/10.1190/1.3278525 where an
%   error of 0.3% still only gives d1=1.00235, whereas Holberg/Kindelan
%   would give a larger error d1=1.003.
%   ## Liu (2014): https://doi.org/10.1093/gji/ggu032
%   The optimization function is Equation 12. Note that the error here is
%   fixed at 0 for k=0, and only grows larger afterwards. It is based on
%   minimizing the dispersion relation for a wavenumber range. Using a
%   simple search-function I try to minimize the Liu coefficients within
%   the group-velocity error function as well.
%   There are many Least-Squares schemes to find finite-difference
%   coefficients, the one by Liu (2014) just happens to work fast.
%
% Written by Erik Koene, erik.koene@erdw.ethz.ch.
% ETH Zurich, 15-May-2018. 
% Copyright goes to the respective authors of the papers and their publishers. 
% 

% =========================================================================
% TAYLOR COEFFICIENTS
% =========================================================================
taylor = FD_taylor(L);

% =========================================================================
% HOLBERG COEFFICIENTS
% =========================================================================
[holberg,k] = FD_holberg(L,grerr);
disp('Optimization Holberg done')

% =========================================================================
% KINDELAN COEFFICIENTS
% =========================================================================
[kindelan,k_extrema] = FD_kindelan(L,grerr);
disp('Optimization Kindelan done')

% =========================================================================
% MITTET COEFFICIENTS
% =========================================================================
[mittet] = FD_mittet(L,grerr);
disp('Optimization Mittet done')

% =========================================================================
% LIMIT LIU COEFFICIENTS WITHIN GROUP VELOCITY ERROR RANGE
% =========================================================================
liu = FD_liu(L,1,grerr);
disp(['Optimization Liu done'])

% =========================================================================
% Plot wavenumber performance
% =========================================================================
groupvelerror = @(d,k) (sum( d .* 2.*([1:L+1/2]-1/2)    * cos( ([1:L+1/2]'-1/2)*k ),1) -1); % Error function
figure(2)
% subplot(2,1,1)
plot( k, groupvelerror(taylor,   k), ...
      k, groupvelerror(holberg,  k), ...
      k, groupvelerror(kindelan, k), ...
      k, groupvelerror(mittet,   k), ...
      k, groupvelerror(liu,      k),'k.' )
% Error bounds
hold on
plot( [xlim],[grerr   grerr],'k--')
plot( [xlim],[-grerr -grerr],'k--')
hold off
legend('Taylor', 'Holberg (1987)','Kindelan (1990)','Mittet (2000)','Liu (2014)')

grid on
grid minor
title('Wavenumber response')
xlabel('k\Deltax')
ylabel('Error in group velocity')

% DERIVATIVE
% derivativeofthat = @(d,k) (sum( -d .* 2.*([1:L+1/2]-1/2).^2    * sin( ([1:L+1/2]'-1/2)*k ),1)); % Error function
% subplot(2,1,2)
% plot( k, derivativeofthat(taylor,   k), ...
%       k, derivativeofthat(holberg,  k), ...
%       k, derivativeofthat(kindelan, k), ...
%       k, derivativeofthat(mittet,   k), ...
%       k, derivativeofthat(liu,      k),'k.' )
% % Error bounds
% hold on
% plot( [xlim],[grerr   grerr],'k--')
% plot( [xlim],[-grerr -grerr],'k--')
% hold off
% legend('Taylor', 'Holberg (1987)','Kindelan (1990)','Mittet (2000)','Liu (2014)')
% 
% grid on
% grid minor
% title('Wavenumber response')
% xlabel('k\Deltax')
% ylabel('Error in group velocity')



% =========================================================================
% Compute some data
% =========================================================================
FDCoeffType = {'Taylor';'Holberg (1987)';'Kindelan (1990)';'Mittet (2000)';'Liu (2014)'};
CriticalWavenumber = [k(find(groupvelerror(taylor,   k)<-grerr*1.0001,1)); ...
                      k(find(groupvelerror(holberg,  k)<-grerr*1.0001,1)); ...
                      k(find(groupvelerror(kindelan, k)<-grerr*1.0001,1)); ...
                      k(find(groupvelerror(mittet,   k)<-grerr*1.0001,1)); ...
                      k(find(groupvelerror(liu,      k)<-grerr*1.0001,1));];
GridPointsPerWavelengthForAccuracy = pi*2 ./ CriticalWavenumber;
table1 = table(FDCoeffType,CriticalWavenumber,GridPointsPerWavelengthForAccuracy)

end