function [c] = FD_taylor( L )
% FD_TAYLOR computes the 'Taylor' finite-difference coefficients
%   [c] = FD_taylor(L) returns the coefficients for Taylor expansion, equal
%   to dividing Lagrange coefficients by L.
%
%   For staggered derivatives of type 'f(x+1/2)-f(x-1/2)'.
%
%   Written by Erik Koene, erik.koene@erdw.ethz.ch.
%   ETH Zurich, 07-May-2018.
%   Using snippet from:
%   https://ccrma.stanford.edu/~jos/Interpolation/Matlab_Code_Lagrange_Fractional.html

    L = 2*L-1;
    delay=L/2;
    n = 0:L;
    h = ones(1,L+1);
    for k = 0:L
        index = find(n ~= k);
        h(index) = h(index) *  (delay-k)./ (n(index)-k);
    end
    
    c = h(ceil(delay)+1:end) * 2 ./ (1:2:L);