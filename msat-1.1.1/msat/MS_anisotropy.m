% MS_ANISOTROPY - Simple measures of anisotropy 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate the degree of anisotropy of an elasticity matrix. 
%
% [ uA, ... ] = MS_anisotropy( C, ... )
%
% Usage: 
%     [ uA ] = MS_anisotropy( C )                    
%         Return the Universal Elastic Anisotropy Index of Ranganathan
%         and Ostoja-Starzewski (2008). Valid for any elasticity matrix,
%         uA is zero for an isotropic case and increases for increasing 
%         anisotropy.
%
%     [ uA, lmA ] = MS_anisotropy( C )
%         Also return the general measure of anisotropy proposed by
%         Ledbetter and Miglion (2006). This the the ratio of the fastest
%         and slowest squared shear wave velocity over all propogation
%         and polarization directions. Equal to one in the isotropic case,
%         increases with increasing anisotropy.
%
%     [ uA, lmA, zA ] = MS_anisotropy( C )
%         Also return the Zenner (1948) measure of anisotropy. This is
%         only valid for cubic crystals (NaN is returned if C does not 
%         represent a cubic crystal). zA is 1 for an isotropic case and
%         increases or decreases with incresing anisotropy.
%
%     [ uA, lmA, zA, cbA ] = MS_anisotropy( C )
%         Also return the Chung-Buessem (1967) anisotropy index. This
%         is a single valued measure of anisotropy derived from zA. Like
%         uA, this is zero for an isotropic case and increases for increasing 
%         anisotropy. Only valid for matricies representing cubic crystals.
%
%     [ uA, lmA, ... ] = MS_anisotropy( C, n )
%         Set the number of random directions to sample for the calculation
%         of lmA. Defaults to 1000, which seems to give results accurate to
%         two decimal places. Ledbetter and Miglion (2006) use 10000 which 
%         gives results reproducable to three decimal places and a 
%         noticable slow down.
%
% Notes:
%     These measures of anisotropy are independent of orientation. However,
%     the test for cubic symmetry assumes the matrix is in an ideal
%     orention. Use MS_AXES to reorentate the imput matrix for the general 
%     case. MS_NORMS can be used to provide an alternate measure of
%     anisotropy. Ledbetter and Miglion (2006) claim lmA is identcal to zA
%     for cubic cases but Ranganathan and Ostoja-Starzewski (2008) point 
%     out cases where zA < 1 while lmA > 1 by construction. 
%
% References:
%     Zenner, C. (1948) Elasticity and Anelasticiy of Metals. University 
%     of Chicago. 
%
%     Chung, D. H. and W. R. Buessem (1967) Journal of Applied Physics 
%     vol.38 p.5 
%
%     Ledbetter, H. and A. Miglion (2006) "A general elastic-anisotropy 
%     measure" Journal of Applied Physics vol.100 art.num.063516
%     http://dx.doi.org/10.1063/1.2338835
%
%     Ranganathan, S. I. and M. Ostoja-Starzewski (2008) "Universal Elastic
%     Anisotropy Index" Physical Review Letters vol.101 art.num.055504.
%     http://dx.doi.org/10.1103/PhysRevLett.101.055504
%
% See also: MS_POLYAVERAGE, MS_NORMS, MS_AXES, MS_PHASEVELS

% Copyright (c) 2011, James Wookey and Andrew Walker
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [ uA, lmA, zA, cbA] = MS_anisotropy( C, varargin )

    numprop = 1000;
    % process the optional arguments
    if nargin==2
        numprop = varargin{1};
    end
    
    uA = UniversalAnisotropyIndex(C);
    if (nargout == 1), return, end % Skip lmA if we don't need it
    lmA = LedbetterAnisotropyIndex(C, numprop);
    if (isCubic( C ))
        zA = ZenerAnisotropyIndex(C);
        cbA = (3*(zA-1)^2) / (3*(zA-1)^2 + 25*zA);
    else
        zA = NaN;
        cbA = NaN;
    end
end

function i = isCubic ( C )

    i = (C(1,1) == C(2,2)) ...
      & (C(1,1) == C(3,3)) ...
      & (C(1,2) == C(1,3)) ...
      & (C(1,2) == C(2,3)) ...
      & (C(4,4) == C(5,5)) ...
      & (C(4,4) == C(6,6)) ...
      & (C(1,4) == 0.0) ...
      & (C(1,5) == 0.0) ...
      & (C(1,6) == 0.0) ...
      & (C(2,4) == 0.0) ...
      & (C(2,5) == 0.0) ...
      & (C(2,6) == 0.0) ...
      & (C(3,4) == 0.0) ...
      & (C(3,5) == 0.0) ...
      & (C(3,6) == 0.0) ...
      & (C(4,5) == 0.0) ...
      & (C(4,6) == 0.0) ...
      & (C(5,6) == 0.0);
end

function zA = ZenerAnisotropyIndex( C )

    zA = C(4,4)*2.0 / (C(1,1) - C(1,2)); 

end 

function uA = UniversalAnisotropyIndex( C )

    [~, ~, B_v, G_v, B_r, G_r] = MS_polyaverage( C );
    uA = (5*(G_v/G_r) + (B_v/B_r)) - 6;
    
end

function lmA = LedbetterAnisotropyIndex( C, n )

    points = randomSphere(n);
    % Density should not matter - we ratio velocities.
    [~,~,vs1,vs2,~] = MS_phasevels(C,3000,points(:,2),points(:,1));
    lmA = (max(vs1))^2 / (min(vs2))^2; 
    
end 

function points = randomSphere( n )
    % n randomly distributed points on a sphere (polar co-ordinates, 
    % in degrees). Points(:,1) are azimuths, points(:,2) are inclinations.
    points = rand(n,2);
    points(:,1) = points(:,1)*pi*2.0;
    points(:,2) = acos(2.0*points(:,2)-1);
    points = points*(180.0/pi);
    points(:,2) = points(:,2) - 90.0;
end 