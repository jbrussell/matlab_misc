% MS_BUILD_ISOTROPIC - Create elasticity matrix from pairs of isotropic moduli.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Given two isotropic moduli create an elasticity matrix for an isotropic 
%     material and optionally return all other moduli.
%
%  %  [C, ...]=MS_build_isotropic('modname1', mod1, 'modname2', mod2')
%
% Usage: 
%     [C]=MS_build_isotropic('modname1', mod1, 'modname2', mod2')                    
%          Given two isotropic moduli mod1 and mod2 named as strings 
%          'modname1' and 'modname2' return an elasticity matrix C.
%          Permitted moduli are:
%              'lam' - first lame parameter, in GPa
%              'mu' - shear modulus (second Lame parameter), in GPa
%              'K' - bulk modulus, in GPa.
%              'E' - Young's modulus, in GPa.
%              'nu' - Poisson's ratio, dimensionless.
%              'M' - 'P-wave modulus', in GPa.
%
%     [C, K, E, lam, mu, nu, M]=MS_build_isotropic('modname1', mod1, 'modname2', mod2') 
%          Also return all 6 moduli.
%
% Notes:
%     Not all combinations of moduli are valid / implemented
%
% References:
%     Stacey F. D. and Davis, P. M. (2008) Physics of the Earth 4th ed. 
%     Appendix D and Section 10.2.
%
% See also: MS_POLYAVERAGE

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

function [ C, K, E, lam, mu, nu, M ] = MS_build_isotropic( varargin )

    % Get Lam\'e constants from input arguments.
    in_args = char(lower(strcat(varargin(1) , varargin(3))));
    switch in_args
        case 'lammu'
            lam = varargin{2};
            mu = varargin{4};
        case 'mulam'
            lam = varargin{4};
            mu = varargin{2};
        case 'emu'
            mu = varargin{4};
            lam = (mu*(varargin{2}-2*mu))/(3*mu-varargin{2});
        case 'mue'
            mu = varargin{2};
            lam = (mu*(varargin{4}-2*mu))/(3*mu-varargin{4});
        case 'klam'
            lam = varargin{4};
            mu = 3*(varargin{2}-lam)/2;
        case 'lamk'
            lam = varargin{2};
            mu = 3*(varargin{4}-lam)/2;
        case 'kmu'
            lam = varargin{2}-(2/3)*varargin{4};
            mu = varargin{4};
        case 'muk'
            lam = varargin{4}-(2/3)*varargin{2};
            mu = varargin{2};
        case 'lamnu'
            lam = varargin{2};
            mu = (lam*(1-2*varargin{4}))/(2*varargin{4});
        case 'nulam'
            lam = varargin{4};
            mu = (lam*(1-2*varargin{2}))/(2*varargin{2});
        case 'munu'
            mu = varargin{2};
            lam = 2*mu*varargin{4}/(1-2*varargin{4});
        case 'numu'
            mu = varargin{4};
            lam = 2*mu*varargin{2}/(1-2*varargin{2});
        case 'enu'
            lam = (varargin{2}*varargin{4})/...
                ((1+varargin{4})*(1-2*varargin{4}));
            mu = varargin{2}/(2*(1+varargin{4}));
        case 'nue'
            lam = (varargin{4}*varargin{2})/...
                ((1+varargin{2})*(1-2*varargin{2}));
            mu = varargin{4}/(2*(1+varargin{2}));
        case 'knu'
            lam = 3*varargin{2}*varargin{4}/(1+varargin{4});
            mu = 3*varargin{2}*(1-2*varargin{4})/(2*(1+varargin{4}));
        case 'nuk'
            lam = 3*varargin{4}*varargin{2}/(1+varargin{2});
            mu = 3*varargin{4}*(1-2*varargin{2})/(2*(1+varargin{2}));
        case 'ke'
            lam = (3*varargin{2}*(3*varargin{2}-varargin{4}))/...
                (9*varargin{2}-varargin{4});
            mu = (3*varargin{2}*varargin{4})/(9*varargin{2}-varargin{4});
        case 'ek'
            lam = (3*varargin{4}*(3*varargin{4}-varargin{2}))/...
                (9*varargin{4}-varargin{2});
            mu = (3*varargin{4}*varargin{2})/(9*varargin{4}-varargin{2});
        case 'mmu'
            mu = varargin{4};
            lam = varargin{2} - 2*mu;
        case 'mum'
            mu = varargin{2};
            lam = varargin{4} - 2*mu;
        otherwise
            error('MS:buildiso:wrongargs', ...
                'Cannot construct C from these arguments')
    end
   

    M = 2*mu+lam;
    E = (mu*((3*lam)+(2*mu)))/(lam+mu);
    K = lam + (2/3)*mu;
    nu = lam/(2*(lam+mu));
    
    % build C
    C = [M   lam lam 0.0 0.0 0.0; ...
         lam M   lam 0.0 0.0 0.0; ...
         lam lam M   0.0 0.0 0.0; ...
         0.0 0.0 0.0 mu  0.0 0.0; ...
         0.0 0.0 0.0 0.0 mu  0.0; ...
         0.0 0.0 0.0 0.0 0.0 mu ];

     
end