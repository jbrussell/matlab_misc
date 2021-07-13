% MS_TI - generate elastic constants for a vertically transverse isotropic
%          medium from specified parameter sets. Symmetry is in the 3-axis. 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
%  [C]=MS_TI( list_of_parameters , parameter_set_string )
%
%  where parameter_set_string defines the set which precede it: 
%
%'thomsen' (default)
%~~~~~~~~~~~~~~~~~~~
%
%   [C]=MS_TI(vp,vs,rh,eps,gam,del) -or-
%   [C]=MS_TI(vp,vs,rh,eps,gam,del,'thomsen')
%
%   Inputs: 
%       rh  : Density (kg/m2) 
%       vp  : km/s (vertical) 
%       vs  : km/s (vertical)
%       eps, gam, del : Dimensionless
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Given Thomsen (1986) parameters for a weakly anisotropic VTI medium, 
%   return the elasticity matrix
%
%'panning'
%~~~~~~~~~
%
%   [C]=MS_TI(vp,vs,rh,xi,phi,eta,'panning')
%
%   Inputs: 
%       rh  : Density (kg/m2) 
%       vp  : km/s (isotropic average)
%       vs  : km/s (isotropic average)
%       xi, phi, eta : Dimensionless anisotropy parameters
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Calculates the elastic tensor for a VTI medium from average Vp and Vs,
%   and anisotropic parameters xi, phi and eta (see, e.g., Panning and 
%   Romanowicz, 2006). Derivation only valid if eta = 1 and phi = 1.
%
%'global'
%~~~~~~~~~
%
%   [C]=MS_TI(vp,vs,rh,xi,phi,eta,'global')
%
%   Inputs: 
%       rh  : Density (kg/m2) 
%       vp  : km/s (Voigt isotropic average)
%       vs  : km/s (Voigt isotropic average)
%       xi, phi, eta : Dimensionless anisotropy parameters
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Calculates the elastic tensor for a VTI medium from the Voigt average
%   Vp and Vs and anisotropic parameters xi, phi and eta (see, e.g., Babuska
%   and Cara, 1991)
%
%'love'
%~~~~~~
%
%   [C]=MS_TI(A,C,L,N,F,'love')
%
%   Inputs: 
%        A,C,L,N,F : Love parameters (GPa)
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Calculates the elastic tensor for a VTI medium from Love (1927) parameters.
%
%'anderson'
%~~~~~~~~~~
%
%   [C]=MS_TI(Vpv, Vph, Vp45, Vsv, Vsh, rho, 'anderson')
%
%   Inputs:
%        Vpv    : Velocity of P-wave along 3-axis
%        Vph    : Velocity of P-wave normal to 3-axis
%        Vp45   : Velocity of P-wave at 45 degrees to 3-axis
%        Vsv    : Velocity of S-wave along 3-axis (or normal to 3-axis
%                 polarised normal to 1-2 plane)
%        Vsh    : Velocity of S-wave normal 3-axis polarised
%                 parallel to 1-2 plane
%        rho    : Density
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Calculates the elastic tensor for a VTI medium from Anderson's (1961)
%   parameters. Velocities in km/s, density in kg/m^3.
%
%References
%~~~~~~~~~~
%
%     Anderson, D. L. (1961) "Elastic wave propagation in layered
%         anisotropic media" Journal of Geophysical Research 66:2953 - 2963
%
%     Thomsen, L. (1986) "Weak elastic anisotropy" Geophysics 
%         vol.51 pp.1954-1966
%
%     Mark Panning and Barbara Romanowicz (2006) A three-dimensional radially 
%        anisotropic model of shear velocity in the whole mantle. Geophysical  
%        Journal International v167, 361–379. 
%        doi: 10.1111/j.1365-246X.2006.03100.x
%
%     Babuska, V. and Cara, M. (1991). Seismic Anisotropy in the Earth. Kluwer 
%        Academic, Boston.
%
%     Love, A.E.H., (1927). A Treatise on the Theory of Elasticity, 
%        Cambridge Univ. Press, Cambridge.
%
% See also: MS_iso, MS_elasticDB

% Copyright (c) 2011-2013, James Wookey and Andrew Walker
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

function [C]=MS_TI(varargin)
   
% see if a parameter set type was selected. 
if ischar(varargin{nargin})
   pset = varargin{nargin} ;
   ncheck = nargin-1;
else
   pset = 'thomsen' ;
   ncheck = nargin ;
end   

% check that all other inputs are scalars
for icheck=1:ncheck
   if ~isscalar(varargin{icheck})
      error('MS:TI:BadScalarInput',...
      'MS_TI requires scalar inputs') ;
   end   
end

%  call the appropriate routine.

switch lower(pset)
%-------------------------------------------------------------------------------
   case {'thomsen', 'thom'}
      if length(varargin)~=7 & length(varargin)~=6 % need to check 
         error('MS:TI:ThomsenWrongArgs', ...
         'Thomsen (1986) requires 6 input parameters.') ;
      end
      vp = varargin{1} ; vs = varargin{2} ; rh = varargin{3} ; 
      eps = varargin{4} ; gam = varargin{5} ; del = varargin{6} ;
      [C]=MS_thomsen(vp,vs,rh,eps,gam,del) ;
%-------------------------------------------------------------------------------
   case {'panning'}
      if length(varargin)~=7 % need to check 
         error('MS:TI:PanningWrongArgs', ...
         'Panning (2006) requires 6 input parameters.') ;
      end
      vp = varargin{1} ; vs = varargin{2} ; rh = varargin{3} ; 
      xi = varargin{4} ; phi = varargin{5} ; eta = varargin{6} ;
      [C]=MS_panning(vp,vs,rh,xi,phi,eta) ;
%-------------------------------------------------------------------------------
   case {'global'}
      if length(varargin)~=7 % need to check
         error('MS:TI:GlobalWrongArgs', ...
         'global requires 6 input parameters.') ;
      end
      vp = varargin{1} ; vs = varargin{2} ; rh = varargin{3} ; 
      xi = varargin{4} ; phi = varargin{5} ; eta = varargin{6} ;
      [C]=MS_global(vp,vs,rh,xi,phi,eta) ;
%-------------------------------------------------------------------------------
   case {'love'}
      if length(varargin)~=6 % need to check 
         error('MS:TI:LoveWrongArgs', ...
         'Love (1926) requires 5 input parameters.') ;
      end
      A = varargin{1} ; C_love = varargin{2} ; L = varargin{3} ; 
      N = varargin{4} ; F = varargin{5} ; 
      [C]=MS_love(A,C_love,L,N,F) ;
%-------------------------------------------------------------------------------
   case {'anderson'}
      if length(varargin)~=7 % need to check 
         error('MS:TI:AndersonWrongArgs', ...
         'Anderson (1961) requires 6 input parameters.') ;
      end
      Vpv = varargin{1} ; Vph = varargin{2}; Vp45 = varargin{3};
      Vsv = varargin{4} ; Vsh = varargin{5}; rho = varargin{6};
      [C]=MS_anderson(Vpv, Vph, Vp45, Vsv, Vsh, rho);
%--------------------------------------------------------------------------
   otherwise
      error('MS:TI:UnknownParSet', ...
         'Specified parameter set is not supported.') ;
%-------------------------------------------------------------------------------
end % of switch
   
% check resulting matrix.
MS_checkC(C) ;
   
end

function [CC]=MS_love(A,C,L,N,F)
   CC= [ A      A-2.*N F  0  0  0  ; ...
         A-2.*N A      F  0  0  0  ; ...
         F      F      C  0  0  0  ; ...
         0      0      0  L  0  0  ; ...
         0      0      0  0  L  0  ; ...
         0      0      0  0  0  N ] ;
end

function [C]=MS_panning(vp,vs,rh,xi,phi,eta)

%  convert to m/s
   vp=vp*1e3;
   vs=vs*1e3;

   if eta ~= 1.0
       warning('MS:TIeta', ['Warning, the derivation of the ' ...
           '"panning" mode in MS_TI assumes eta = 1. Use "global" for ' ...
           'correct general treatment']);
   end
   
   if phi ~= 1.0
       warning('MS:TIphi', ['Warning, the derivation of the ' ...
           '"panning" mode in MS_TI assumes phi = 1. Use "global" for ' ...
           'correct general treatment']);
   end
   
%  note: xi and phi are defined oppositely; i.e.: 
%  xi = vsv^2/vsh^2 and phi = vph^2/vpv^2
   vsv = sqrt((3.*vs.^2)./(2+xi)) ;
   vsh = sqrt(xi.*vsv.^2) ;
   
   vph = sqrt((5.*vp.^2)./(4+phi)) ;
   vpv = sqrt(phi.*vph.^2) ;
   
   C11 = vph.^2.*rh ; % A
   C33 = vpv.^2.*rh ; % C
   C44 = vsv.^2.*rh ; % L
   C66 = vsh.^2.*rh ; % N
   
   C12 = C11-2.*C66 ;
   C13 = eta.*(C11-2.*C44) ; % F
   
   C22 = C11 ;
   C55 = C44 ;
   
   C = [C11 C12 C13  0   0   0  ; ...
        C12 C22 C13  0   0   0  ; ...
        C13 C13 C33  0   0   0  ; ...
         0   0   0  C44  0   0  ; ...
         0   0   0   0  C55  0  ; ...
         0   0   0   0   0  C66 ] ;
   
   C = C./1e9 ; % convert to GPa

end

function [Cglobal]=MS_global(vp,vs,rho,xi,phi,eta)
%  Output the elastic tensor given a set of radial anisotropy parameters
%  as used typically in global seismology.  Average velocities are given by:
%        15*rho*<Vp>^2 = 3*C + (8 + 4*eta)*A + 8*(1 - eta)*L
%        15*rho*<Vs>^2 =   C + (1 - 2*eta)*A + (6 + 4*eta)*L + 5*N
%     vp:   Voigt average P wave velocity
%     vs:   Voigt average shear wave velocity
%     rho:  Density
%     xi:   (Vsh^2/Vsv^2) of horizontal waves
%     phi:  (Vpv^2/Vph^2)
%     eta:  C13/(C11 - 2*C44)

% convert to m/s
   vp=vp*1e3;
   vs=vs*1e3;

   L = 15.*rho.*((3.*phi + 8. + 4.*eta).*vs.^2 - (phi + 1. - 2.*eta).*vp.^2) ...
         ./ ((6. + 4.*eta + 5.*xi).*(3.*phi + 8. + 4.*eta) ...
            - 8.*(phi + 1. - 2.*eta).*(1. - eta)) ;

   A = (15.*rho.*vp.^2 - 8.*(1. - eta).*L) ./ (3.*phi + 8. + 4.*eta) ;

   F = eta.*(A - 2.*L) ;
   C = phi.*A ;
   N = xi.*L ;
   C12 = A - 2.*N ;

   Cglobal = [ A  C12 F 0 0 0; ...
			  C12  A  F 0 0 0; ...
			   F   F  C 0 0 0; ...
			   0   0  0 L 0 0; ...
			   0   0  0 0 L 0; ...
			   0   0  0 0 0 N] ;

   Cglobal = Cglobal./1e9 ; % convert to GPa

end

function [C]=MS_thomsen(vp,vs,rh,eps,gam,del)

%  convert to m/s
   vp=vp*1e3;
   vs=vs*1e3;

   C=zeros(6,6) ;
   C(3,3) = vp*vp ; % Eq 9a in Thomsen paper.
   C(4,4) = vs*vs ; % 9b
   C(6,6) = C(4,4)*(2.0*gam +1.0) ; % 8b
   C(1,1) = C(3,3)*(2.0*eps +1.0) ; % 8a
   
   btm = 2.0*C(4,4) ;
   term = C(3,3) - C(4,4) ;
   ctm = C(4,4)*C(4,4) - (2.0*del*C(3,3)*term + term*term) ;
   dsrmt = (btm*btm - 4.0*ctm) ;
   
	if dsrmt < 0
		error('MS:TI:ThomsenHiVS',...
		   'S-velocity too high or delta too negative for Thomsen routine.') ;
	end
   
   C(1,3) = -btm/2.0 + sqrt(dsrmt)/2.0 ;
   
   C(1,2) = C(1,1) - 2.0*C(6,6) ; 
   C(2,3) = C(1,3) ;
   C(5,5) = C(4,4) ;
   C(2,2) = C(1,1) ;

   % make symmetrical
   for i=1:6
      for j=i:6
         C(j,i) = C(i,j) ;
      end
   end

%  convert to GPa
   C = C.*rh./1e9 ;

end

function [C]=MS_anderson(Vpv, Vph, Vp45, Vsv, Vsh, rho)

    % Convert V to m/s
    Vpv = Vpv*1e3;
    Vph = Vph*1e3;
    Vp45 = Vp45*1e3;
    Vsv = Vsv*1e3;
    Vsh = Vsh*1e3;
    
    C = zeros(6);
    C(1,1) = Vph^2 * rho;
    C(3,3) = Vpv^2 * rho;
    C(4,4) = Vsv^2 * rho;
    C(1,2) = C(1,1) - (Vsh^2 * 2.0 * rho);
    C(6,6) = 0.5*(C(1,1)-C(1,2));
    C(1,3) = sqrt( ((2.0 * rho * Vp45^2) - ...
                   (0.5*(C(1,1)+C(3,3)+2*C(4,4))))^2 - ...
                   (0.5*(C(1,1)-C(3,3))^2)) - C(4,4);
    % Fill in the gaps.
    C(1,2) = 0.0 ; % So we can expand
    C = MS_expand(C, 'vti');        
    % Convert C from Pa to GPa
    C = C./1e9;

end
