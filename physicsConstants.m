function pc = physicsConstants()
% PHYSICSCONSTANTS Define physics constants (2019) to their known accuracy.
%
%   [pc] = physicsConstants() 
%                returns a pc structure containing the constants.
%                The default units are SI units kg, m, s, C, J, etc. 
%
%   See also symunit.
% Dr. Jodi Christiansen - Cal Poly San Luis Obispo

pc = struct('c',299792458, ... % m/s <-- 15 sig figs.
    'G', 6.67430e-11, ... % m^3/kg/s^2 <-- 4 sig figs. 
    'g', 9.80665, ... % m/s^2 <-- 15 sig figs.
    'e', 1.602176634e-19, ... % C <-- 9 sig figs.
    'h', 6.62607015e-34, ... % J*s <-- 6 sig figs.
    'h_eV', 0, ... % eV*s <-- 6 sig figs.
    'hbar', 1.054571817e-34, ... % J*s <-- exact* 
    'hbar_eV', 6.582119569e-19 , ... % eV*s <-- exact*
    'hc', 0, ... % see below
    'hc_eV', 0, ... 
    'muB', 0, ...
    'epsilonNot', 8.8541878128e-12, ... % F/m = farad/meter = C^2/(N m^2) <-- 9 sig figs
    'muNot', 4*pi*1.00000000055e-7, ... % N/A^2 <-- 10 sig figs
    'K', 0, ... % 1/(4*pi*epsilonNot) <-- ?? sig figs.
    'K_eV', 0, ... % K but in eVmC^-2
    'Na', 6.02214076e23, ... %  mol^-1 (Avogadro's number) <-- exact
    'kb', 1.38064849e-23, ... % J/K. <-- 6 sig figs
    'me', 9.1093837015e-31, ...  % kg <-- 9 sig figs
    'u',  1.66053906660e-27, ... % kg <-- 9 sig figs
    'mp', 0, ... % kg <-- 9 sig figs
    'mn', 0, ... % kg <-- 9 sig figs
    'me_MeV',0 , ... % kg <-- 9 sig figs
    'mp_MeV',0 , ... 
    'mn_MeV',0);

pc.hc = pc.h*pc.c; % J*m <-- 6 sig figs.
pc.K = 1/4/pi/pc.epsilonNot;  % Vm/C <-- 
pc.K_eV = 14.3996e-10/pc.e^2; %eVm/C^2
pc.muB = pc.e*pc.hbar/2/pc.me; % J/T
pc.mp = 1.007276466621*pc.u;
pc.mn = 1.00866491595*pc.u;

Joules2eV = 1/pc.e;
pc.h_eV = pc.h*Joules2eV;
pc.hbar_eV = pc.hbar*Joules2eV;
pc.hc_eV = pc.hc*Joules2eV;
pc.me_MeV = pc.me*pc.c^2*Joules2eV/1e6;
pc.mp_MeV = pc.mp*pc.c^2*Joules2eV/1e6;
pc.mn_MeV = pc.mn*pc.c^2*Joules2eV/1e6;

end

% units = struct('c', 'm/s', ...
%         'G', 'm^3/kg/s^2', ...
%         'g', 'm/s^2', ...
%         'h', 'J*s', ...
%         'hbar', 'J*s', ...
%         'hc', 'J*m', ...
%         'muB', 'J/T', ...
%         'epsilonNot','C^2/(N m^2) = farad/meter', ...
%         'muNot', 'N/A^2', ...
%         'K', 'Vm/C = Nm^2/C^2', ...
%         'me', 'kg', ...
%         'mp', 'kg', ...
%         'Na', '1/mol', ...
%         'kb', 'J/K');
%     
%   Joules2eV = 1/pc.e;
%   switch nargin
%     case 1   % Assume 'MeV'
%        pc.h = pc.h*Joules2eV;
%        units.h = 'eV*s';
%        pc.hbar = pc.hbar*Joules2eV;
%        units.hbar = 'eV*s';
%        pc.hc = pc.hc*Joules2eV;
%        units.hc = 'eV*m';
%        pc.me = pc.me*pc.c^2*Joules2eV;
%        units.me = 'eV';
%        pc.mp = pc.mp*pc.c^2*Joules2eV;
%        units.mp = 'eV';
%        pc.K = pc.K*Joules2eV;
%        units.K = 'eVm/C^2';
%        pc.kb = pc.kb*Joules2eV;
%        units.kb = 'eV/k';
%     otherwise
%         % nothing needed, assume units are Joules.
