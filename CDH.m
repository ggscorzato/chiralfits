function [mpiFV fpiFV report] = CDH(ampiV,afpiV,aF0,aLamb1,aLamb2,aLamb3,aLamb4,a_fm,L,print,rev,parm)
% When rev=1 (-1), compute the finite (infinite) volume am_pi and af_pi from the infinite (finite) volume ones.
% rev=1 corresponds to the formulae in hep-lat/0503014. Equations and Table numbers refer to that paper.
% The expansion parameter is 1/aF0, which must be set to 1/af_pi if these formula are to be used beyond LO.
% It is possible to change the default value of the parameters \tilda{r}_i (i=1,6) with the variable param. 
% The variables mpiFV fpiFV ampiV,afpiV,aF0, aLamb_i are in lattice units. 
% "a_fm" is the lattice spacing in fm and it is used only where is necessary (not at LO).
% L is the number of points in one spatial direction (we assume the spatial volume V=L^3)

gg=[2-pi/2, pi/4-0.5, 0.5-pi/8, 3*pi/16 - 0.5]; % Eq.(62)
mm=[6 12 8 6 24 24 0 12 30 24 24 8 24 48 0 6 48 36 24 24]; % Tab.1
Omm=length(mm);
N=16 * pi^2;
amrho_phys=a_fm*770/197.3; % physical mass of the rho in lattice units

lb1=2*log(aLamb1./ampiV); % related to those of tab.2a by Eq.(53)
lb2=2*log(aLamb2./ampiV);
lb3=2*log(aLamb3./ampiV);
lb4=2*log(aLamb4./ampiV);
lpi=2*log(ampiV./amrho_phys);
rt =[-1.5  3.2 -4.2 -2.5  3.8 1.0]; % tab. 2b
rtr=[-1.5  3.2 -4.2 -2.5  1.0 0.1]; 
if(isempty(parm))
  parm=zeros(1,6);
end
rt = rt + parm .* rtr;
rt=rt';

% choose the quantity (in future I can extend to m_K, etc...)
M_P = ampiV;
if(~isempty(afpiV))
F_P = afpiV;
end
%%%%%%%%
xi_P = (M_P ./(4 * pi * aF0)).^2; % Eq.(10)

for jj=1:length(ampiV)
 lambda_pi = ampiV(jj) .* L(jj); % Eq.(11)
 z = sqrt([1:Omm]')*lambda_pi; %argument of unctions in Eq.(50). sqrt(n) comes from Eq.(26-27)
 B0=2*besselk(1,z)./z; % Eq. (51) one further denominator z comes from Eq.(26-27)
 B2=2*besselk(2,z)./(z.^2); % "
 mmB0(jj)=mm*B0; % remaining factor from Eq.(26-27) and sum, ...
 mmB2(jj)=mm*B2; % which I can already do since all the dependence on n is here
end

% simplifyed S's: Eq.(59)
S4mpi=(13/3)*gg(1) * mmB0 - (1/3)*(40*gg(1) + 32*gg(2) + 26*gg(3))*mmB2;
S6mpi=0; % they say it is negligeable and do not provide a simplified expression.
S4fpi=(1/6)*(8*gg(1) - 13*gg(2))* mmB0 - (1/3)*(40*gg(1) - 12*gg(2) - 8*gg(3) - 13*gg(4))*mmB2;

% Eq. (49) and (54)
I2mpi = -mmB0;
I4mpi = mmB0.*(-55/18 + 4*lb1 + 8/3*lb2 - 5/2*lb3 -2*lb4) + ...
	mmB2.*(112/9 - (8/3)*lb1 - (32/3)*lb2) + S4mpi;
I6mpi = mmB0.*(10049/1296 - 13/72*N + 20/9*lb1 - 40/27*lb2 - 3/4*lb3 - 110/9*lb4 +...
	       - 5/2*lb3.^2 - 5*lb4.^2 +...
               + lb4.*(16*lb1 + 32/3*lb2 - 11*lb3) +...
	       + lpi.*(70/9*lpi + 12*lb1 + 32/9*lb2 - lb3 + lb4 + 47/18) +...
	       + 5*rt(1) + 4*rt(2) + 8*rt(3) + 8*rt(4) + 16*rt(5) + 16*rt(6)) +...
	mmB2.*(3476/81 - 77/288*N + 32/9*lb1 + 464/27*lb2 + 448/9*lb4 +...
	       - 32/3.*lb4.*(lb1+4*lb2) + lpi.*(100/9*lpi + 8/3*lb1 + 176/9*lb2 - 248/9) +...
	       - 8*rt(3) - 56*rt(4) - 48*rt(5) + 16*rt(6))+...
	S6mpi;

if(~isempty(afpiV))
I2fpi = -2*mmB0;
I4fpi = ...
    mmB0.*(-7/9 + 2*lb1 + (4/3)*lb2 - 3*lb4) + ...  % -77
    mmB2.*(112/9 - (8/3)*lb1 -(32/3)*lb2) + ...     % -35
    S4fpi;   % -12
I6fpi = 0; % NNLO not available for fpi.
end

% Eq. (26-27). The sum on n is already done
Rmpi= - (xi_P/2) .* (ampiV./M_P) .* (I2mpi + xi_P .* I4mpi + xi_P.^2 .* I6mpi);
if(~isempty(afpiV))
Rfpi=   (xi_P)   .* (afpiV./F_P) .* (I2fpi + xi_P .* I4fpi + xi_P.^2 .* I6fpi);
end

mpiFV= ampiV .* (1 + rev* Rmpi);
if(~isempty(afpiV))
fpiFV= afpiV .* (1 + rev* Rfpi);
end


% print out some further information (i.e. the contribuition of each order)
if (print == 1)
mlo=- (xi_P/2) .* (ampiV./M_P) .* (I2mpi);
mnlo=- (xi_P/2) .* (ampiV./M_P) .* (xi_P .* I4mpi);
mnnlo=- (xi_P/2) .* (ampiV./M_P) .* (xi_P.^2 .* I6mpi);
flo= (xi_P) .* (afpiV./F_P) .* (I2fpi);
fnlo= (xi_P) .* (afpiV./F_P) .* (xi_P .* I4fpi);
['percent fin V corr: Mpi_LO   Mpi_NLO   Mpi_NNLO    Fpi_LO   Fpi_NLO']
[mlo', mnlo', mnnlo', flo', fnlo']
%report=[(ampiV.*L)', ampiV', Rmpi', sqrt(2)*afpiV', Rfpi'];
end