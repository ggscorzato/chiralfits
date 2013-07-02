function LEC_out = rescaleLEC(N_in,N_out,mu_in,mu_out,LEC_in,comb)
% converts LEC_in from SU(N_in) at scale mu_in to SU(N_out) at scale mu_out
% comb is a row vector describing your combination (e.g. 2L_6+L_8-> com=[0 0 0 0 0 2 0 1 0 0].)
% Scales have to be in MeV.
% Refs. are [DGH VI.3.19] for the change of scale and [DGH VI.2.9] for the change of Nf.


%%%%%%%%% 1. first change the scale 

if(N_in==2)
  c= 1/(32 * pi^2)*[1/12,  1/6, 0, 1/8, 1/4,   3/32, 0,    0, 1/6, -1/6];
elseif(N_in==3)
  c= 1/(32 * pi^2)*[3/32, 3/16, 0, 1/8, 3/8, 11/144, 0, 5/48, 1/4, -1/4];
else
  % and for general Nf ??
  %c= 1/(128 * pi^2)*[?,?,0,1/2,Nf/2,(Nf^2+2)/(4*Nf^2),0,(Nf^2-4)/(4*Nf),?,?];
  'i can do only SU(2) or SU(3)'
  return;
end

cc= comb*c';
LEC_out=LEC_in - cc *log((mu_out/mu_in)^2);



%%%%%%%% 2. then convert to new Nf 

if(N_in~=N_out)
  %          1  2  3  4  5  6  7  8  9  10
  base(:,1)=[2  0  1  0  0  0  0  0  0  0 ]';
  base(:,2)=[0  0  0  2  1  0  0  0  0  0 ]';
  base(:,3)=[0  0  0  0  0  2  0  1  0  0 ]';
  base(:,4)=[0  0  0  1  0 -1 -9 -3  0  0 ]';
  base(:,5)=[0  1  0  0  0  0  0  0  0  0 ]';
  base(:,6)=[0  0  0  0  0  0  0  0  1  0 ]';
  base(:,7)=[0  0  0  0  0  0  0  0  0  1 ]';
  
  newcomb= pinv(base) * comb';

  MH=548.8; % MeV
  MK=495; % MeV
  F0=93; % MeV
  nuK=1/(384*pi^2) * (log((MK/mu_out)^2) + 1);
  nuH=1/(384*pi^2) * (log((MH/mu_out)^2) + 1);

  shifts(1)=-nuK/4;
  shifts(2)=-3/2 * nuK;
  shifts(3)=-3/4 * nuK - 1/12 * nuH;
  shifts(4)= 3/2 * nuK + F0^2/(24*MH^2) + 5/(1152 *pi^2)*log((MH/mu_out)^2);
  shifts(5)=-nuK/4;
  shifts(6)=-nuK;
  shifts(7)=+nuK;

  if((N_in==3)&&(N_out==2))
    LEC_out=LEC_out+shifts*newcomb;
  end
  if((N_in==2)&&(N_out==3))
    LEC_out=LEC_out-shifts*newcomb;
  end
end

% 4: 2 3 3 2