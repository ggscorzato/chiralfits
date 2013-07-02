function [mpi_ fpi_ BK_] = FSE_BV(mpi,fpi,BK,m11,m1S,mSS,aF0,L,rev)

xi12_m1S=xi_BV(1/2,L(1),sqrt(m1S));
xi12_m11=xi_BV(1/2,L(1),sqrt(m11));
xi32_m11=xi_BV(3/2,L(1),sqrt(m11));

x_M=1./(16*aF0.^2);
x_f=1./(16*aF0.^2);
x_B=1./(4*aF0.^2);

R_M=x_M .* (5 - mSS./m11 + 2*L(1).*sqrt(m11).*(-1 + mSS./m11)) .* xi12_m11;
degnontrivlim=((m11-mSS).*(2*L(1)*sqrt(m11)-1)./(2*m11)) .* xi12_m11;
R_f = x_f .* ((mSS-m11).*xi32_m11 - 8 .* xi12_m1S + degnontrivlim);
R_B = - x_B .* (xi12_m11 - mSS .* xi32_m11); 

%mpi_=FSE_GL(mpi,fpi,aF0,L,0,rev);
mpi_=mpi.* (1 + rev*R_M);

if(~isempty(fpi))
fpi_ = fpi .* (1 + rev*R_f);
end
if(~isempty(BK))
BK_ = BK .* (1 + rev*R_B);
end

function xi = xi_BV(s,L,M)
prefact=3*sqrt(pi)/(gamma(s)*(2*pi)^(3/2));
xi=prefact * exp(-M.*L) .* (2*M.^2).^(3/2-s)./(M.*L).^(2-s);
