function [mpiFV fpiFV report] = FSE_GL(mpiV,fpiV,aF0,L,print,rev)
mm=[6 12 8 6 24 24 0 12 30 24 24 8 24 48 0 6 48 36 24 24];
mm4=4*mm;
Nf=2;
x=(mpiV./(4*pi*aF0)).^2;

lambda_pi = mpiV .* L;
z = sqrt([1:20]')*lambda_pi;
g1=mm4 * (besselk(1,z)./z);
%g1=mm4 * (sqrt(pi./(2*z)) .* exp(-z) ./ z);
Rmpi=x .* g1 /(2*Nf);
Rfpi=- x .* g1 *Nf/2;
mpiFV = mpiV .* (1 + rev*Rmpi);
if(~isempty(fpiV))
fpiFV = fpiV .* (1 + rev*Rfpi);
end

if (print == 1)
['percent fin V corr: Mpi_LO  Fpi_LO']
[Rmpi',  Rfpi']
report=[(mpiV.*L)', mpiV', Rmpi', sqrt(2)*fpiV', Rfpi'];
end