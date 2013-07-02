function [z_estim, dz, chi2_P, dev, a_lat] = fit_fpi_mpi_mixed(search,myestim,parm)

  % function for fitting m_pi and f_pi from lattice data (using Finite volume chpt formula).
  % author: Luigi Scorzato
  % last updated 31/07/2007 (used for fits of mixed action data).

addpath /Users/luigi/Work/lib/Matlab
global lec lec_flag out_flag onlyren Nf finV fpi_exp_FSE nnlo assumemaxtwist
%%% SCALE (only enters at higher orders)
r0_fm=0.454; % fm
m_rho_phys_ifm=770 / 197.3; % fm^-1

%%%%%%%%%%%%%%%%%%  CHOOSE DATA
load data%.uni
tab = data;

% initial guess and unknown coefficients from hep-ph/0311023: 
%[86.2e-3 0.12 1.2 0.59 1.25]*1000/197.3  fm^-1

%         F_0, 2B_0/Z_P, Lambda1, Lambda2, Lambda3, Lambda4, kM, kF, r_a  r_ZP
lec =  [0.4369  57.2738 0.6082   6.0821   2.9904   6.3355    0   0   1.0  1.0 ];
range= [0.0025  10.0    0.2027   0.3041   5.0684   0.7096    3   3   0.1  0.1 ];

suppa=1e-6; % suppress errors on a_lat (since at the end we will multiply by a again to obtain an independent results)
with_corr=0;
nnlo=0;
onlyren=0; % 0: f_pi(m) and m_pi(m); 1: f_pi(m_pi).
setfinV=1; % 0: inf vol formulae; 1: Gasser-Leutwyler; 2: G-L with 1/Fpi expans; 3: Colangelo-Duerr-Haefli
assumemaxtwist=1;
brute_av=0;
Nf=2;
DBETA=[3.9];
DMU= [0 0.06];
DLAT=[24];
L2Plot=24;


% q_flag = f_pi,m_pi,m_pcac,a_lat
q_flag=[1 1 1 1];

if(onlyren==0)
%  lec_flag=[1 1 0 0 1 1 0 0 1 1];
  lec_flag=[1 1 0 0 1 1 0 0 0 0];
  out_flag = [1 1 0 0];
  plot_flag =  [1 1 0 0];
elseif(onlyren==1)
  lec_flag=[1 0 0 0 0 1 0 0 0 0];
  out_flag = [1 0 0 0];
  plot_flag =  [0 0 0 1];
end
if(length(DBETA)==1)
  lec_flag(9)=0;
  lec_flag(10)=0;
end
if(isempty(parm))
  parm=zeros(size(lec));
end
lec= lec + parm .* range;

%%%%%%%%%%%%%%%%%% INTERPRET DATA
% vol, run, \beta, \kappa_{ud}, \kappa_{cs}, a\mu_{ud}, a\mu_{cs}, N_{conf}
% r_0/a, am_{\pi}, am_{\rho}, am_{K}, am_{D},  am_{Ds}
% r_0m_{\pi}, (r_0m_{\pi})^2, r_0m_{K}, r_0m_{D}, r_0m_{Ds} 
% m_{\pi}/m_{\rho}, m_{\pi}/m_{K}, am_{\chi}^{PCAC}, r_0m_{\chi}^{PCAC}, aF_{\pi}, r_0F_{\pi}

indph=1;
indL=2;
indbeta=3;
indkud=4;
indmu=5;
indr0 =6;
indmpi =8;
indfpi=10;
indmpcacx=12;
indcorrfm=16;
%%%%%%%%%%%%%%%%%%%% SELECT DATA

tab=tab(find(tab(:,indmpi) < Inf ),:);
tab=tab(find(tab(:,indbeta) <= DBETA(end) ),:);
tab=tab(find(tab(:,indbeta) >= DBETA(1) ),:);
tab=tab(find(tab(:,indL) <= DLAT(end) ),:);
tab=tab(find(tab(:,indL) >= DLAT(1) ),:);
tab=tab(find(tab(:,indmu) <= DMU(end) ),:);
tab=tab(find(tab(:,indmu) >= DMU(1) ),:);
ONLYPHASE=0;
finV=setfinV;

%%%%%%%%%%%%%%%%%%%% lists
lisb =sort(tab(:,indbeta));
lisp =sort(tab(:,indph));
lismu =sort(tab(:,indmu ));

lb =lisb([find(diff(lisb )); end]);
lp =lisp([find(diff(lisp )); end]);
lmu= lismu( [find(diff(lismu )); end]);

%%%%%%%%%%%%%%%%%%%% VARIABLES NAMES
beta=tab(:,indbeta)';
phase=tab(:,indph)';
kud=tab(:,indkud)';
mu=tab(:,indmu)';

LL=tab(:,indL)';
r0=tab(:,indr0)';
dr0=tab(:,indr0+1)';

m_pi=tab(:,indmpi)';
dm_pi=tab(:,indmpi+1)';
ZA=1.0;
dZA=0.0;
m_pcac_a=(-1)*phase.*abs(tab(:,indmpcacx)') .* ZA;
dm_pcac_a= tab(:,indmpcacx+1)'; %abs(m_pcac_a) .*sqrt(((tab(:,indmpcacx+1)')./abs(tab(:,indmpcacx)')).^2+(dZA./ZA).^2)
f_pi=tab(:,indfpi)'/sqrt(2);
df_pi=tab(:,indfpi+1)'/sqrt(2);

m_pi2=m_pi.^2;
dm_pi2=2.*m_pi.*dm_pi;

%% we want the lattice spacing (and Z's) to be functions of beta
for ib=1:length(lb)
  list=find(tab(:,indbeta)==lb(ib))';
  list_num=find(tab(:,indr0) < Inf)';
  list_num=intersect(list,list_num);
  if (brute_av==1)
    a_lat(list)=mean(r0_fm ./r0(list_num)); % lattice spacing in fm
    da_lat(list)=(a_lat(list)).*mean(dr0(list_num)./r0(list_num))*suppa;
  else
    [pol, str]=polyfit(mu(list_num).^2,r0(list_num),1);
    [res, dres]= polyval(pol,0,str);
    if (dres == Inf || dres < 1e-8 )
      dres=max(dr0(list_num));
    end
    a_lat(list)=r0_fm ./ res;
    da_lat(list)=a_lat(list).*(dres./res)*suppa;
  end
end

scaleMeV=1; % 
%a_iMeV(list) = 1/(197.3) * a_lat(list)  % a in MeV^-1


tempx=[f_pi;m_pi;m_pcac_a;a_lat];
tempdx=[df_pi;dm_pi;dm_pcac_a;da_lat];
x=tempx(find(q_flag==1),:);
dx=tempdx(find(q_flag==1),:);

t=[mu;beta;LL];

%%%%%%%%%%%%%%%% Covariance Matrix
clear BrittLuecke
for i=1:size(tab,1)
  sigma(:,:,i)= diag([dx(:,i).^2]);
  if(with_corr==1)
    corrfm=tab(:,indcorrfm)';
    sigma(2,3,i)=corrfm(i)*sqrt(sigma(2,2,i)*sigma(3,3,i));
    sigma(3,2,i)=sigma(2,3,i);
  end
end


%%%%%%%%%%%%%%%%% MIN SEARCH and COMPUTE ERRORS with BrittLuecke
if(search==1)
  model=@(p)chptdev(p,x,t,sigma);
  if(isempty(myestim))
    lec_start= lec(find(lec_flag))
  else
    lec_start= myestim';
  end
  opt=optimset('MaxFunEvals',10000,'MaxIter',10000); % default 1000
  [estim,ff,ef] = fminsearch(model, lec_start,opt);
  z_estim=estim';
  if (ef<1)
    ef
  end
  [dev dz chi2_P DXG DZG WW AA] = BrittLuecke(x,t,z_estim,sigma,@chptd,out_flag,4);
else
  if(isempty(myestim))
    z_estim=lec(find(lec_flag))';
  else
    z_estim=myestim;
  end
  [dev dz chi2_P DXG DZG WW AA] = BrittLuecke(x,t,z_estim,sigma,@chptd,out_flag,4);
end

a_fm_out = a_lat(1);
chi2_P;
%y=z_estim; %convert(z_estim,dz,a_lat(1),1);
%dy=dz;
%y=z_estim;
%dy=dz;
%%%%%%%%%%%%%%%%%%%%%%% PLOTS
%%%% colors
ic=1;
indexbeta(ic)=lb(ic);mycol(ic,:)=[ 0, 0, 0]; % 
mywid(ic)=1  ;psty(ic,:)='x';lsty(ic,:)=' -';ic=ic+1;
if (length(lb)>1)
 indexbeta(ic)=lb(ic);mycol(ic,:)=[ 0, 0, 0]; % 
 mywid(ic)=1  ;psty(ic,:)='+';lsty(ic,:)=' -';ic=ic+1;
end

%%%% plots
for ib=1:length(lb)
    jb=find(indexbeta==lb(ib));
    list=find((tab(:,indbeta)==lb(ib)) );
    if(length(list)>0)
      clear x1 t1 y1

      pstart=0.0005;
      pend=max(abs(m_pi(list)));
      npnts=100;
      dp=(pend-pstart)/npnts;
      x1(2,:)=[pstart:dp:pend];
      pend=max(abs(  mu(list)));
      dp=(pend-pstart)/npnts;
      t1(1,:)=[pstart:dp:pend];

      x1(3,:)=0*ones(1,size(x1(2,:),2));
      x1(1,:)=lec(1)*ones(1,size(x1(2,:),2));
      x1(4,:)=a_lat(list(1))*ones(1,size(x1(2,:),2));
      t1(3,:)=L2Plot*ones(1,size(x1(2,:),2));
      t1(2,:)=beta(list(1))*ones(1,size(x1(2,:),2));

% print finite volume informations
now_f=a_lat(list(1))*z_estim(1);
now_l1=a_lat(list(1))*lec(3);
now_l2=a_lat(list(1))*lec(4);
if(onlyren==1)
now_b=0;
now_l3=a_lat(list(1))*lec(5);
now_l4=a_lat(list(1))*z_estim(2);
elseif(onlyren==0)
now_b=a_lat(list(1))*z_estim(2);
now_l3=a_lat(list(1))*z_estim(3);
now_l4=a_lat(list(1))*z_estim(4);
end
%      [m1 f1]=CDH(m_pi,f_pi,f_pi,now_l1,now_l2,now_l3,now_l4,a_lat,LL,1,-1,parm(1:6)); %CDH,1/F
%      [m2 f2]=FSE_GL(m_pi,f_pi,f_pi,LL,1,-1); % GL,1/F
%      [m3 f3]=FSE_GL(m_pi,f_pi,now_f,LL,1,-1); % GL,1/F0

%      [m_pi', m1', m2', m3']
%      [f_pi', f1', f2', f3']*sqrt(2)

%      XX=now_b/(4 * pi * now_f).^2

      finV=0;
      y1iv = chpteq(x1,t1,z_estim);
      finV=setfinV;
      y1 = chpteq(x1,t1,z_estim);

      mpiovfpi_ph=135/92.42;
      if(onlyren==0)
	mpiovfpi=y1iv(2,:)./y1iv(1,:);
	amu_pi=interp1(mpiovfpi,t1(1,:),mpiovfpi_ph,'linear')
	t1_(1,1)=amu_pi;
	x1_(2,1)=nan; % mpi is irrelevant if onlyren=0
      elseif(onlyren==1)
	mpiovfpi=x1(2,:)./y1iv(1,:);
	ampi_pi=interp1(mpiovfpi,x1(2,:),mpiovfpi_ph,'linear')
	afpi_pi=ampi_pi/mpiovfpi_ph;
	t1_(1,1)=nan; % mu is irrelevant if onlyren=1
	x1_(2,1)=ampi_pi;
      end
      t1_(2,1)=beta(list(1));
      t1_(3,1)=L2Plot;
      x1_(1,1)=nan; % fpi is irrelevant if finv=0
      x1_(3,1)=0; %set m_pcac to 0.
      x1_(4,1)=a_lat(list(1));
      finV=0;
      y1_ = chpteq(x1_,t1_,z_estim);
      afpi_ph=y1_(1);
      a_pi=afpi_ph/(92.4/197.3)







      if (out_flag(2)==1)
	yy(2,:)=y1(2,:).^2;
	yyiv(2,:)=y1iv(2,:).^2;
      end

      if(ONLYPHASE==0)
	cpl=[1:length(x1(3,:))];
      elseif(ONLYPHASE==-1)
	cpl=[ceil(length(x1(3,:))/2)+1:length(x1(3,:))];
      elseif(ONLYPHASE==+1)
	cpl=[1:ceil(length(x1(3,:))/2)];
      end


      if(plot_flag(1)==1)
	figh=figure(1);
	hold on;
	box on;
	XF=t1(1,cpl)*scaleMeV;
	YF=yy(2,cpl)*(scaleMeV^2)./XF;
	XD=mu(list)*scaleMeV;
	dXD=0*mu(list)*scaleMeV;
	YD=m_pi(list).^2*(scaleMeV^2)./XD;
	dYD=2*m_pi(list).*dm_pi(list)*(scaleMeV^2)./XD;
	ylabel('(am_{\pi})^2/(a\mu)','FontSize',20)
	xlabel('a\mu','FontSize',20)
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle','--','LineWidth',mywid(jb));
	errorbarxy(XD,YD,dXD,dYD,[],[],mycol(jb,:),mycol(jb,:),psty(jb,:),[],mywid(jb),[])

	figh=figure(1);
	hold on;
	box on;
	XF=t1(1,cpl)*scaleMeV;
	YF=yyiv(2,cpl)*(scaleMeV^2)./XF;
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle',lsty(jb,:),'LineWidth',mywid(jb));


	figh=figure(10);
	hold on;
	box on;
	XF=t1(1,cpl)*scaleMeV;
	YF=yy(2,cpl)*(scaleMeV^2);
	XD=mu(list)*scaleMeV;
	dXD=0*mu(list)*scaleMeV;
	YD=m_pi(list).^2*(scaleMeV^2);
	dYD=2*m_pi(list).*dm_pi(list)*(scaleMeV^2);
	ylabel('(am_{\pi})^2','FontSize',20)
	xlabel('a\mu','FontSize',20)
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle','--','LineWidth',mywid(jb));
	errorbarxy(XD,YD,dXD,dYD,[],[],mycol(jb,:),mycol(jb,:),psty(jb,:),[],mywid(jb),[])

	figh=figure(10);
	hold on;
	box on;
	XF=t1(1,cpl)*scaleMeV;
	YF=yyiv(2,cpl)*(scaleMeV^2);
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle',lsty(jb,:),'LineWidth',mywid(jb));
      end

      if(plot_flag(2)==1)
	figh=figure(2);
	hold on;
	box on;
	XF=t1(1,cpl)*scaleMeV;
	YF=y1(1,cpl)*scaleMeV;
	XD=mu(list)*scaleMeV;
	dXD=0*mu(list)*scaleMeV;
	YD=f_pi(list)*scaleMeV;
	dYD=df_pi(list)*scaleMeV;
	ylabel('af_{\pi}','FontSize',20)
	xlabel('a\mu','FontSize',20)
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle','--','LineWidth',mywid(jb));
	errorbarxy(XD,YD,dXD,dYD,[],[],mycol(jb,:),mycol(jb,:),psty(jb,:),[],mywid(jb),[])

	figh=figure(2);
	hold on;
	box on;
	XF=t1(1,cpl)*scaleMeV;
	YF=y1iv(1,cpl)*scaleMeV;
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle',lsty(jb,:),'LineWidth',mywid(jb));
      end

      if(plot_flag(4)==1)
	figh=figure(4);
	hold on;
	box on;
	ylabel('af_{\pi}','FontSize',20)
	xlabel('(am_{\pi})^2','FontSize',20)

	figh=figure(4);
	hold on;
	XD=m_pi(list).^2*scaleMeV^2;
	dXD=2*dm_pi(list).*(m_pi(list))*scaleMeV^2;
%	XD=m_pi(list)*scaleMeV^2;
%	dXD=dm_pi(list).*scaleMeV^2;
	YD=f_pi(list)*scaleMeV;
	dYD=df_pi(list)*scaleMeV;
	errorbarxy(XD,YD,dXD,dYD,[],[],mycol(jb,:),mycol(jb,:),psty(jb,:),[],mywid(jb),[])

	figh=figure(4);
	hold on;
	XF=x1(2,cpl).^2*scaleMeV^2;
%	XF=x1(2,cpl)*scaleMeV^2;
	YF=y1(1,cpl)*scaleMeV;
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle','--','LineWidth',mywid(jb));

	figh=figure(4);
	hold on;
	XF=x1(2,cpl).^2*scaleMeV^2;
%	XF=x1(2,cpl)*scaleMeV;
	YF=y1iv(1,cpl)*scaleMeV;
	plot(XF,YF,'Color',mycol(jb,:),'LineStyle',lsty(jb,:),'LineWidth',mywid(jb));

      end




    end
end
%%%%%%%%%%%%%%%%%%%%%%% Functions
%%%%%%%%%% chptdev
function dev = chptdev(params,x,t,sigma)
global lec lec_flag out_flag
lec(find(lec_flag)) = params(:);
z=lec(find(lec_flag))';

dev = BrittLuecke(x,t,z,sigma,@chptd,out_flag,1);

%%%%%%%%%% chptd
function d = chptd(x1,t1,z1)
global out_flag
y = x1(find(out_flag==1),:);
d=chpteq(x1,t1,z1)-y;

%%%%%%%%%% chpteq
function y = chpteq(x1,t1,z1)
global lec lec_flag out_flag onlyren Nf finV fpi_exp_FSE nnlo assumemaxtwist
lec(find(lec_flag))=z1';
F0    = lec(1); % f0 
B0ZP  = lec(2); % 2 B0/ZP
lamb1 = lec(3); % Lambda1
lamb2 = lec(4); % Lambda2
lamb3 = lec(5); % Lambda3
lamb4 = lec(6); % Lambda4
kM    = lec(7);
kF    = lec(8);
r_a   = lec(9);
r_ZP  = lec(10);
za=1;

mu_ = t1(1,:);
b_  = t1(2,:);
LL_ = t1(3,:);
f_pi_ = x1(1,:);
m_pi_ = x1(2,:);
m_pcac_ = za.*x1(3,:);
a_ = x1(4,:);
if(b_ < 4.00)
%  a_ = r_a *a_;
  B0ZP= r_ZP * B0ZP;
end

%% finite -> infinite volume

if(onlyren==1)
if(finV==0)
  mpi2IV=m_pi_.^2;
elseif(finV==1)
  fexp=a_.*F0;
  M_ = FSE_GL(m_pi_,[],fexp,LL_,0,-1);
  mpi2IV=M_.^2;
elseif(finV==2)
  fexp=f_pi_;
  M_ = FSE_GL(m_pi_,[],fexp,LL_,0,-1);
  mpi2IV=M_.^2;
elseif(finV==3)
  fexp=f_pi_;
  M_ = CDH(m_pi_,[],fexp,a_*lamb1,a_*lamb2,a_*lamb3,a_*lamb4,a_,LL_,0,-1,[]);
  mpi2IV=M_.^2;
end
end

%% pcac correction to the quark mass

if(assumemaxtwist==1)
  Bm2 = (a_ .* B0ZP) .* mu_; 
else
  Bm2 = (a_ .* B0ZP) .* sqrt(abs(m_pcac_).^2 + mu_.^2); %  ASSUMING d ZP/da =0 !!!!!!
end

%% f_pi(m_pi)

if(onlyren==1)
  x=mpi2IV ./ (4 * pi .* a_ .* F0).^2;
  l4 = log((a_ .* lamb4).^2./(mpi2IV));
  fpi_th  = a_ .* F0 .* (1 + x .* l4);
  if (nnlo==1)
    l1 = log((a_ .* lamb1).^2./(mpi2IV));
    l2 = log((a_ .* lamb2).^2./(mpi2IV));
    l3 = log((a_ .* lamb3).^2./(mpi2IV));
    lF2= - 1/720 * (23 + 14 * l1 + 16 * l2 + 6 * l3 - 6 * l4).^2;
    fpi_th = fpi_th + a_ .* F0 .* x.^2 .* (-1/2 * l3 + lF2 + 1/2 * l3 .* l4 + kF);
  end
  mpi_th2=mpi2IV;
end

%%  f_pi(m) and m_pi(m)

if(onlyren==0)
  x = Bm2 ./ (4 * pi .* a_ .* F0).^2;
  l3 = log((a_ .* lamb3).^2./(Bm2));
  l4 = log((a_ .* lamb4).^2./(Bm2));
  mpi_th2 = Bm2 .* (1 - (x/2) .* l3);
  fpi_th  = a_ .* F0 .* (1 + x .* l4);
  if (nnlo==1)
    l1 = log((a_ .* lamb1).^2./(Bm2));
    l2 = log((a_ .* lamb2).^2./(Bm2));
    lM=1/51 * (28 * l1 + 32 * l2 - 9 * l3 + 49);
    lF=1/30 * (14 * l1 + 16 * l2 + 6 * l3 - 6 * l4 + 23);
    mpi_th2=mpi_th2 + Bm2 .* x.^2 .* (17/8 * lM.^2 + kM);
    fpi_th=fpi_th + a_ .* F0 .* x.^2 .* (-5/4 * lF.^2 + kF);
  end
end

%% infinite -> finite volume

if(finV==1)
  Fexp=a_ .* F0;
  [mpi_th fpi_th] = FSE_GL(sqrt(mpi_th2),fpi_th,Fexp,LL_,0,1);
  mpi_th2=mpi_th.^2;
elseif(finV==2)
  Fexp=fpi_th;
  [mpi_th fpi_th] = FSE_GL(sqrt(mpi_th2),fpi_th,Fexp,LL_,0,1);
  mpi_th2=mpi_th.^2;
elseif(finV==3)
  Fexp=fpi_th;
  [mpi_th fpi_th] = CDH(sqrt(mpi_th2),fpi_th,Fexp,a_*lamb1,a_*lamb2,a_*lamb3,a_*lamb4,a_,LL_,0,1,[]);
  mpi_th2=mpi_th.^2;
end

temp=[fpi_th; sqrt(mpi_th2); m_pcac_; a_];
y=temp(find(out_flag==1),:);

