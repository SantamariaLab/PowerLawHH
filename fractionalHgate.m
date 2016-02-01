
%% this funcion is used in fractional neuron integration. It integrates the fractional derivative and  the  voltage v at each time t.

function out=fractionalHgate(NetProp,V,t,alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncells=NetProp.Ncells;
dt=NetProp.dt;
Cm=NetProp.Cm;
v0=NetProp.v0;
vrest=NetProp.vrest;
gK=NetProp.gK;
gNa=NetProp.gNa;
gL=NetProp.gL;
EK=NetProp.EK;
ENa=NetProp.ENa;
EL=NetProp.EL;
m=NetProp.m;
h=NetProp.h;
n=NetProp.n;

Namp=NetProp.Noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=V;%vrest.*ones(length(t),Ncells);
dVdt=zeros(length(t)-1,Ncells);
dVdtnow =dVdt(1,:);


% mV=m*ones(length(t),Ncells);
% hV=h*ones(length(t),Ncells);
nV=h*ones(length(t),Ncells);

VMemory=zeros(length(t)-1,Ncells);
Ngatememory=zeros(length(t)-1,Ncells);

I_HionicV=0*ones(length(t),Ncells);
INaV=0*ones(length(t),Ncells);
IKV=0*ones(length(t),Ncells);
ILV=0*ones(length(t),Ncells);

Noise=Namp*randn(length(t),Ncells);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The weight for the  voltage memory trace of the fractional drivative for
% calculated here for the total time t for faster simulation
NN=length(t);
nn=1:NN-1;
WCoet=(NN+1-nn).^(1-alpha)-(NN-nn).^(1-alpha);

% Iionic=@(v,m,h,n,gbarL,gbarNa,gbarK,vrest,Ena,Ek)((gL*(v-EL)+...
%         gNa*m^3*h*(v-Ena)+...
%         gK*n^4*(v-Ek)));
   
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preT=1;%:length(t)-1;

%it integrates the first step as classical HH
for a=preT;
%     %[v(a+1,1),m,h,n,dVdtnow]=HH_RKfull(v(a,1),m,h,n,gL,gNa,gK,EL,ENa,EK,Iinj(a),Cm,dt,v0); 
%     mV(a+1,:) = m;
%     hV(a+1,:) = h;
%     nV(a+1,:) = n;
%     dVdt(a,:)=dVdtnow;
%     VMemory(a,:)=0;
%     Ngatememory(a,:)=0;
    
   % nf=@(x)(alphan*(1-x)-betan*x);
    nf=@(Hh,V,v0)((0.07*exp(-(V-v0)/20))*(1-Hh)-(1./(exp(3-0.1*(V-v0))+1))*Hh);

    k1=nf(nV(a),V(a),v0);
    k2=nf(nV(a)+(dt/2)*k1,V(a),v0);
    k3=nf(nV(a)+(dt/2)*k2,V(a),v0);
    k4=nf(nV(a)+dt*k3,V(a),v0);
    nV(a+1)=nV(a)+(dt/6)*(k1+2*k2+2*k3+k4);
end


kr = dt^alpha*gamma(2-alpha);     %  the kernel   from the fractional derivative and  weighted  the markovian term

% mf=@(Mm,V)(((2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1))*(1-Mm)-(4.*exp(-(V-v0)./18))*Mm);
% hf=@(Hh,V)((0.07*exp(-(V-v0)/20))*(1-Hh)-(1./(exp(3-0.1*(V-v0))+1))*Hh);
% f=@(v,m,h,n,I)(-(1/Cm)*(gL*(v-EL)+gNa*m^3*h*(v-ENa)+gK*n^4*(v-EK))+I);
% 

DeltaN=diff(nV,1,1);
for a=(preT(end)+1):length(t)-1
    
    
%         %%%%% The weight of the memory trace
%         WCoe=WCoet(end-a+2:end);  % The weight for the  voltage memory trace of the fractional drivative  at each  tiime t        
%         %%% Memory trace for  for gating variable n
%         DeltaN =nV(2:a,:)-nV(1:a-1,:);
%     NgateMemory=WCoe*DeltaN;
%     Ngatememory(a)= NgateMemory;
%nf=@(Hh,V,v0)((0.07*exp(-(V-v0)/20))*(1-Hh)-(1./(exp(3-0.1*(V-v0))+1))*Hh);
    
    DeltaN(a-1)=diff(nV(a-1:a,:),1,1); 
    NgateMemory=WCoet((NN-a+1):(NN-1))*DeltaN(1:a-1,:);
    Ngatememory(a)= NgateMemory;
    
    alphan=(0.07*exp(-(V(a)-v0)/20));
    betan=(1./(exp(3-0.1*(V(a)-v0))+1));
    
    n = kr*(alphan*(1-nV(a))-betan*nV(a))+nV(a)-NgateMemory ;
    
    nV(a+1,:) = n;
    
%     if ~(a*dt-round(a*dt))
%         toc
%         plot(nV)
%         a*dt
%         tic
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
out.v=v;
%out. VMemory= VMemory;
out.t=t;
%out.fraccalcu=fraccalcu;
%out.I_Hionic=0;%I_Hionic;
%out.INa=INa;
%out.IK=IK;
%out.IL=IL;
%out.mV=mV;
out.nV=nV;
%out.hV=hV;

out.Ngatememory= Ngatememory;
%out.dVdt= dVdt;
end


function [vnew,mnew,hnew,dVdtnew]=HH_RK(vold,mold,h1old,nold,dt,Iapp)

v0=-65;
V=vold;
alpham=(2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1);
betam=4.*exp(-(V-v0)./18);
alphah=0.07*exp(-(V-v0)/20);
betah=1./(exp(3-0.1*(V-v0))+1);
mf=@(Mm)(alpham*(1-Mm)-betam*Mm);
hf=@(Hh)(alphah*(1-Hh)-betah*Hh);


%%%%%%%%
k1=mf(mold);
k2=mf(mold+(dt/2)*k1);
k3=mf(mold+(dt/2)*k2);
k4=mf(mold+dt*k3);
mnew=mold+(dt/6)*(k1+2*k2+2*k3+k4);

k1=hf(h1old);
k2=hf(h1old+(dt/2)*k1);
k3=hf(h1old+(dt/2)*k2);
k4=hf(h1old+dt*k3);
hnew=h1old+(dt/6)*(k1+2*k2+2*k3+k4);

%%%%%%%%%%%

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1); % channel reversal potentials: mV
Cm=1;

f=@(v,m,h,n,I)(-(1/Cm)*(gL*(v-EL)+gNa*m^3*h*(v-ENa)+gK*n^4*(v-EK))+I);


k1=f(vold,mold,h1old,nold,Iapp);
k2=f(vold+(dt/2)*k1,mold,h1old,nold,Iapp);
k3=f(vold+(dt/2)*k2,mold,h1old,nold,Iapp);
k4=f(vold+dt*k3,mold,h1old,nold,Iapp);
vnew=vold+(dt/6)*(k1+2*k2+2*k3+k4);

dVdtnew=k1;

end





function [vnew,mnew,hnew,nnew,dVdtnew]=HH_RKfull(vold,mold,h1old,nold,gbarL,gbarNa,gbarK,vrest,Ena,Ek,Im,Cm,dt,alpha)

v0=-65;
V=vold;
alpham=(2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1);
betam=4.*exp(-(V-v0)./18);
alphah=0.07*exp(-(V-v0)/20);
betah=1./(exp(3-0.1*(V-v0))+1);
alphan=(0.1-0.01*(V-v0))./(exp(1-0.1*(V-v0))-1);
betan=0.125.*exp(-(V-v0)./80);


mf=@(Mm)(alpham*(1-Mm)-betam*Mm);
hf=@(Hh)(alphah*(1-Hh)-betah*Hh);
nf=@(Nna)(alphan*(1-Nna)-betan*Nna);
f=@(v,m,h,n,gbarL,gbarNa,gbarK,vrest,Ena,Ek,I,Cm)(-(1/Cm)*(gbarL*(v-vrest)+...
        gbarNa*m^3*h*(v-Ena)+...
        gbarK*n^4*(v-Ek))+...
        I);
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=f(vold,mold,h1old,nold,gbarL,gbarNa,gbarK,vrest,Ena,Ek,Im,Cm);
k2=f(vold+(dt/2)*k1,mold,h1old,nold,gbarL,gbarNa,gbarK,vrest,Ena,Ek,Im,Cm);
k3=f(vold+(dt/2)*k2,mold,h1old,nold,gbarL,gbarNa,gbarK,vrest,Ena,Ek,Im,Cm);
k4=f(vold+dt*k3,mold,h1old,nold,gbarL,gbarNa,gbarK,vrest,Ena,Ek,Im,Cm);
vnew=vold+(dt/6)*(k1+2*k2+2*k3+k4);
dVdtnew =k1; 


k1=mf(mold);
k2=mf(mold+(dt/2)*k1);
k3=mf(mold+(dt/2)*k2);
k4=mf(mold+dt*k3);
mnew=mold+(dt/6)*(k1+2*k2+2*k3+k4);

k1=hf(h1old);
k2=hf(h1old+(dt/2)*k1);
k3=hf(h1old+(dt/2)*k2);
k4=hf(h1old+dt*k3);
hnew=h1old+(dt/6)*(k1+2*k2+2*k3+k4);


k1=nf(nold);
k2=nf(nold+(dt/2)*k1);
k3=nf(nold+(dt/2)*k2);
k4=nf(nold+dt*k3);
nnew=nold+(dt/6)*(k1+2*k2+2*k3+k4);


end