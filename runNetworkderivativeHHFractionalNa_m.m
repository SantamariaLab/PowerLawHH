%% this funcion is used in fractional neuron integration. It integrates the fractional derivative and  the  voltage v at each time t.

function out=runNetworkderivativeHHFractionalNa_m(NetProp,Iinj,t,alpha)

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

%=================================
% variables in vecotr form 
v=vrest.*ones(length(t),Ncells);

mV=m*ones(length(t),Ncells);
hV=h*ones(length(t),Ncells);
nV=n*ones(length(t),Ncells);

VMemory=zeros(length(t)-1,Ncells);
Ngatememory=zeros(length(t)-1,Ncells);
Mgatememory=zeros(length(t)-1,Ncells);

I_HionicV=0*ones(length(t),Ncells);
INaV=0*ones(length(t),Ncells);
IKV=0*ones(length(t),Ncells);
ILV=0*ones(length(t),Ncells);


%============================================================
% The weight for the   memory trace of the fractional drivative for
% calculated here for the total time t for faster simulation
NN=length(t);
nn=1:NN-1;
WCoet=(NN+1-nn).^(1-alpha)-(NN-nn).^(1-alpha);

% Iionic=@(v,m,h,n,gbarL,gbarNa,gbarK,vrest,Ena,Ek)((gL*(v-EL)+...
%         gNa*m^3*h*(v-Ena)+...
%         gK*n^4*(v-Ek)));
   
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preT=1;

for a=preT
    [v(a+1,1),m,h,n]=HH_RKfull(v(a,1),m,h,n,gL,gNa,gK,EL,ENa,EK,Iinj(a),Cm,dt,v0);
   
    mV(a+1,:) = m;
    hV(a+1,:) = h;
    nV(a+1,:) = n;

     VMemory(a,:)=0;
    Ngatememory(a,:)=0;
     
    
       Hgatememory(a,:)=0;
        Mgatememory(a,:)=0;
end
%========================================================================================

kr = dt^alpha*gamma(2-alpha);     %  the kernel   from the fractional derivative that  weighted  the markovian term

for a=(preT(end)+1):length(t)-1   
    
   %Iionic =gL*(v(a,1)-EL)+gNa*m^3*h*(v(a,1)-ENa)+gK*n^4*(v(a,1)-EK);
           
    [v(a+1,1),h,n]=HH_RK(v(a,1),m,h,n,gL,gNa,gK,EL,ENa,EK,Iinj(a),Cm,dt,v0);
   
    
    
    %[m,h,n]=HH_RK(v(a,1),m,h,n,dt,1);
    
   % mV(a+1,:) = m;
    hV(a+1,:) = h;
    nV(a+1,:) = n;
    
    
       
    
    %%%%% The weight of the memory trace  
    WCoe=WCoet(end-a+2:end);  % The weight for the  voltage memory trace of the fractional drivative  at each  tiime t
    
    
  % voltage memory trace  
%     TeDi=v(2:a,:)-v(1:a-1,:); % Delta V (using all past values of V)  of  the  voltage memory trace of the fractional drivative  at each  tiime t
%     VoltMemory=WCoe*TeDi;  %  voltage memory trace 
%      VMemory(a,:)= VoltMemory;
%     fraccalcu=VoltMemory-v(a,1);  %  The fraction derivative
    
    %v(a+1,:) = (kr/Cm)*(Iinj(a)- Iionic) + v(a,1) - VoltMemory;
    
   %  v(a+1,:) = (dt/Cm)*(Iinj(a)- Iionic) + v(a,1) ;
     
     
    
    %%% Memory trace for  for gating variable n
    
%    DeltaN =nV(2:a,:)-nV(1:a-1,:);
%    NgateMemory=WCoe*DeltaN;
%    Ngatememory(a,:)= NgateMemory;
%       
%    alphan=(0.1-0.01*(v(a,1)-v0))./(exp(1-0.1*(v(a,1)-v0))-1);
%    betan=0.125.*exp(-(v(a,1)-v0)./80);
%   
%    
%     n = kr*(alphan*(1-nV(a,1))-betan*nV(a,1))+nV(a,1)-NgateMemory;
%    
%     nV(a+1,:) = n;
   


   DeltaM =mV(2:a,:)-mV(1:a-1,:);
   MgateMemory=WCoe*DeltaM;
   Mgatememory(a,:)= MgateMemory;
   
   
  
   
 alpham=(2.5-0.1.*(v(a,1)-v0))./(exp(2.5-0.1.*(v(a,1)-v0))-1);
betam=4.*exp(-(v(a,1)-v0)./18);
 
  
   
    m = kr*(alpham*(1-mV(a,1))-betam*mV(a,1))+mV(a,1)-MgateMemory;
   
    mV(a+1,:) = m;



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
out.mV=mV;
out.nV=nV;
out.hV=hV;

%out.Ngatememory= Ngatememory;
out.Mgatememory= Mgatememory;
end


%function [vnew,mnew,hnew]=HH_RK(vold,mold,h1old,nold,dt,Iapp)
function [vnew,hnew,nnew]=HH_RK(vold,mold,h1old,nold,gL,gNa,gK,EL,ENa,EK,Iapp,Cm,dt,v0)




V=vold;

alphah=0.07*exp(-(V-v0)/20);
betah=1./(exp(3-0.1*(V-v0))+1);
hf=@(Hh)(alphah*(1-Hh)-betah*Hh);




alphan=(0.1-0.01*(V-v0))./(exp(1-0.1*(V-v0))-1);
betan=0.125.*exp(-(V-v0)./80);
nf=@(Nna)(alphan*(1-Nna)-betan*Nna);

 


%%%%%%%%
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




% k1=hf(h1old);
% k2=hf(h1old+(dt/2)*k1);
% k3=hf(h1old+(dt/2)*k2);
% k4=hf(h1old+dt*k3);
% hnew=h1old+(dt/6)*(k1+2*k2+2*k3+k4);

%%%%%%%%%%%



f=@(v,m,h,n,I)(-(1/Cm)*(gL*(v-EL)+gNa*m^3*h*(v-ENa)+gK*n^4*(v-EK))+I);
    
 
 k1=f(vold,mold,h1old,nold,Iapp);
k2=f(vold+(dt/2)*k1,mold,h1old,nold,Iapp);
k3=f(vold+(dt/2)*k2,mold,h1old,nold,Iapp);
k4=f(vold+dt*k3,mold,h1old,nold,Iapp);
vnew=vold+(dt/6)*(k1+2*k2+2*k3+k4);

end





function [vnew,mnew,hnew,nnew]=HH_RKfull(vold,mold,h1old,nold,gbarL,gbarNa,gbarK,EL,Ena,Ek,Im,Cm,dt,v0)


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



f=@(v,m,h,n,gbarL,gbarNa,gbarK,vrest,Ena,Ek,I,Cm)(-(1/Cm)*(gbarL*(v-EL)+...
        gbarNa*m^3*h*(v-Ena)+...
        gbarK*n^4*(v-Ek))+...
        I);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=f(vold,mold,h1old,nold,gbarL,gbarNa,gbarK,EL,Ena,Ek,Im,Cm);
k2=f(vold+(dt/2)*k1,mold,h1old,nold,gbarL,gbarNa,gbarK,EL,Ena,Ek,Im,Cm);
k3=f(vold+(dt/2)*k2,mold,h1old,nold,gbarL,gbarNa,gbarK,EL,Ena,Ek,Im,Cm);
k4=f(vold+dt*k3,mold,h1old,nold,gbarL,gbarNa,gbarK,EL,Ena,Ek,Im,Cm);
vnew=vold+(dt/6)*(k1+2*k2+2*k3+k4);


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
