
%% this funcion is used in fractional neuron integration. It integrates the fractional derivative and  the  voltage v at each time t.

function out=fractionalNgate(NetProp,V,t,alpha)

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


mV=m*ones(length(t),Ncells);
hV=h*ones(length(t),Ncells);
nV=zeros(length(t),Ncells);
nV(1)=n;

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
    
    alphan=(0.1-0.01*(v(a)-v0))./(exp(1-0.1*(v(a)-v0))-1);
    betan=0.125.*exp(-(v(a)-v0)./80);
    nf=@(Nna)(alphan*(1-Nna)-betan*Nna);
    k1=nf(nV(a));
    k2=nf(nV(a)+(dt/2)*k1);
    k3=nf(nV(a)+(dt/2)*k2);
    k4=nf(nV(a)+dt*k3);
    nV(a+1)=nV(a)+(dt/6)*(k1+2*k2+2*k3+k4);
 
end


kr = dt^alpha*gamma(2-alpha);     %  the kernel   from the fractional derivative and  weighted  the markovian term

% mf=@(Mm,V)(((2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1))*(1-Mm)-(4.*exp(-(V-v0)./18))*Mm);
% hf=@(Hh,V)((0.07*exp(-(V-v0)/20))*(1-Hh)-(1./(exp(3-0.1*(V-v0))+1))*Hh);
% f=@(v,m,h,n,I)(-(1/Cm)*(gL*(v-EL)+gNa*m^3*h*(v-ENa)+gK*n^4*(v-EK))+I);
% 
DeltaN=diff(nV,1,1);
DeltaN(preT+1)=0;
%DeltaD=DeltaN';
for a=(preT(end)+1):length(t)-1
    
%     %%%%% The weight of the memory trace
%     WCoe=WCoet(end-a+2:end);  % The weight for the  voltage memory trace of the fractional drivative  at each  tiime t
%     %%% Memory trace for  for gating variable n
%     DeltaN =nV(2:a,:)-nV(1:a-1,:);
    %%%%% The weight of the memory trace
   % WCoe=WCoet((NN-a+1):(NN-1));  % The weight for the  voltage memory trace of the fractional drivative  at each  tiime t
    %%% Memory trace for  for gating variable n
    %dummy= [diff(nV(a-1:a,:),1,1); DeltaN];        
    %NgateMemory=WCoe*DeltaN;
    %efficient v1
    %DeltaN(a-1,:)=diff(nV(a-1:a,:),1,1); 
    %NgateMemory=WCoet((NN-a+1):(NN-1))*DeltaN(1:a-1,:);
    
%     tic
%     DeltaD=circshift(DeltaD,-1);
%     DeltaD(NN-1,:)=diff(nV(a-1:a,:),1,1); 
%     NgateMemory=WCoet*DeltaD;
%     clock1=toc;
    
    
    DeltaN(a-1,:)=diff(nV(a-1:a,:),1,1);
    %tic
    NgateMemory=WCoet((NN-a+1):(NN-1))*DeltaN(1:a-1,:);
    %clock2=toc;
    
%     DeltaD(a-1,:)=diff(nV(a-1:a));
%     tic
%     NgateMemory=sum(WCoet((NN-a+1):(NN-1)).*DeltaD(1:a-1));
%     clock1=toc;
%     disp(['clock 1 ' num2str(clock1) ' clock 2 ' num2str(clock2)]);
    
    Ngatememory(a)= NgateMemory;
    
    alphan=(0.1-0.01*(v(a,1)-v0))./(exp(1-0.1*(v(a,1)-v0))-1);
    betan=0.125.*exp(-(v(a,1)-v0)./80);
    
    n = kr*(alphan*(1-nV(a))-betan*nV(a))+nV(a)-NgateMemory ;
    
    nV(a+1,:) = n;
    
%     if ~(a*dt-round(a*dt))
%         %toc
%         %plot(nV)
%         a*dt
%         %tic
%     end
    
    
    
end
out.v=v;
out.t=t;
out.nV=nV;
out.Ngatememory= Ngatememory;
end

