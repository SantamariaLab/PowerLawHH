% Simulation of power-law dynamic gate variables in the Hodgkin-Huxley
% model using fractional order derivatives.
% Teka W, Stockton D, Santamaria F. "Power-law dynamics of membrane
% conductances increase spiking diversity in a Hdgkin-Huxley model" PLoS
% Computational Biology, in press, 2016.
% If you use this software please reference our paper. 

% Similarly to the Tekaetaal_SignleGate, sections 7-9 generate a few action potentials. 
%Sections 10-12 analyze these data which corresponds to Figure 2. 
% The other figures were generate by running sections 7-9 for very long
% periods of time. They require super-computer resources and multiple days
%of simulation time. You can request the data files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The full Hodgkin-Huxley model with individual gates having power-law
% dynamics. We model only a few miliseconds to classify the spike shapes.
% Sections 7, 8, and 9 run the models for power-law N, H, and M. Section
% 10, 11, 12 analyze the respective data.
% We provide the data files in case you don't want to wait several days for
% the simulations to run. 
% In sections 7-10 we vary the current (0-24 nA) and the values of the fractional
% order derivative (0.2-1.0). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 7 - Modeling Spiking N
clear
homedir='.'; %this is your directory
cd(homedir)
tstop  = 100; %ms
tstep  = 1e-3; %ms
mVinit = -65; %mV

Ncells = 1; %only modeling one cell
Cm=1;  %microF/cm^2
dt =tstep;
t=0:dt:tstop; %ms
v0 = mVinit*ones([1 Ncells]);              % mV  initial value
vrest=v0(1,1);% mV  the resting  potential.

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2 -- (micro A/mV)/cm^2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1);  % channel reversal potentials: mV
m=0.0529;h=0.596; n=0.3177; % Initial vlaues:  steady state for 0 input n=0.3177, m=0.0529, h=0.596, v=-64.9997

NetProp.Ncells=Ncells;
NetProp.dt=dt;
NetProp.Cm=Cm;
NetProp.v0=v0;
NetProp.vrest=vrest;
NetProp.gK=gK;
NetProp.gNa=gNa;
NetProp.gL=gL;
NetProp.EK=EK;
NetProp.ENa=ENa;
NetProp.EL=EL;
NetProp.m=m;
NetProp.h=h;
NetProp.n=n;
NetProp.Noise=0;


Iamp_1=[1:12];
Iamp_2=[13:24];
Iamp=[Iamp_1 Iamp_2];

I=ones(length(t),Ncells);
I(1:10/dt,Ncells)=0;
f2safe='HH_fracN_firingrateVI';
save(f2safe)

matlabpool('open',12) %use if you have parallel matlab
tic
parfor b=1:length(Iamp)
    out02(b)=runNetworkderivativeHHFractionalPotassium(NetProp,Iamp(b)*I,t,0.2);
end
toc
fprintf('Done with 02 \n')
save(f2safe,'out02','-append');
clear out02

tic
parfor b=1:length(Iamp)    
    out04(b)=runNetworkderivativeHHFractionalPotassium(NetProp,Iamp(b)*I,t,0.4);
end
fprintf('Done with 04 \n')
save(f2safe,'out04','-append');
clear out04
toc

tic
parfor b=1:length(Iamp)    
    out06(b)=runNetworkderivativeHHFractionalPotassium(NetProp,Iamp(b)*I,t,0.6);
end
fprintf('Done with 06 \n')
save(f2safe,'out06','-append');
clear out06
toc

tic
parfor b=1:length(Iamp)
    out08(b)=runNetworkderivativeHHFractionalPotassium(NetProp,Iamp(b)*I,t,0.8);
    
end
fprintf('Done with 08 \n')
save(f2safe,'out08','-append');
clear out08
toc

% we don't run this because we used the one from the H gate. This is the
% classic Hodgkin-Huxley model
% parfor b=1:length(Iamp)
%     out10(b)=runNetworkderivativeHHFractionalPotassium(NetProp,Iamp(b)*I,t,1.0);
%     
% end
% fprintf('done\n');
% save(f2safe,'out10','-append');
% clear out10

matlabpool('close')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 8- Modeling spiking fractional M
clear
homedir='/home/fidel/Wondimu';
cd(homedir)
tstop  = 100; %ms
tstep  = 1e-3; %ms
mVinit = -65; %mv

Ncells = 1;
Cm=1;  %microF/cm^2
dt =tstep;
t=0:dt:tstop; %ms
v0 = mVinit*ones([1 Ncells]);              % mV  initial value
vrest=v0(1,1);% mV  the resting  potential.

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2 -- (micro A/mV)/cm^2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1);  % channel reversal potentials: mV
m=0.0529;h=0.596; n=0.3177; % Initial vlaues:  steady state for 0 input n=0.3177, m=0.0529, h=0.596, v=-64.9997


NetProp.Ncells=Ncells;
NetProp.dt=dt;
NetProp.Cm=Cm;
NetProp.v0=v0;
NetProp.vrest=vrest;
NetProp.gK=gK;
NetProp.gNa=gNa;
NetProp.gL=gL;
NetProp.EK=EK;
NetProp.ENa=ENa;
NetProp.EL=EL;
NetProp.m=m;
NetProp.h=h;
NetProp.n=n;
NetProp.Noise=0;

Iamp_1=[1:12];
Iamp_2=[13:24];
Iamp=[Iamp_1 Iamp_2];



I=ones(length(t),Ncells);
I(1:10/dt,Ncells)=0;
f2safe='HH_fracM_firingrateVI';
save(f2safe)
matlabpool('open',12)
tic
parfor b=1:length(Iamp)
    out02(b)=runNetworkderivativeHHFractionalNa_m(NetProp,Iamp(b)*I,t,0.2);
end
toc
save(f2safe,'out02','-append');
clear out02
fprintf('Done with 02 \n')
parfor b=1:length(Iamp)    
    out04(b)=runNetworkderivativeHHFractionalNa_m(NetProp,Iamp(b)*I,t,0.4);
end
fprintf('Done with 04 \n')

save(f2safe,'out04','-append');
clear out04

parfor b=1:length(Iamp)    
    out06(b)=runNetworkderivativeHHFractionalNa_m(NetProp,Iamp(b)*I,t,0.6);
end
fprintf('Done with 06 \n')
save(f2safe,'out06','-append');
clear out06
parfor b=1:length(Iamp)
    out08(b)=runNetworkderivativeHHFractionalNa_m(NetProp,Iamp(b)*I,t,0.8);
end
fprintf('Done with 08 \n')
save(f2safe,'out08','-append');
clear out08

% parfor b=1:length(Iamp)
%     out10(b)=runNetworkderivativeHHFractionalNa_m(NetProp,Iamp(b)*I,t,1.0);
%     
% end
% 
% fprintf('done\n');
% save(f2safe,'out10','-append');
% clear out10

matlabpool('close')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Section 9 - Modeling Spiking fractional H
clear
homedir='.';
cd(homedir)
tstop  = 100;
tstep  = 1e-3;
mVinit = -65;

Ncells = 1;
Cm=1;  %microF/cm^2
dt =tstep;
t=0:dt:tstop; %ms
v0 = mVinit*ones([1 Ncells]);              % mV  initial value
vrest=v0(1,1);% mV  the resting  potential.

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2 -- (micro A/mV)/cm^2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1);  % channel reversal potentials: mV
m=0.0529;h=0.596; n=0.3177; % Initial vlaues:  steady state for 0 input n=0.3177, m=0.0529, h=0.596, v=-64.9997


NetProp.Ncells=Ncells;
NetProp.dt=dt;
NetProp.Cm=Cm;
NetProp.v0=v0;
NetProp.vrest=vrest;
NetProp.gK=gK;
NetProp.gNa=gNa;
NetProp.gL=gL;
NetProp.EK=EK;
NetProp.ENa=ENa;
NetProp.EL=EL;
NetProp.m=m;
NetProp.h=h;
NetProp.n=n;
NetProp.Noise=0;

alphavalue=[1 0.8 0.4 0.2 ];
Iamp_1=[1:12];
Iamp_2=[13:24];
Iamp=[Iamp_1 Iamp_2];


b=1;

I=ones(length(t),Ncells);
I(1:10/dt,Ncells)=0;
f2safe='HH_fracH_firingrateVI';
save(f2safe)
matlabpool('open',12)
tic
parfor b=1:length(Iamp)
    out02(b)=runNetworkderivativeHHFractionalNa_h(NetProp,Iamp(b)*I,t,0.2);
end
toc
save(f2safe,'out02','-append');
clear out02
fprintf('Done with 02 \n')
parfor b=1:length(Iamp)    
    out04(b)=runNetworkderivativeHHFractionalNa_h(NetProp,Iamp(b)*I,t,0.4);
end
fprintf('Done with 04 \n')

save(f2safe,'out04','-append');
clear out04

parfor b=1:length(Iamp)    
    out06(b)=runNetworkderivativeHHFractionalNa_h(NetProp,Iamp(b)*I,t,0.6);
end
fprintf('Done with 06 \n')
save(f2safe,'out06','-append');
clear out06
parfor b=1:length(Iamp)
    out08(b)=runNetworkderivativeHHFractionalNa_h(NetProp,Iamp(b)*I,t,0.8);
end
fprintf('Done with 08 \n')
save(f2safe,'out08','-append');
clear out08

parfor b=1:length(Iamp)
    out10(b)=runNetworkderivativeHHFractionalNa_h(NetProp,Iamp(b)*I,t,1.0);
    
end

fprintf('done\n');
save(f2safe,'out10','-append');
clear out10

matlabpool('close')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 10 - Analyze full HH model with power-law N
clear
homedir='.';
cd(homedir)

load HH_fracN_firingrateVI 
%the alpha=1 is the same for all the simulations, it is the classic HH
%model
load HH_fracH_firingrateVI out10 
fs={'02','04','06','08','10'};
V2u=[1:length(Iamp)];

for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    spSh{cc1}=[];cc2=1;
    for b=V2u %go through the currents
        t=o2a(b).t;
        mTh=(o2a(b).v>-20); %detect if there's a spike
        sp=find(diff(mTh)==-1);
        spN(b,cc1)=numel(sp);
        fr(b,cc1)=(numel(sp))/90e-3;
        if numel(sp)>1 %record the first evoked spike for each input current
            if (sp(2)+7000)<=length(o2a(b).v)
                [~,mpr]=max(o2a(b).v(sp(2)+[-2500:7000]));
                maxPeak=sp(1)+mpr-2000;
                spSh{cc1}(:,cc2)=o2a(b).v(maxPeak+[-2500:7000])';
                cc2=cc2+1;
            end
        end
    end
end

clf

% Plot the spike shape of the second generated action potential for all
% values of the fractional exponent at the input current values just above 
% rheobase.

cv=['bgrcm'];
etav=[0.2:0.2:1]';
ppv=[1:5];
tsp=[-2500:7000]*dt;
for a=1:5
    subplot(4,4,1)
    %h=plot(tsp,mean(spSh{a},2),cv(a));
    h=plot(tsp,spSh{a}(:,1),cv(a));
    formatFig(h,1)
    hold on
    axis([-2 7 -80 40])
end

%Plot the current at which there was at least 2 action potentials
%generated.
subplot(4,4,5)
spTh=(spN.*~(spN==1));
for a=1:5
    dummy=find(spTh(:,a));
    Ithv(a)=Iamp(dummy(1));
    plot(etav(a),Ithv(a),[cv(a) 'o'])
    hold on
end
h=plot([0.2:0.2:1],Ithv,'k--');
formatFig(h,1)
axis([0 1 0 15]);
nst.etavI=etav;
nst.ITh=Ithv;


% now the action potential threshold using phase plane
pcV=[9:16];
for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    cc2=1;
    for b=V2u
        t=o2a(b).t;
        mTh=(o2a(b).v>-20);
        sp=find(diff(mTh)==-1);
        spN(b,cc1)=numel(sp);
        fr(b,cc1)=(numel(sp))/90e-3;
        if numel(sp)>1
            data=o2a(b).v(sp(2)-4000:end);
            v=data(1:end-1);
            dvdt=diff(data)/1e-3;
            w20=((dvdt<21) .* (dvdt>19));
            vthw=(v<-30);
            th1=(w20.*vthw);
            be=find(diff(th1));
            be2=diff(reshape(be,2,length(be)/2));
            be3=be(1:2:end)+be2';
            vthV=v(be3);
            dvthV=dvdt(be3);
            
            figure(1)
            subplot(4,4,9);%pcV(cc1))
           % h=plot(v(1:10:end),dvdt(1:10:end),'k');
            h=plot(spSh{cc1}(1:50:(end-50),1),diff(spSh{cc1}(1:50:end,1),1,1)/dt/50,cv(cc1));
            formatFig(h,1)
            hold on
            plot(vthV,dvthV,'rs')
            mvthV(cc2)=mean(vthV);
            cc2=cc2+1;
        end
    end
    meanspTh(cc1)=mean(mvthV);
end
axis([-80 50 -200 700])

%Plot the voltage threshold as a function of the fractional order
%derivative
subplot(4,4,13)
for a=1:5
    h=plot(etav(a),meanspTh(a),[cv(a) 'o']);
    hold on
end
h=plot(etav,meanspTh,'k');
formatFig(h,1)
axis([0 1.0 -55 -40])
nst.etav=etav;
nst.meanspTh=meanspTh;

%save  VThCompare nst -append %just when you want to compare to the other
%analyses.


print -depsc2 FracN_APshape.eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 11 - Analyse full Hodgkin-Huxley model with power-law M gate

clear
homedir='.';
cd(homedir)

load HH_fracM_firingrateVI
load HH_fracH_firingrateVI out10
%remember that the simulation at fractional order 0.2 were unstable
fs={'04','06','08','10'};
V2u=[1:length(Iamp)];

for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    spSh{cc1}=[];cc2=1;
    for b=V2u
       t=o2a(b).t;
       mTh=(o2a(b).v>0);
       sp=find(diff(mTh)==-1);
       spN(b,cc1)=numel(sp);
       fr(b,cc1)=(numel(sp))/90e-3;
       if numel(sp)>1
         [~,mpr]=max(o2a(b).v(sp(2)+[-2000:7000]));
         maxPeak=sp(2)+mpr-2000;
         spSh{cc1}(:,cc2)=o2a(b).v(maxPeak+[-2000:7000])';    
         cc2=cc2+1;
       end
    end
end

clf
cv=['bgrcm'];
etav=[0.4:0.2:1]';
ppv=[1:5];
tsp=[-2000:7000]*dt;
for a=1:4
    subplot(4,4,1)
    h=plot(tsp,mean(spSh{a},2),cv(a));
    formatFig(h,1)
    hold on
    axis([-2 7 -80 50])
end

subplot(4,4,5)
spTh=(spN.*~(spN==1));
for a=1:4
    dummy=find(spTh(:,a));
    Ithv(a)=Iamp(dummy(1));
    plot(etav(a),Ithv(a),[cv(a) 'o'])
    hold on
end
h=plot([0.4:0.2:1],Ithv,'k-');
formatFig(h,1)
axis([0 1 0 15]);
mst.etavI=etav;
mst.ITh=Ithv;


% now the action potential threshold using phase plane
pcV=[9:16];
for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    cc2=1;
    for b=V2u
        t=o2a(b).t;
        mTh=(o2a(b).v>-20);
        sp=find(diff(mTh)==-1);
        spN(b,cc1)=numel(sp);
        fr(b,cc1)=(numel(sp))/90e-3;
        if numel(sp)>1
            data=o2a(b).v(sp(2)-4000:end);
            v=data(1:end-1);
            dvdt=diff(data)/1e-3;
            w20=((dvdt<21) .* (dvdt>19));
            vthw=(v<-30);
            th1=(w20.*vthw);
            be=find(diff(th1));
            be2=diff(reshape(be,2,length(be)/2));
            be3=be(1:2:end)+be2';
            vthV=v(be3);
            dvthV=dvdt(be3);
            
            figure(1)
            subplot(4,4,9);%pcV(cc1))
           % h=plot(v(1:10:end),dvdt(1:10:end),'k');
            h=plot(spSh{cc1}(1:10:(end-10),1),diff(spSh{cc1}(1:10:end,1),1,1)/dt/10,cv(cc1));
            formatFig(h,1)
            hold on
            plot(vthV,dvthV,'rs')
            mvthV(cc2)=mean(vthV);
            cc2=cc2+1;
        end
    end
    meanspTh(cc1)=mean(mvthV);
end
axis([-80 50 -200 700])

subplot(4,4,13)
for a=1:4
    h=plot(etav(a),meanspTh(a),[cv(a+1) 'o']);
    hold on
end
h=plot(etav,meanspTh,'k');
formatFig(h,1)
axis([0 1.0 -55 -40])
mst.etav=etav;
mst.meanspTh=meanspTh;
%save  VThCompare mst -append 

print -depsc2 FracM_APshape.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 12 - Analize full Hodgkin-Huxley model with power-law H gate
clear
homedir='.';
cd(homedir)

load HH_fracH_firingrateVI
fs={'02','04','06','08','10'};
V2u=[1:length(Iamp)];


for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    spSh{cc1}=[];cc2=1;
    for b=V2u
       t=o2a(b).t;
       mTh=(o2a(b).v>0);
       sp=find(diff(mTh)==-1);
       spN(b,cc1)=numel(sp);
       fr(b,cc1)=(numel(sp))/90e-3;
       if numel(sp)>0
         [~,mpr]=max(o2a(b).v(sp(1)+[-2000:7000]));
         maxPeak=sp(1)+mpr-2000;
         spSh{cc1}(:,cc2)=o2a(b).v(maxPeak+[-2000:7000])';    
         cc2=cc2+1;
       end
    end
end

clf
cv=['bgrcm'];
etav=[0.2:0.2:1]';
ppv=[1:5];
tsp=[-2000:7000]*dt;
for a=1:5
    subplot(4,4,1)
    h=plot(tsp,mean(spSh{a},2),cv(a));
    formatFig(h,1)
    hold on
    axis([-2 7 -80 40])
end

subplot(4,4,5)
spTh=(spN.*~(spN==1));
for a=2:5
    dummy=find(spTh(:,a));
    Ithv(a)=Iamp(dummy(1));
    plot(etav(a),Ithv(a),[cv(a) 'o'])
    hold on
end
h=plot(etav,Ithv,'k--');
formatFig(h,1)
axis([0 1 0 15]);
hst.etavI=etav;
hst.ITh=Ithv;


% now the action potential threshold using phase plane
pcV=[9:16];
for cc1=2:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    cc2=1;
    for b=V2u
        t=o2a(b).t;
        mTh=(o2a(b).v>0);
        sp=find(diff(mTh)==-1);
        spN(b,cc1)=numel(sp);
        fr(b,cc1)=(numel(sp))/90e-3;
        if numel(sp)>1
            data=o2a(b).v(sp(2)-4000:end);
            v=data(1:end-1);
            dvdt=diff(data)/1e-3;
            w20=((dvdt<21) .* (dvdt>19));
            vthw=(v<-30);
            th1=(w20.*vthw);
            be=find(diff(th1));
            be=be(1:2*floor(length(be)/2));
            be2=diff(reshape(be,2,length(be)/2));
            be3=be(1:2:end)+be2';
            vthV=v(be3);
            dvthV=dvdt(be3);
          
            figure(1)
            subplot(4,4,9);%pcV(cc1))
           % h=plot(v(1:10:end),dvdt(1:10:end),'k');
            h=plot(spSh{cc1}(1:10:(end-10),1),diff(spSh{cc1}(1:10:end,1),1,1)/dt/10,cv(cc1));
            formatFig(h,1)
            hold on
            plot(vthV,dvthV,'rs')
            mvthV(cc2)=mean(vthV);
            cc2=cc2+1;
        end
    end
    meanspTh(cc1)=mean(mvthV);
end
axis([-80 50 -200 700])
subplot(4,4,13)

for a=2:5
    h=plot(etav(a),meanspTh(a),[cv(a) 'o']);
    hold on
end
h=plot(etav,meanspTh,'k');
formatFig(h,1)
axis([0 1.0 -55 -40])
hst.etav=etav;
hst.meanspTh=meanspTh;
save  VThCompare hst -append 
print -depsc2 FracH_APshape.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparing thresholds
%This makes a graph that compares the current and voltage thresholds from
%sections 8-10

load VThCompare
clf

subplot(4,4,13)
for a=1:5
    h=plot(hst.etav(a),hst.ITh(a),[cv(a) 'o']);
    hold on
    h=plot(nst.etav(a),nst.ITh(a),[cv(a) 'o']);
end
for a=1:4
    h=plot(mst.etav(a),mst.ITh(a),[cv(a+1) 'o']);
hold on
    
end
h=plot(hst.etav,hst.ITh,'k');
h=plot(mst.etav,mst.ITh,'k');
h=plot(nst.etav,nst.ITh,'k');
formatFig(h,1)


subplot(4,4,14)
for a=1:5
    h=plot(nst.etav(a),nst.meanspTh(a),[cv(a) 'o']);
    hold on
end
for a=2:5
    h=plot(hst.etav(a),hst.meanspTh(a),[cv(a) 'o']);
end

for a=1:4
    h=plot(mst.etav(a),mst.meanspTh(a),[cv(a+1) 'o']);
hold on
    
end

h=plot(hst.etav(2:end),hst.meanspTh(2:end),'k');
h=plot(mst.etav,mst.meanspTh,'k');
h=plot(nst.etav,nst.meanspTh,'k');

formatFig(h,1)
h2=get(h,'parent');
set(h2(1),'xtick',[0.2 0.4 0.8 1])

print -depsc2 ThrehsoldCompare.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

