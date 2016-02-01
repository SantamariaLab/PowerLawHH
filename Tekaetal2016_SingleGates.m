% Simulation of power-law dynamic gate variables in the Hodgkin-Huxley
% model using fractional order derivatives.
% Teka W, Stockton D, Santamaria F. "Power-law dynamics of membrane
% conductances increase spiking diversity in a Hdgkin-Huxley model" PLoS
% Computational Biology, in press, 2016.
% If you use this software please reference our paper. 

%There are 6 sections in this file that compute the simulations and
%analyze the results for Figures 1 and 2 of our paper.
%Section 1-3 simulate voltage clamp on the individual gates having 
%power-law dynamcis. They correspond to figure 1 in the paper. 
% Section 4-6 analyze the data produced by section 1-3.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation section 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 1
% Simulationg power-law dynamics in the N gate of the Hodgking-Huxley model
% corresponding to the potassium conductance.

clear
homedir='.'; %Your home direcoty
cd(homedir)

f2save='fractoinalNgateSweep'; %Name of file

tstop   = 110; %ms
tstep   = 1e-3; %ms
mVinit  = -65; %mV
Iinjamp = -45; %nA


Ncells  = 1; %Model only one cell
Cm      = 1;  %microF/cm^2
dt      = tstep;
t       = 0:dt:tstop; %ms
v0      = mVinit*ones([1 Ncells]); % mV  initial value
vrest   = v0(1,1);% mV  the resting  potential.

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2 -- (micro A/mV)/cm^2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1);  % channel reversal potentials: mV
m=0.0529;h=0.596; n=0.3177; % Initial vlaues:  steady state for 0 input n=0.3177, m=0.0529, h=0.596, v=-64.9997

%Setup the modeling structure NetProp
NetProp.Ncells= Ncells;
NetProp.dt    = dt;
NetProp.Cm    = Cm;
NetProp.v0    = v0;
NetProp.vrest = vrest;
NetProp.gK    = gK;
NetProp.gNa   = gNa;
NetProp.gL    = gL;
NetProp.EK    = EK;
NetProp.ENa   = ENa;
NetProp.EL    = EL;
NetProp.m     = m;
NetProp.h     = h;
NetProp.n     = n; 
NetProp.Noise = 0; %noise added to the input current

alphavalue=[1 0.8 0.4 0.2 ];

%This is the sweep variable. 
Vamp=[-110:10:120];

V=zeros(length(t),length(Vamp));
for a=1:length(Vamp)
    V(:,a)=Vamp(a)*ones(length(t),Ncells);
    V(1:30/dt,a)=mVinit;
end

save(f2save); %Save all the structures used to run the simulations

%if you have parallel matlab then use this, otherwise comment out next line
matlabpool('open',12) 

parfor b=1:length(Vamp)
    out02(b)=fractionalNgate(NetProp,V(:,b),t,0.2); 
end
save(f2save,'out02','-append')
clear out02;
fprintf('N donee 02\n');

parfor b=1:length(Vamp)
    out04(b)=fractionalNgate(NetProp,V(:,b),t,0.4);
end
save(f2save,'out04','-append')
clear out04;
fprintf('N done 04\n');

parfor b=1:length(Vamp)
    out06(b)=fractionalNgate(NetProp,V(:,b),t,0.6);
end
save(f2save,'out06','-append')
clear out06;
fprintf('N done 06\n');

parfor b=1:length(Vamp)
    out08(b)=fractionalNgate(NetProp,V(:,b),t,0.8);
end
save(f2save,'out08','-append')
clear out08;
fprintf('N done 08\n');

parfor b=1:length(Vamp)
    out10(b)=fractionalNgate(NetProp,V(:,b),t,1); % main output
end
save(f2save,'out10','-append')
clear out10;
fprintf('N done 10\n');

matlabpool('close') %only if you have parallel Matlab
save fractoinalNgateSweep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2 - Power-law dynamics of the H gate (Sodium conductance)

clear
homedir='/home/fidel/Wondimu';
cd(homedir)
tstop  = 130;
tstep  = 1e-3;
mVinit = -65;
Iinjamp = -45;
alphav = 1;

Ncells = 1;
Cm=1;  %microF/cm^2
dt =tstep;
t=0:dt:tstop; %ms
v0 = mVinit*ones([1 Ncells]);              % mV  initial value
vrest=v0(1,1);% mV  the resting  potential.

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2 -- (micro A/mV)/cm^2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1);  % channel reversal potentials: mV
m=0.0529;h=0.596; n=0.3177; % Initial vlaues:  steady state for 0 input n=0.3177, m=0.0529, h=0.596, v=-64.9997

Iinjamplitude=Iinjamp; % nA (or microA/cm2) injected current amplitude,
V=Iinjamplitude*ones(length(t),Ncells);
%V(1:60/dt,Ncells)=mVinit;

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

Vamp=[-110:5:-90 -88:2:-40 -35:5:60];
f2save='fractoinalHgateSweep';

V=zeros(length(t),length(Vamp));
for a=1:length(Vamp)
    V(:,a)=Vamp(a)*ones(length(t),Ncells);
    V(1:30/dt,a)=mVinit;
end


save(f2save);
matlabpool('open',12)

tic
parfor b=1:length(Vamp)
    out02(b)=fractionalHgate(NetProp,V(:,b),t,0.2); % main output
end
save(f2save,'out02','-append')
clear out02;
fprintf('H donee 02\n');
toc

tic
parfor b=1:length(Vamp)
    out04(b)=fractionalHgate(NetProp,V(:,b),t,0.4);
end
save(f2save,'out04','-append')
clear out04;
fprintf('H donee 04\n');
toc

tic
parfor b=1:length(Vamp)
    out06(b)=fractionalHgate(NetProp,V(:,b),t,0.6);
end
save(f2save,'out06','-append')
clear out06;
fprintf('H donee 06\n');
toc

tic
parfor b=1:length(Vamp)
    out08(b)=fractionalHgate(NetProp,V(:,b),t,0.8);
end
save(f2save,'out08','-append')
clear out08;
fprintf('H donee 08\n');
toc

tic
parfor b=1:length(Vamp)
    out10(b)=fractionalHgate(NetProp,V(:,b),t,1); % main output
end
save(f2save,'out10','-append')
clear out10;
fprintf('H donee 10\n');
toc

matlabpool('close')

save fractoinalHgateSweep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3 - Power-law dynamics of the M gate (Sodium conductance)

clear
homedir='/home/fidel/Wondimu';
cd(homedir)
tstop  = 1130;
tstep  = 0.5e-3;
mVinit = -65;
Iinjamp = -45;
alphav = 1;

Ncells = 1;
Cm=1;  %microF/cm^2
dt =tstep;
t=0:dt:tstop; %ms
v0 = mVinit*ones([1 Ncells]);              % mV  initial value
vrest=v0(1,1);% mV  the resting  potential.

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2 -- (micro A/mV)/cm^2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1);  % channel reversal potentials: mV
m=0.0529;h=0.596; n=0.3177; % Initial vlaues:  steady state for 0 input n=0.3177, m=0.0529, h=0.596, v=-64.9997

Iinjamplitude=Iinjamp; % nA (or microA/cm2) injected current amplitude,
V=Iinjamplitude*ones(length(t),Ncells);
V(1:20/dt,Ncells)=mVinit;

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

Vamp=[-110:5:60];
f2save='fractoinalMgateSweep';

V=zeros(length(t),length(Vamp));
for a=1:length(Vamp)
    V(:,a)=Vamp(a)*ones(length(t),Ncells);
    V(1:30/dt,a)=mVinit;
end


save(f2save);

matlabpool('open',12)

% These simulation were unstable and were not performed for the project.
% See text for details.
% parfor b=1:length(Vamp)
%     out02(b)=fractionalMgate(NetProp,V(:,b),t,0.2); % main output
% end
% save(f2save,'out02','-append')
% clear out02;
% fprintf('M donee 02\n');
% toc

tic
parfor b=1:length(Vamp)
    out04(b)=fractionalMgate(NetProp,V(:,b),t,0.4);
end
save(f2save,'out04','-append')
clear out04;
fprintf('M donee 04\n');
toc

tic
parfor b=1:length(Vamp)
    out06(b)=fractionalMgate(NetProp,V(:,b),t,0.6);
end
save(f2save,'out06','-append')
clear out06;
fprintf('M donee 06\n');
toc

tic
parfor b=1:length(Vamp)
    out08(b)=fractionalMgate(NetProp,V(:,b),t,0.8);
end
save(f2save,'out08','-append')
clear out08;
fprintf('M donee 08\n');
toc

tic
parfor b=1:length(Vamp)
    out10(b)=fractionalMgate(NetProp,V(:,b),t,1); % main output
end
save(f2save,'out10','-append')
clear out10;
fprintf('M donee 10\n');
toc

matlabpool('close')
save fractoinalMgateSweep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After running the simulations of voltage clamp for each gate we then 
% analyze the reults.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 4 analyze fractional N kinetics
clear
homedir='.';%type your directory
cd(homedir)

load fractoinalNgateSweep

%dual exponential function to fit the gate response.
f=fittype('ninf*(1-exp(-x*tau))+n2*(1-exp(-x*tau2))'); 
f0=fitoptions(f);
g=[];
cc1=1;
%the suffix for the simulations with different values of fractional
%exponent
fs={'02','04','06','08','10'};
V2u=[1:length(Vamp)]; %Vamp was defined in the particular simulations


%All this section does is to automate the fitting process by selecting
%appropriate bounds. If the gate decreased value from the initial condition
%the the min-max bound would be different if the gate value increased for a
%different voltage. 
for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    cc2=1;
    for b=V2u
        g(cc2).tr=(o2a(b).t(logical(o2a(b).t>=30))-30)';
        d2f=o2a(b).nV(logical(o2a(b).t>=30));
        g(cc2).nfinf=d2f(end);
        g(cc2).d2f=d2f-d2f(1);
        g(cc2).d2fraw=d2f;
        if cc1==5
            f0.Start=[(g(cc2).d2f(end))*[1 1] 1 1 ];
        else
            f0.Start=[(g(cc2).d2f(end))*[0.1 0.5 ] 0.002 1 ];
        end
        
        if (max(g(cc2).d2f))>0 && min(g(cc2).d2f)==0
            f0.Upper=[abs(g(cc2).d2f(end))*[1.5 1.5 ] 40 40 ];
            f0.Lower=[0*[1.5 1.5 ] 0 0];
        elseif (max(g(cc2).d2f))>0 && min(g(cc2).d2f)<0
            f0.Upper=[0*[1.5 1.5 ] 40 40 ];
            f0.Lower=[min(g(cc2).d2f)*[1.5 1.5 ] 0 0];
        elseif (max(g(cc2).d2f))==0
            f0.Upper=[0*[1.5 1.5 ] 40 40 ];
            f0.Lower=[-abs(g(cc2).d2f(end))*[1.5 1.5 ] 0 0];
        else 
            disp('WHAaaaat!')
        end

        
        g(cc2).fL=fit(g(cc2).tr,...
            g(cc2).d2f,f,f0);
        tauVdummy=1./[g(cc2).fL.tau g(cc2).fL.tau2];
        ref2E=find(g(cc2).tr>min(tauVdummy));
        if isempty(ref2E)
            ref2E=length(g(cc2).tr);
        end
         [b Vamp(b)]
        clf
        plot(g(cc2).tr,(g(cc2).d2f),'k',...
            g(cc2).tr(ref2E(1)),g(cc2).d2f(ref2E(1)),'r*',...
            g(cc2).tr,g(cc2).fL(g(cc2).tr),'r');
        drawnow
        
        cc2=cc2+1;
    end
    gAll(cc1).g=g;%the selected fit.
end

%Analyze and plot data
clf

%Plot the value of the N gate vs time for V=30 and values of fractional 
%order 0.2, 0.4, 0.6, 0.8, 1.0
subplot(4,4,1)
cla
vamp2p=15;
Vamp(vamp2p)
for a=1:length(gAll)
    trEx4(:,a)=gAll(a).g(vamp2p).tr;
    d2a4r(:,a)=gAll(a).g(vamp2p).d2fraw;
    d2a4(:,a)=gAll(a).g(vamp2p).d2f;
    fL(a).f=gAll(a).g(vamp2p).fL;
    trEx12(:,a)=gAll(a).g(vamp2p).tr;
    d2a12(:,a)=gAll(a).g(vamp2p).d2fraw;
end
h=plot(trEx4(:,1),(d2a4r),'-');
formatFig(h,1)
axis([0 30 0.3 1]);

%Log-log plot of the same data and fitting - Not used in the paper
%Same data as in subplot(4,4,1)
subplot(4,4,5)
cla
h=loglog(trEx4(:,1),(d2a4));
formatFig(h,1)
cla
cv='bgrmc';
for a=1:size(d2a4,2)
    [logfit(a).f, logfit(a).g]=fit(log10(trEx4(2:500,1)),log10(d2a4(2:500,a)),'poly1');
    
    h=plot(log10(trEx4(:,1)),log10(d2a4(:,a)),cv(a));
    hold on
    fracexpfig(a)=logfit(a).f.p1;
end
formatFig(h,1)
axis([-3.2 2 -1 0])


%do an adaptive fit for all the log traces!
for a=1:length(gAll)
    for b=1:length(gAll(a).g)
        g2f=gAll(a).g(b);
        out=adaptiveLogFit(g2f,0.95);
        ftM(a,b)=out.time;% get the average from this.
        t2int(a,b)=10.^((log10(gAll(5).g(b).d2fraw(end))-out.f.p2)./out.f.p1); %he last point with alpha =1
        
    end
end

%plot the exponent of the log fit to the voltage response. This corresponds
% to the fractional order of the derivative - Not used in the paper.
subplot(4,4,9)
h=plot(fracexpfig,'o--k');
formatFig(h,1)

%Plot the value of N_inf. In the paper we call this the long term response
%fuction
subplot(4,4,2)
cla
for cc1=1:length(gAll);
    tg=gAll(cc1).g;
    for b=1:length(tg)
        mginf(cc1,b)=tg(b).nfinf;
    end
end
h=plot(Vamp(V2u),mginf','-');
formatFig(h,1);
axis([-110 20 0 1])

%now calculate the slope at the half max value - Not used in the paper
for a=1:size(mginf,1)
    
    V05=(find(diff((mginf(a,:)<0.5))));
    dX=mginf(a,V05+[-5:5]);
    dV=Vamp(V05+[-5:5]);
    flXinf(a).f=fit(dV',dX','poly1');
    slopeXinf(a)=flXinf(a).f.p1;
   % hold on
   % plot(Vamp,flXinf(a).f(Vamp),'k')
end
subplot(4,4,6)
h=plot([0.2:0.2:1],slopeXinf,'o--k');
formatFig(h,1)


%Now calculate the value of Tau_N. Here we fit a slow and a fast time
%constant. When the fractional order is 1 then both time constant are the
%same and equal to the classic Hodgkin-Huxley model
for cc1=1:length(gAll);
    tg=gAll(cc1).g;
    for b=1:length(tg)
        tauV=sort(1./[tg(b).fL.tau tg(b).fL.tau2 ]);
        tauL(cc1,b)=tauV(1);
        tauE(cc1,b)=tauV(2);
        mginf(cc1,b)=tg(b).nfinf;
    end
end

subplot(4,4,10)
cla
v2p=Vamp([1:end] );
tau2p=tauL(1:end,[1:end])';
tau2E=tauE(:,[1:end])';
h=plot(v2p,tau2p);
formatFig(h,1);
 hold on
 h=plot(v2p,tau2E,'--');
 formatFig(h,1);
axis([-110 120 0 35]);

print -depsc Frac_N_Analysis.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 5 -analysis of power-law dynamics H gate
clear
homedir='.';%type your directory
cd(homedir)

load fractoinalHgateSweep

f=fittype('ninf*(1-exp(-x*tau))+n2*(1-exp(-x*tau2))');
f0=fitoptions(f);
g=[];
cc1=1;
fs={'02','04','06','08','10'};
V2u=[1:length(Vamp)];

for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    cc2=1;
    for b=V2u
        g(cc2).tr=(o2a(b).t(logical(o2a(b).t>=30))-30)';
        d2f=o2a(b).nV(logical(o2a(b).t>=30));
        g(cc2).nfinf=d2f(end);
        g(cc2).d2f=d2f-d2f(1);
        g(cc2).d2fraw=d2f;
        if cc1==5
            f0.Start=[(g(cc2).d2f(end))*[1 1] 1 1 ];
        else
            f0.Start=[(g(cc2).d2f(end))*[0.1 0.5 ] 0.002 1 ];
        end
        
        if (max(g(cc2).d2f))>0 && min(g(cc2).d2f)==0
            f0.Upper=[abs(g(cc2).d2f(end))*[1.5 1.5 ] 40 40 ];
            f0.Lower=[0*[1.5 1.5 ] 0 0];
        elseif (max(g(cc2).d2f))>0 && min(g(cc2).d2f)<0
            f0.Upper=[0*[1.5 1.5 ] 40 40 ];
            f0.Lower=[min(g(cc2).d2f)*[1.5 1.5 ] 0 0];
        elseif (max(g(cc2).d2f))==0
            f0.Upper=[0*[1.5 1.5 ] 40 40 ];
            f0.Lower=[-abs(g(cc2).d2f(end))*[1.5 1.5 ] 0 0];
        else 
            disp('WHAaaaat')
        end

        
        g(cc2).fL=fit(g(cc2).tr,...
            g(cc2).d2f,f,f0);
        tauVdummy=1./[g(cc2).fL.tau g(cc2).fL.tau2];
        ref2E=find(g(cc2).tr>min(tauVdummy));
        if isempty(ref2E)
            ref2E=length(g(cc2).tr);
        end
        clf
        plot(g(cc2).tr,(g(cc2).d2f),'k',...
            g(cc2).tr(ref2E(1)),g(cc2).d2f(ref2E(1)),'r*',...
            g(cc2).tr,g(cc2).fL(g(cc2).tr),'r');
        drawnow
%        input('r')
        
        cc2=cc2+1;
    end
    gAll(cc1).g=g;
end

clf
%Plot the value of the H gate vs time for V=-70  and values of fractional 
%order 0.2, 0.4, 0.6, 0.8, 1.0
Vamp(15)
subplot(4,4,1)
cla
for a=1:length(gAll)
    trEx4(:,a)=gAll(a).g(15).tr;
    d2a4r(:,a)=gAll(a).g(15).d2fraw;
    d2a4(:,a)=gAll(a).g(15).d2f;
    fL(a).f=gAll(a).g(15).fL;
    trEx12(:,a)=gAll(a).g(15).tr;
    d2a12(:,a)=gAll(a).g(15).d2fraw;
end
h=plot(trEx4(:,1),(d2a4r),'-');
formatFig(h,1)
axis([0 100 0.6 0.8]);

subplot(4,4,5)
cla
%Log-log plot of the same data and fitting - Not used in the paper
%Same data as in subplot(4,4,1)
cv='bgrmc';
for a=1:size(d2a4,2)
    logfit(a).f=fit(log10(trEx4(2:500,1)),log10(d2a4(2:500,a)),'poly1');
    
    h=plot(log10(trEx4(:,1)),log10(d2a4(:,a)),cv(a));
    hold on
    fracexpfig(a)=logfit(a).f.p1;
end
axis([-2 2 -2 -0.5])

%do an adaptive log fit for all the traces! Not used in the paper
for a=1:length(gAll)
    for b=1:length(gAll(a).g)
        g2f=gAll(a).g(b);
        out=adaptiveLogFit(g2f,0.95);
        ftM(a,b)=out.time;% get the average from this.
    end
end

%plot the exponent of the log fit to the voltage response - Not used in the
%paper
subplot(4,4,9)
h=plot(fracexpfig,'o--k');
formatFig(h,1)

%Plot the value of H_inf. In the paper we call this the long term response
%fuction
subplot(4,4,2)
cla

for cc1=1:length(gAll);
    tg=gAll(cc1).g;
    for b=1:length(tg)
        mginf(cc1,b)=tg(b).nfinf;
    end
end

h=plot(Vamp(V2u),mginf','-');
formatFig(h,1);
axis([-110 20 0 1])

%now calculate the slope at the half max value - Not used in the paper
for a=1:size(mginf,1)   
    V05=(find(diff((mginf(a,:)<0.5))));
    dX=mginf(a,V05+[-5:5]);
    dV=Vamp(V05+[-5:5]);
    flXinf(a).f=fit(dV',dX','poly1');
    slopeXinf(a)=flXinf(a).f.p1;
end
subplot(4,4,6)
h=plot([0.2:0.2:1],slopeXinf,'o--k');
formatFig(h,1)

%Now calculate the value of Tau_H. Here we fit a slow and a fast time
%constant. When the fractional order is 1 then both time constant are the
%same and equal to the classic Hodgkin-Huxley model
%h tau
for cc1=1:length(gAll);
    tg=gAll(cc1).g;
    for b=1:length(tg)
        tauV=sort(1./[tg(b).fL.tau tg(b).fL.tau2 ]);
        tauL(cc1,b)=tauV(1);
        tauE(cc1,b)=tauV(2);
        mginf(cc1,b)=tg(b).nfinf;
    end
end

subplot(4,4,10)
cla
v2p=Vamp([1:end] );
tau2p=tauL(1:end,[1:end])';
tau2E=tauE(:,[1:end])';
h=plot(v2p,tau2p);
formatFig(h,1);
 hold on
 h=plot(v2p,tau2E,'--');
 formatFig(h,1);
axis([-110 20 0 50]);

print -depsc Frac_H_Analysis.eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 6 - Analysis of power-law dynamics M gate

clear
load fractoinalMgateSweep

f=fittype('ninf*(1-exp(-x*tau))'); % here we used a single time constant

f0=fitoptions(f);
g=[];
cc1=1;
fs={'02','04','06','08','10'};
V2u=[1:length(Vamp)];

for cc1=1:length(fs);
    eval(['o2a=out' fs{cc1} ';']);
    cc2=1;
    for b=V2u
        if sum(isnan(o2a(b).nV))==0
            g(cc2).tr=(o2a(b).t(logical(o2a(b).t>=20))-20)';
            d2f=o2a(b).nV(logical(o2a(b).t>=20));
            g(cc2).nfinf=d2f(end);
            g(cc2).d2f=d2f-d2f(1);
            g(cc2).d2fraw=d2f;
            if cc1==5
                f0.Start=[(g(cc2).d2f(end))*[1] 1  ];
            else
                f0.Start=[(g(cc2).d2f(end))*[0.1 ] 0.2  ];
            end
            
            if sum(g(cc2).d2f)==0
                f0.Upper = [1  40];
                f0.Lower = [-1  0 ];
            elseif (max(g(cc2).d2f))>0 && min(g(cc2).d2f)==0
                f0.Upper=[abs(g(cc2).d2f(end))*[ 1.5 ]  40 ];
                f0.Lower=[0*[1.5 ]  0];
            elseif (max(g(cc2).d2f))>0 && min(g(cc2).d2f)<0
                f0.Upper=[0*[ 1.5 ]  40 ];
                f0.Lower=[min(g(cc2).d2f)*[ 1.5 ]  0];
            elseif (max(g(cc2).d2f))==0
                f0.Upper=[0*[ 1.5 ]  40 ];
                f0.Lower=[-abs(g(cc2).d2f(end))*[ 1.5 ]  0];
            else
                disp('WHAaaaat')
            end
            
            if sum(f0.Upper(1:2)==f0.Lower(1:2))==2
                f0.Upper(1:2)=0.2;
                f0.Lower(1:2)=-0.2;
               
            end
            
            g(cc2).fL=fit(g(cc2).tr,...
                g(cc2).d2f,f,f0);
            tauVdummy=1./[g(cc2).fL.tau ];%g(cc2).fL.tau2];
            ref2E=find(g(cc2).tr>min(tauVdummy));
            if isempty(ref2E)
                ref2E=length(g(cc2).tr);
            end
            [b Vamp(b)]
            clf
            plot(g(cc2).tr,(g(cc2).d2f),'k',...
                g(cc2).tr(ref2E(1)),g(cc2).d2f(ref2E(1)),'r*',...
                g(cc2).tr,g(cc2).fL(g(cc2).tr),'r');
            %        input('r') %just in case you can to see the results
        else 
            g(cc2).tr=[];
            g(cc2).nfinf=[];
            g(cc2).d2f=[];
            g(cc2).d2fraw=[];
            g(cc2).fL=[];
        end 
        cc2=cc2+1;
    end
    gAll(cc1).g=g;
end

clf

%Plot the value of the M gate vs time for V=-55  and values of fractional 
%order  0.4, 0.6, 0.8, 1.0

subplot(4,4,1)
cla
for a=1:length(gAll)
   if ~isempty(gAll(a).g(12).tr)
    trEx4(:,a)=gAll(a).g(12).tr;
    d2a4r(:,a)=gAll(a).g(12).d2fraw;
    d2a4(:,a)=gAll(a).g(12).d2f;
    fL(a).f=gAll(a).g(12).fL;
    trEx12(:,a)=gAll(a).g(12).tr;
    d2a12(:,a)=gAll(a).g(12).d2fraw;
   end
end
h=plot(trEx4(:,1),(d2a4r),'-');
formatFig(h,1)
axis([0 20 0 0.2]);

subplot(4,4,5)
cla
h=loglog(trEx4(:,1),(d2a4));
formatFig(h,1)
cla

%Log-log plot of the same data and fitting - Not used in the paper
%Same data as in subplot(4,4,1)
cv='bgrmc';
for a=1:size(d2a4,2)
    logfit(a).f=fit(log10(trEx4(2:500,1)),log10(d2a4(2:500,a)),'poly1');
    
    h=plot(log10(trEx4(:,1)),log10(d2a4(:,a)),cv(a));
    hold on
   % h=plot(log10(trEx4(2:end,1)),(logfit(a).f(log10(trEx4(2:end,1)))),cv(a));
    fracexpfig(a)=logfit(a).f.p1;
end
formatFig(h,1)
axis([-2 1 -1.5 -0.8])

%plot the exponent of the log fit to the voltage response - Not used in the
%paper.
for a=1:length(gAll)
    for b=1:length(gAll(a).g)
        g2f=gAll(a).g(b);
        out=adaptiveLogFit(g2f,0.95);
        ftM(a,b)=out.time;% get the average from this.
    end
end

%plot the exponent of the log fit to the voltage response
subplot(4,4,9)
h=plot(fracexpfig,'o--k');
formatFig(h,1)

%Plot the value of M_inf. In the paper we call this the long term response
%fuction
subplot(4,4,2)
cla
for cc1=1:length(gAll);
    tg=gAll(cc1).g;
    for b=1:length(tg)
        if ~isempty(tg(b).nfinf)
          mginf(cc1,b)=tg(b).nfinf;
        else
            mginf(cc1,b)=0;
        end
    end
end
%some fits got confused because the values were aroung zero, so we deleted them
h=plot(Vamp([1:14 16:end]),mginf(:,[1:14 16:end])','-');
formatFig(h,1);
axis([-110 60 0 1])

%now calculate the slope at the half max value - Not used in the paper

for a=1:size(mginf,1)
    
    V05=(find(diff((mginf(a,:)<0.5))));
    dX=mginf(a,V05(1)+[-5:5]);
    dV=Vamp(V05(1)+[-5:5]);
    flXinf(a).f=fit(dV',dX','poly1');
    slopeXinf(a)=flXinf(a).f.p1;
   % hold on
   % plot(Vamp,flXinf(a).f(Vamp),'k')
end
subplot(4,4,6)
h=plot([0.2:0.2:1],slopeXinf,'o--k');
formatFig(h,1)
axis([0 1 0 0.04])

%Now calculate the value of Tau_M. Here we fit a slow and a fast time
%constant. When the fractional order is 1 then both time constant are the
%same and equal to the classic Hodgkin-Huxley model.
for cc1=1:length(gAll);
    tg=gAll(cc1).g;
    for b=1:length(tg)
        if ~isempty(tg(b).nfinf)
        tauV=sort(1./[tg(b).fL.tau]);% tg(b).fL.tau2 ]);
        tauL(cc1,b)=tauV(1);
        tauE(cc1,b)=tauV(1);
        else
            tauL(cc1,b)=0;
            tauE(cc1,b)=0;
        end
        
    end
end

subplot(4,4,10)
cla
v2p=Vamp([1:9 12:14 16:end] );
tau2p=tauL(1:end,[1:9 12:14 16:end])';
tau2E=tauE(:,[1:9 12:14 16:end])';
h=plot(v2p,tau2p);
formatFig(h,1);
hold on
h=plot(v2p,tau2E,'-');
formatFig(h,1);
axis([-110 120 0 0.6]);


print -depsc Frac_M_Analysis.eps

