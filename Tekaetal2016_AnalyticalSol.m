% Simulation of power-law dynamic gate variables in the Hodgkin-Huxley
% model using fractional order derivatives.
% Teka W, Stockton D, Santamaria F. "Power-law dynamics of membrane
% conductances increase spiking diversity in a Hdgkin-Huxley model" PLoS
% Computational Biology, in press, 2016.
% If you use this software please reference our paper. 

%This script computes the analytical solution to the the gate variable
%equation with constant voltage input. The results are used to validate the
%numerical integration of the fractional order derivative used for the
%paper. See Figure 1 of the paper.
%% Section 13
% Compute the analytical solution of power-law N gate and compare to 
%numerical solution
clear
load fractoinalNgateSweep
t=0:1e-3:80;
v0=-65;

c1=1;
for V=Vamp
    alphan=(0.1-0.01*(V-v0))./(exp(1-0.1*(V-v0))-1);
    betan=0.125.*exp(-(V-v0)./80);
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    x0=0.0529;
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        xt(c,:)=xinf+(0.3177-xinf)*mlf(eta,1,-(t.^eta)/taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=30)&(tt<=110)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=30)&(tt<=110));
    n04=out04(t2a).nV((tt>=30)&(tt<=110));
    n06=out06(t2a).nV((tt>=30)&(tt<=110));
    n08=out08(t2a).nV((tt>=30)&(tt<=110));
    n10=out10(t2a).nV((tt>=30)&(tt<=110));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    
    n_inf(:,c1)=xt(:,end);
    nc_inf(:,c1)=[n02(end) n04(end) n06(end) n08(end) n10(end)]';
    err(c1)=(size(xt,2)*5)^-1*sqrt(sum(...
        (n02-xt(1,:)').^2+...
        (n04-xt(2,:)').^2+(n06-xt(3,:)').^2+...
        (n08-xt(4,:)').^2+(n10-xt(5,:)').^2));
    c1=c1+1;
    clf
    plot(t,xt','.k')
    xlim([0 2])
    hold on
    
    %t2a=12;
    plot(ntr,n02,'r',ntr,n04,'g',...
        ntr,n06,'.b',ntr,n08,'.m',...
        ntr,n10,'.y')

    drawnow 
    V
%     input('r')
end


%% Section 14
% Compute the analytical solution of power-law M gate and compare to 
%numerical solution

clear
load fractoinalMgateSweep
t=0:5e-4:20;
v0=-65;

c1=1;
for V=Vamp
    
    %V=-55;
    alphan=((2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1));
    betan=(4.*exp(-(V-v0)./18));
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    x0=0.0529;
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        %xt(c,:)=xinf+(0.0529-xinf)*MittagLeffler(-(t./taux).^eta,eta,1);
        xt(c,:)=xinf+(0.0529-xinf)*mlf(eta,1,-(t.^eta)/taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=20)&(tt<=40)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=20)&(tt<=40));
    n04=out04(t2a).nV((tt>=20)&(tt<=40));
    n06=out06(t2a).nV((tt>=20)&(tt<=40));
    n08=out08(t2a).nV((tt>=20)&(tt<=40));
    n10=out10(t2a).nV((tt>=20)&(tt<=40));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    
    err(c1)=(size(xt,2)*5)^-1*sqrt(sum(...
        (n02-xt(1,:)').^2+...
        (n04-xt(2,:)').^2+(n06-xt(3,:)').^2+...
        (n08-xt(4,:)').^2+(n10-xt(5,:)').^2));
    c1=c1+1;
    clf
    plot(t,xt','.k')
    xlim([0 2])
    hold on
    
    %t2a=12;
    plot(ntr,n02,'r',ntr,n04,'g',...
        ntr,n06,'.b',ntr,n08,'.m',...
        ntr,n10,'.k')

    drawnow 
    V
%     input('r')
end




%% Section 15
% Compute the analytical solution of power-law H gate and compare to 
%numerical solution

clear
load fractoinalHgateSweep


t=0:1e-3:20;
v0=-65;
V=00;
Vzero=0;

c1=1;
%this for loop goes throug the same V values used for the numerical simulations
for V=Vamp
    alphan=(0.07*exp(-(V-v0)/20));
    betan=(1./(exp(3-0.1*(V-v0))+1));
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        %this is the solution for a given eta and time.
        xt(c,:)=xinf+(0.596-xinf)*mlf(eta,1,-(t.^eta)/taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=30)&(tt<=50)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=30)&(tt<=50));
    n04=out04(t2a).nV((tt>=30)&(tt<=50));
    n06=out06(t2a).nV((tt>=30)&(tt<=50));
    n08=out08(t2a).nV((tt>=30)&(tt<=50));
    n10=out10(t2a).nV((tt>=30)&(tt<=50));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    
    err(c1)=(size(xt,2)*5)^-1*sqrt(sum(...
        (n02-xt(1,:)').^2+...
        (n04-xt(2,:)').^2+(n06-xt(3,:)').^2+...
        (n08-xt(4,:)').^2+(n10-xt(5,:)').^2));
    c1=c1+1;
    clf
    plot(t,xt','.k')
    xlim([0 2])
    hold on
    
    %t2a=12;
    plot(ntr,n02,'r',ntr,n04,'g',...
        ntr,n06,'.b',ntr,n08,'.m',...
        ntr,n10,'.k')

    drawnow 
    V
%     input('r')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% samples 
clear
load fractoinalNgateSweep
t=0:1e-3:80;
v0=-65;

c1=1;
for V=Vamp(15)
    %V=-55;
    alphan=(0.1-0.01*(V-v0))./(exp(1-0.1*(V-v0))-1);
    betan=0.125.*exp(-(V-v0)./80);
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    x0=0.0529;
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        xt(c,:)=xinf+(0.3177-xinf)*mlf(eta,1,-(t.^eta)./taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=30)&(tt<=110)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=30)&(tt<=110));
    n04=out04(t2a).nV((tt>=30)&(tt<=110));
    n06=out06(t2a).nV((tt>=30)&(tt<=110));
    n08=out08(t2a).nV((tt>=30)&(tt<=110));
    n10=out10(t2a).nV((tt>=30)&(tt<=110));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;        
    c1=c1+1;    
end
clf
subplot(4,4,1)
h=plot(t,xt');
formatFig(h,1)
hold on
r=1:1500:length(ntr);
h=plot(ntr(r),n02(r),'.b',ntr(r),n04(r),'.g',...
    ntr(r),n06(r),'.r',ntr(r),n08(r),'.c',...
    ntr(r),n10(r),'.m');
formatFig(h,1)
axis([0 30 0.3177 1])

clear
load fractoinalMgateSweep
t=0:5e-4:40;
v0=-65;

c1=1;
for V=Vamp(25)
    alphan=((2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1));
    betan=(4.*exp(-(V-v0)./18));
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    x0=0.0529;
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        xt(c,:)=xinf+(0.0529-xinf)*mlf(eta,1,-(t.^eta)/taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=20)&(tt<=60)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=20)&(tt<=60));
    n04=out04(t2a).nV((tt>=20)&(tt<=60));
    n06=out06(t2a).nV((tt>=20)&(tt<=60));
    n08=out08(t2a).nV((tt>=20)&(tt<=60));
    n10=out10(t2a).nV((tt>=20)&(tt<=60));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    c1=c1+1;
end
subplot(4,4,2)
h=plot(t,xt');
formatFig(h,1)
hold on
r=1:1500:length(ntr);
h=plot(ntr(r),n02(r),'.b',ntr(r),n04(r),'.g',...
    ntr(r),n06(r),'.r',ntr(r),n08(r),'.c',...
    ntr(r),n10(r),'.m');
formatFig(h,1)
axis([0 30 0.529 1])

clear
load fractoinalHgateSweep
t=0:1e-3:100;
v0=-65;
c1=1;
for V=Vamp(15)
    alphan=(0.07*exp(-(V-v0)/20));
    betan=(1./(exp(3-0.1*(V-v0))+1));
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        xt(c,:)=xinf+(0.596-xinf)*mlf(eta,1,-(t.^eta)/taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=30)&(tt<=130)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=30)&(tt<=130));
    n04=out04(t2a).nV((tt>=30)&(tt<=130));
    n06=out06(t2a).nV((tt>=30)&(tt<=130));
    n08=out08(t2a).nV((tt>=30)&(tt<=130));
    n10=out10(t2a).nV((tt>=30)&(tt<=130));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    
    c1=c1+1;
end
subplot(4,4,3)
h=plot(t,xt');
formatFig(h,1)
hold on
r=1:1500:length(ntr);
h=plot(ntr(r),n02(r),'.b',ntr(r),n04(r),'.g',...
    ntr(r),n06(r),'.r',ntr(r),n08(r),'.c',...
    ntr(r),n10(r),'.m');
formatFig(h,1)
axis([0 30 0.596 0.8])
print -depsc2 FracGate_Traces_TheoryNum.eps
save FracGate_Traces_TheoryNum
%% Xinf
clear
load fractoinalNgateSweep
t=0:1e-3:80;
v0=-65;

c1=1;
for V=Vamp
    %V=-55;
    alphan=(0.1-0.01*(V-v0))./(exp(1-0.1*(V-v0))-1);
    betan=0.125.*exp(-(V-v0)./80);
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    x0=0.0529;
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        xt(c,:)=xinf+(0.3177-xinf)*mlf(eta,1,-(t(end).^eta)./taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=30)&(tt<=110)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=30)&(tt<=110));
    n04=out04(t2a).nV((tt>=30)&(tt<=110));
    n06=out06(t2a).nV((tt>=30)&(tt<=110));
    n08=out08(t2a).nV((tt>=30)&(tt<=110));
    n10=out10(t2a).nV((tt>=30)&(tt<=110));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    
    n_inf(:,c1)=xt(:,end);
    nc_inf(:,c1)=[n02(end) n04(end) n06(end) n08(end) n10(end)]';
    c1=c1+1;
end
clf
subplot(4,4,1)
h=plot(Vamp,n_inf);
formatFig(h,1)
hold on
h=plot(Vamp,nc_inf,'.');
formatFig(h,1)
axis([-100 110 0 1])
%print -depsc2 FracN_Ninf_TheoryNum.eps

clear
load fractoinalMgateSweep
t=0:5e-4:40;
v0=-65;

c1=1;
for V=Vamp
    
    %V=-55;
    alphan=((2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1));
    betan=(4.*exp(-(V-v0)./18));
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    x0=0.0529;
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        %xt(c,:)=xinf+(0.0529-xinf)*MittagLeffler(-(t./taux).^eta,eta,1);
        xt(c,:)=xinf+(0.0529-xinf)*mlf(eta,1,-(t(end).^eta)/taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=20)&(tt<=60)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=20)&(tt<=60));
    n04=out04(t2a).nV((tt>=20)&(tt<=60));
    n06=out06(t2a).nV((tt>=20)&(tt<=60));
    n08=out08(t2a).nV((tt>=20)&(tt<=60));
    n10=out10(t2a).nV((tt>=20)&(tt<=60));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    m_inf(:,c1)=xt(:,end);
    mc_inf(:,c1)=[n02(end) n04(end) n06(end) n08(end) n10(end)]';
    c1=c1+1;
end
subplot(4,4,2)
h=plot(Vamp([1:14 16:end]),m_inf(:,[1:14 16:end]));
formatFig(h,1)
hold on
h=plot(Vamp([1:14 16:end]),mc_inf(:,[1:14 16:end]),'.');
formatFig(h,1)
axis([-100 60 0 1])

clear
load fractoinalHgateSweep
t=0:1e-3:100;
v0=-65;
V=00;
Vzero=0;

c1=1;
for V=Vamp
    alphan=(0.07*exp(-(V-v0)/20));
    betan=(1./(exp(3-0.1*(V-v0))+1));
    xinf=alphan./(alphan+betan);
    taux=1./(alphan+betan);
    xt=[];
    c=1;
    for eta=0.2:0.2:1;
        xt(c,:)=xinf+(0.596-xinf)*mlf(eta,1,-(t(end).^eta)/taux,6);
        c=c+1;
    end
    t2a=find(Vamp==V);
    tt=out04(t2a).t;
    nt=tt(logical((tt>=30)&(tt<=130)));
    ntr=nt-nt(1);
    n02=out02(t2a).nV((tt>=30)&(tt<=130));
    n04=out04(t2a).nV((tt>=30)&(tt<=130));
    n06=out06(t2a).nV((tt>=30)&(tt<=130));
    n08=out08(t2a).nV((tt>=30)&(tt<=130));
    n10=out10(t2a).nV((tt>=30)&(tt<=130));
    
    n02(logical(~((n02<=1)&(n02>=0))))=0;
    n04(logical(~((n04<=1)&(n04>=0))))=0;
    n06(logical(~((n06<=1)&(n06>=0))))=0;
    n08(logical(~((n08<=1)&(n08>=0))))=0;
    n10(logical(~((n10<=1)&(n10>=0))))=0;
    h_inf(:,c1)=xt(:,end);
    hc_inf(:,c1)=[n02(end) n04(end) n06(end) n08(end) n10(end)]';
    
    c1=c1+1;

end
subplot(4,4,3)
h=plot(Vamp,h_inf);
formatFig(h,1)
hold on
h=plot(Vamp,hc_inf,'.');
formatFig(h,1)
axis([-110 10 0 1])
print -depsc2 FracGate_Xinf_TheoryNum.eps