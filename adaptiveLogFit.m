function out=adaptiveLogFit(s,r2)

t=s.tr;
d=s.d2fraw;

if isempty(s.tr)
    out.f=[];
    out.g=0;
    out.time=0;
    return
end

tmax=max(t);
epV=[0.5 1 5 10 20 40 60 80 100 200 300 500 600 700 800 900 1000];
epE=epV(logical(epV<=tmax));
for a=1:length(epE); 
   p2f=2:find(diff([t<=epE(a); 0]));
   [f, g]=fit(log10(t(p2f,1)),log10(d(p2f)),'poly1');
   r2v(a)=g.rsquare;
   fs(a).f=f;
   fs(a).g=g;
end

r2test=1;
a=1;
while r2test && a<=length(epE)
    r2test=~(r2v(a)<r2);
    thisP=a;
    a=a+1;
end
out=fs(thisP);
out.time=epE(thisP);