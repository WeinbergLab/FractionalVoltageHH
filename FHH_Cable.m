function FHH_Cable(a, Iapp, T, dt, dn, dur, rhoi, L, dx, save_flag)
if nargin<1
    a = .5;      % fractional-order
    Iapp = 100;  % applied stimulus amplitude
    T = 40;     % total time, in ms
    dt = 1e-3;  % time step, in ms
    dn = 1;    % sampling interval
    dur = .1;   % stimulus duration, in ms  
    rhoi = 0.0354; % k-ohm * cm   % resistivity
    L = 1;  % cm (cable length)
    dx = 0.05;  % cm (discretization)
    save_flag = 1;
end

ls = '--';  % line style
fs = 24;  % font size
plot_flag = 0;

% parameters 
d = 10e-4;  % cm (cable diameter)
% cable conductance g = d/(4*rhoi) % in mS
Nodes = round(L/dx)+1;

% conduction velocity measurement locations
m1 = .25*L; m2 = .75*L;     % cm

% 0 - ode15s (a = 1)
% 1 - euler (a = 1)
% 2 - GL-euler (0 < a <= 1)


int_flag = 0;  
if a ~= 1
    int_flag = 2;
end

V0 = 0; m0 = .0529; h0 = .5961; n0 = .3177;
X0 = [V0*ones(Nodes,1); m0*ones(Nodes,1); h0*ones(Nodes,1); n0*ones(Nodes,1)];

fNa = 1; fK = 1;
tonset = 0;  

tic;
switch int_flag
    case 0
         opt = odeset('maxstep',.05);
        [t,X] = ode15s(@(t,x) HHode(t,x,tonset,dur,Iapp,Nodes,d,rhoi,dx), 0:dt:T, X0,opt);
    case 1
        [t,X] = euler(@(t,x) HHode(t,x,tonset, dur,Iapp,Nodes,d,rhoi,dx), dt, T, X0);
    case 2
        [t,X] =  GL_euler(@(t,x) HHode(t,x,tonset, dur,Iapp, Nodes, d, rhoi,dx,fNa, fK), dt, a, T, X0,Nodes);
end
toc;

% downsample
t = t(1:dn:end);
X = X(1:dn:end,:);
x = 0:dx:L;

[~,k1] = min(abs(x-m1));
[~,k2] = min(abs(x-m2));

V = X(:,1:Nodes); m = X(:,Nodes+1:2*Nodes); 
h = X(:,2*Nodes+1:3*Nodes); n = X(:,3*Nodes+1:4*Nodes);
[INa, IK, IL, Iext] = computeCurrents(t, V, m, n, h, tonset, dur, Iapp, fNa, fK);

V1 = V(:,k1); V2 = V(:,k2);

% [ISI1, t1] = measISI(t,V1);
% [ISI2, t2] = measISI(t,V2);

t1 = t1(1:length(t2));
cv = (m2-m1)./(t2-t1);  % cm/ms

if save_flag
    savename = ['FHH_Cable_a',num2str(round(100*a)),'_Iapp',num2str(Iapp)];
    save(savename,'t','X','x');
end

if plot_flag
    figure(17);
    s = 'k';  Vrest = -80;
    Vm = V + Vrest;
    subplot(3,1,1);
    plot(t,Vm,[s,ls],'linewidth',2); hold on;
    set(gca,'fontsize',fs,'ylim',[min(Vm(:))-5 max(Vm(:))+5]);
    ylabel('V (mV)');
    subplot(3,1,2); 
    plot(t,m,['r',ls],t,h,['b',ls],t,n,['g',ls],'linewidth',2); hold on;
    legend('m','h','n');
    set(gca,'fontsize',fs);
    subplot(3,1,3);
    plot(t, INa,['b',ls],t, IK,['r',ls],t,IL,['g',ls],'linewidth',2); hold on;
    legend('INa','IK','IL');
    set(gca,'fontsize',fs);
    
    figure(18); 
    imagesc(t,x,Vm'); colorbar; set(gca,'fontsize',18);
    xlabel('time (ms)'); ylabel('x (cm)');
end

function [t,X] =  GL_euler(f, h, a, T, X0, Na)
t = (0:h:T)';
X = nan(length(t),length(X0));
X(1,:) = X0;

c = nan(1,length(t)); c(1) = a;
for j = 2:length(t)
    c(j) = (1-(1+a)/j)*c(j-1);
end

for i = 1:length(t)-1
    dX = f(t(i),X(i,:));
    X(i+1,1:Na) = h^a*dX(1:Na) + (memo(X(:,1:Na),c,i+1))'+ X0(1:Na)*h^a*t(i+1)^(-a)/gamma(1-a); 
    X(i+1,Na+1:end) = h*dX(Na+1:end) + X(i,Na+1:end)';
   
end

function yo = memo(x,c,k)
 yo = (c(1:k-1)*x(k-1:-1:1,:));



function [t,X] = euler(f, h, T, X0)
t = (0:h:T)';
X = nan(length(t), length(X0));
X(1,:) = X0;
for i = 1:length(t)-1
    X(i+1,:) = X(i,:)' + h*f(t(i),X(i,:));
end


function dX = HHode(t,X, tonset, dur, amp,Nodes,d,rhoi,dx, fNa, fK)
V = X(1:Nodes);
m = X(Nodes+1:2*Nodes);
h = X(2*Nodes+1:3*Nodes);
n = X(3*Nodes+1:4*Nodes);

stim_sites = zeros(Nodes,1); 
stim_sites(1) = 1;  % stimulus at node 1
Iext = ((t>=tonset)&&t<tonset+dur)*amp*stim_sites;

% parameters
gL = .3;        % mS/cm^2
gK = 36*fK;
gNa = 120*fNa;
C = 1;          % uF/cm^2
EL = 10.6;       % mV
EK = -12;
ENa = 115;

Iion = gNa*m.^3.*h.*(V-ENa) + gK*n.^4.*(V-EK) + gL*(V-EL);

Vvec = zeros(Nodes+2,1);
Vvec(2:end-1) = V;
Vvec(1) = V(1); Vvec(end) = V(end);

dV = (Iext-Iion(:) + d/(4*rhoi)*(Vvec(1:end-2)-2*Vvec(2:end-1)+Vvec(3:end))/dx^2)/C;

dm = am(V).*(1-m)-bm(V).*m;
dh = ah(V).*(1-h)-bh(V).*h;
dn = an(V).*(1-n)-bn(V).*n;

dX = zeros(4*Nodes,1);
dX(1:Nodes) = dV;
dX(1+Nodes:2*Nodes) = dm;
dX(1+2*Nodes:3*Nodes) = dh;
dX(1+3*Nodes:4*Nodes) = dn;

return

function am = am(V)
am = 0.1*(-V+25)./(exp(-0.1*V+2.5)-1);
return

function bm = bm(V)
bm = 4*exp(-V/18);
return

function ah = ah(V)
ah = .07*exp(-V/20);
return

function bh = bh(V)
bh = 1./(exp(-0.1*V+3)+1);
return

function an = an(V)
an = 0.01*(-V+10)./(exp(-0.1*V+1)-1);
return

function bn = bn(V)
bn = .125*exp(-V/80);
return

function [INa, IK, IL, Iext] = computeCurrents(t, V, m, n, h, t0, dur, Iapp, fNa, fK)
% parameters
gL = .3;        % mS/cm^2
gK = 36*fK;
gNa = 120*fNa;
EL = 10.6;       % mV
EK = -12;
ENa = 115;

Iext = ((t>=t0)&t<t0+dur)*Iapp;
INa = gNa*m.^3.*h.*(V-ENa);
IK =  gK*n.^4.*(V-EK);
IL = gL*(V-EL);

return

function [ISI, peaks] = measISI(t,x)

peaks = [];
for i = 2:length(t)-1
   if x(i)>x(i-1) && x(i)>x(i+1)
       peaks = [peaks t(i)];
   end
end

if length(peaks)>1
    ISI = diff(peaks);
else
    ISI = nan;
end
return




