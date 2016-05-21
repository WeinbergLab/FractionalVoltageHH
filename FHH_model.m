function FHH_model(a, Iapp, T, dt, dn, dur, bcl, save_flag, toutput)

if nargin<1
    a = .8;      % fractional-order
    Iapp = 20;  % applied stimulus amplitude
    T = 1e3;     % total time, in ms
    dt = 1e-2;  % time step, in ms
    dn = 10;    % sampling interval
    dur = T;   % stimulus duration, in ms 
    bcl = T;  % stimulus period, in ms
    save_flag = 1;
    toutput = 100;   % interval to output data, if saved

end

% tags for plotting and saving
ls = '-';  % line style
fs = 24;  % font size
plot_flag = 1;


% 0 - ode15s (a = 1)
% 1 - euler (a = 1)
% 2 - GL-euler (0 < a <= 1)

if save_flag
    savename = ['FHH_a',num2str(round(100*a)),'_Iapp',num2str(Iapp)];
else
    savename = [];
end

int_flag = 0;  
if a ~= 1
    int_flag = 2;
end

% initial conditions
V0 = 0; m0 = .0529; h0 = .5961; n0 = .3177;
X0 = [V0; m0; h0; n0];


tonset = 0; 
tic;
switch int_flag
    case 0
         opt = odeset('maxstep',.05);
        [t,X] = ode15s(@(t,x) HHode(t,x,tonset,dur,Iapp,bcl), 0:dt:T, X0,opt);
    case 1
        [t,X] = euler(@(t,x) HHode(t,x,tonset, dur,Iapp,bcl), dt, T, X0);
    case 2
        [t,X] =  GL_euler(@(t,x) HHode(t,x,tonset, dur,Iapp,bcl), dt, a, T, X0, toutput, savename, dn);
end
toc;

% downsample
t = t(1:dn:end);
X = X(1:dn:end,:);

V = X(:,1); m = X(:,2); h = X(:,3); n = X(:,4);
[INa, IK, IL, Iext] = computeCurrents(t, V, m, n, h, tonset, dur, Iapp, bcl);
if save_flag
    save(savename,'t','X','a','Iapp','dur','bcl','Iext');
end

if plot_flag
    s = 'k';  Vrest = -80;
    Vm = V + Vrest;
    subplot(3,1,1);
    plot(t,Vm,[s,ls],'linewidth',2); hold on;
    set(gca,'fontsize',fs,'ylim',[min(Vm)-5 max(Vm)+5]);
    ylabel('V (mV)');
    subplot(3,1,2); 
    plot(t,m,['r',ls],t,h,['b',ls],t,n,['g',ls],'linewidth',2); hold on;
    legend('m','h','n');
    set(gca,'fontsize',fs);
    subplot(3,1,3);
    plot(t, INa,['b',ls],t, IK,['r',ls],t,IL,['g',ls],'linewidth',2); hold on;
    legend('INa','IK','IL');
    set(gca,'fontsize',fs);
end

function [t,X] =  GL_euler(f, h, a, T, X0,toutput, savename, dn)
t = (0:h:T)';
X = nan(length(t),length(X0));
X(1,:) = X0;

c = nan(1,length(t)); c(1) = a;
for j = 2:length(t)
    c(j) = (1-(1+a)/j)*c(j-1);
end

for i = 1:length(t)-1
    dX = f(t(i),X(i,:));
    X(i+1,1) = h^a*dX(1) + memo(X(:,1),c,i+1)+ X0(1)*h^a*t(i+1)^(-a)/gamma(1-a); 
    X(i+1,2:4) = h*dX(2:4) + X(i,2:4)';
    
    if ~mod(t(i),toutput) && ~isempty(savename)
        tsim = t(1:dn:i);
        Xsim = X(1:dn:i,:);
        save(savename,'tsim','Xsim','a');
        
        [a t(i)]
    end
end


function yo = memo(x,c,k)
% iend = min(k-1,1e4);
%yo = sum(c(1:iend)'.*x(k-1:-1:k-1-iend+1));
 yo = sum(c(1:k-1)'.*x(k-1:-1:1));  % more efficient



function [t,X] = euler(f, h, T, X0)
t = (0:h:T)';
X = nan(length(t), length(X0));
X(1,:) = X0;
for i = 1:length(t)-1
    X(i+1,:) = X(i,:)' + h*f(t(i),X(i,:));
end


function dX = HHode(t,X, tonset, dur, amp, bcl)
V = X(1);
m = X(2);
h = X(3);
n = X(4);

Iext = ((mod(t,bcl)>=tonset)&&mod(t,bcl)<tonset+dur)*amp;

% parameters
gL = .3;        % mS/cm^2
gK = 36;
gNa = 120;
C = 1;          % uF/cm^2
EL = 10.6;       % mV
EK = -12;
ENa = 115;

Iion = gNa*m^3*h*(V-ENa) + gK*n^4*(V-EK) + gL*(V-EL);
dV = (Iext-Iion)/C;

dm = am(V)*(1-m)-bm(V)*m;
dh = ah(V)*(1-h)-bh(V)*h;
dn = an(V)*(1-n)-bn(V)*n;

dX = zeros(4,1);
dX(1) = dV;
dX(2) = dm;
dX(3) = dh;
dX(4) = dn;

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

function [INa, IK, IL, Iext] = computeCurrents(t, V, m, n, h, t0, dur, Iapp, bcl)
% parameters
gL = .3;        % mS/cm^2
gK = 36;
gNa = 120;
EL = 10.6;       % mV
EK = -12;
ENa = 115;

Iext = ((mod(t,bcl)>=t0)&mod(t,bcl)<t0+dur)*Iapp;
INa = gNa*m.^3.*h.*(V-ENa);
IK =  gK*n.^4.*(V-EK);
IL = gL*(V-EL);

return




