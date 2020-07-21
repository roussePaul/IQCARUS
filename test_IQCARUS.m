clear all
clc

A = - eye(2);
B = [1;0];
Bu = [1;1];

u = @(t) sin(t);

M = diag([1,1,0,-1]);

E0 = eye(2); E0 = E0*E0';
f0 = [1;2];
g0 = [1];


obj = IQCARUS(A,B,Bu,M, 0, E0,f0,g0);

T = [0,0.3];

obj.sim_until(0.1,u);

figure(1)
clf
hold on

x = E0^-1*f0;

[tt,xx,xq] = obj.pts_section_plane(0.1,100,[1;0],x);
plot(tt,xq)



%%
Mxw = randn(5,2);
B = randn(5,2);
ll = eig(expm([-B;Mxw]*[Mxw',B']))
norm([expm(A)])

%%
n = 2;
Mw = -1; Mwi = Mw^-1;
Mxw = zeros(n,1);
Mx = eye(n);
A = -eye(n);
B = ones(n,1);
%%
n = 2;
Mw = -1; Mwi = Mw^-1;
Mxw = zeros(n,1);
Mx = randn(n); Mx = Mx*Mx';
A = -randn(n);
B = randn(n,1);
%%
AA = [A-B*Mwi*Mxw'      ,-B*Mwi*B';
      -(Mx-Mxw*Mwi*Mxw'),-A'+Mxw*Mwi*B'];

Ei = 0.3*eye(n); Ei = Ei*Ei';

P = eye(n);
N = Ei;

tt = 1.3;
phi = expm(AA*tt);

PN = phi*[P,0*P;N,N];
PNc = mat2cell(PN,[n,n],[n,n]);
P0 = PNc{1,1}; N0 = PNc{2,1};
P1 = PNc{1,2}; N1 = PNc{2,2};

E = N0*P0^-1;
El = (N1 - N0*P0^-1*P1)*P0^-1;

Ell = 2*(N0*P0^-1*P1-N1)*P0^-1*P1*P0^-1;
found = any(eig(Ell)>1e-4) && all(eig(P0)>0) && any(eig(N0)<0);
max(eig(Ell))

figure(1)
clf
hold on

fr = @(r) min(eig([eye(n),zeros(n)]*expm(AA*tt*r)*[P;N]));
xx = linspace(0,1,1000);
rr = arrayfun(fr,xx);
plot(xx,rr)

fr = @(r) min(eig(f_Ell(AA,r,n,P,N)));
xx = linspace(0,3,1000);
rr = arrayfun(fr,xx);
plot(xx,rr)

ylim([-1,0])
%%
ee = 1e-6;
diff_f = @(f,r) (f(r+ee)-f(r))/ee;
fr = @(r) f_Ell(AA,r,n,P,N);
J = diff_f(fr,tt)
eig(J)
%%
dP0 = A*P0-B*Mwi*B'*N0;
dP1 = A*P1-B*Mwi*B'*N1;
dN0 = -Mx*P0-A'*N0;
dN1 = -Mx*P1-A'*N1;

(N0*P0^-1*P1-N1)*P0^-1*B*Mwi*B'*(N0*P0^-1*P1-N1)*P0^-1 +...
(N0*P0^-1*B*Mwi*B'-A') * (N0*P0^-1*P1-N1) * P0^-1*P1*P0^-1 +...
(N0*P0^-1*P1-N1) * P0^-1*P1*P0^-1 * (B*Mwi*B'*N0*P0^-1 - A)

(N0*P0^-1*P1-N1)*P0^-1 * B*Mwi*B' * (N0*P0^-1*P1-N1)*P0^-1
(N0*P0^-1*B*Mwi*B'-A') * (N0*P0^-1*P1-N1) * P0^-1*P1*P0^-1
(N0*P0^-1*P1-N1) * P0^-1*P1*P0^-1 * (B*Mwi*B'*N0*P0^-1 - A)

DD = (N0*P0^-1*P1-N1) * P0^-1*P1*P0^-1
EE = N0*P0^-1

dDD = (N0*P0^-1*P1-N1)*P0^-1 * B*Mwi*B' * (N0*P0^-1*P1-N1)*P0^-1+...
(EE * B*Mwi*B' - A') * DD + DD * (B*Mwi*B' * EE - A)

G = [ P0*pinv(P1)*B*Mwi*B'*pinv(P1)'*P0, (B*Mwi*B' * EE - A);
    (EE * B*Mwi*B' - A'), zeros(2)]

eig([DD,eye(2)]*G*[DD,eye(2)]')

% -> ie: la variation est négative, du coup, on est convex
%%

Ei = randn(2); Ei = Ei*Ei';

ee = 1e-4;

eig(Ei)
comp_E = @(l) comp_E_tmp(l,Ei,phi);


diff_E = @(l) (comp_E(l+ee) -comp_E(l))/ee;
diff2_E = @(l) (diff_E(l+ee) -diff_E(l))/ee;

l = 1;

E2 = diff2_E(l);

[E,El,Ell] = comp_Ell(l,Ei,phi)

min(eig(Ell))
max(eig(Ell))

function E = comp_E_tmp(l,E0,phi)
    P0 = eye(2);
    N0 = E0;

    
    PN = phi*[P0;l*N0];
    PN = mat2cell(PN,[2,2],[2]);
    P0 = PN{1}; N0 = PN{2};

    
    PN = phi*[0*P0;N0];
    PN = mat2cell(PN,[2,2],[2]);
    P1 = PN{1}; N1 = PN{2};

    E = N0*P0^-1;
    El = N1*(P0+l*P1)^-1 - N1*(P0+l*P1)^-1*P1*(P0+l*P1)^-1;
    Ell = N1*P0^-1*(P0+l*P1)*P0^-1 - N1*P0^-1*P1*P0^-1;
end

function [E,El,Ell] = comp_Ell(l,E0,phi)
    P = eye(2);
    N = l*E0;

    
    PN = phi*[P;N];
    PN = mat2cell(PN,[2,2],[2]);
    P0 = PN{1}; N0 = PN{2};

    
    PN = phi*[0*P;N];
    PN = mat2cell(PN,[2,2],[2]);
    P1 = PN{1}; N1 = PN{2};

    E = N0*P0^-1;
    El = N1*P0^-1 - N0*P0^-1*P1*P0^-1;
    Ell = 2*(N0*P0^-1*P1-N1)*P0^-1*P1*P0^-1;
end



function Ell = f_Ell(AA,tt,n,P,N)
    phi = expm(AA*tt);

    PN = phi*[P,0*P;N,N];
    PNc = mat2cell(PN,[n,n],[n,n]);
    P0 = PNc{1,1}; N0 = PNc{2,1};
    P1 = PNc{1,2}; N1 = PNc{2,2};

    E = N0*P0^-1;
    El = (N1 - N0*P0^-1*P1)*P0^-1;

    Ell = 2*(N0*P0^-1*P1-N1)*P0^-1*P1*P0^-1;
end


function El = f_El(AA,tt,n,P,N)
    phi = expm(AA*tt);

    PN = phi*[P,0*P;N,N];
    PNc = mat2cell(PN,[n,n],[n,n]);
    P0 = PNc{1,1}; N0 = PNc{2,1};
    P1 = PNc{1,2}; N1 = PNc{2,2};

    E = N0*P0^-1;
    El = (N1 - N0*P0^-1*P1)*P0^-1;
    Ell = 2*(N0*P0^-1*P1-N1)*P0^-1*P1*P0^-1;
end

