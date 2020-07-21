if ~exist('redef_params')
    order_red = 19;
    reduce_sys = true;
    Tmax = 1.8;
    show_plots = true;
    N_max_add_parab = 0;
    save_sim = true;
end
%%
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('HE7');

D22 = zeros(ny,nu);
P = {A,B1,B;C1,D11,D12;C,D21,D22}

P = ss(A,[B1,B],[C1;C],[D11,D12;D21,D22]);
[K,CL,GAM,INFO] = h2syn(P,ny,nu);

%%
sys = CL;

if reduce_sys
    sys = balred(CL,order_red);
end
sys.A

ee = abs(eig(sys.A));
min(ee),max(ee)
%%

A = sys.A;
nx = size(A,1);
Ix = eye(nx);


Bw = Ix;
B = sys.B;

nu = size(B,2);
nw = size(Bw,2);
zu = zeros(nu);
zw = zeros(nw);
Iw = eye(nw);

mw = -5e3;
M = blkdiag(Ix,zu,mw*Iw);

id_xw = [ones(1,nx),zeros(1,nu),ones(1,nw)]==1;
M_xw = M(id_xw,id_xw);
mm = mat2cell(M_xw,[nx,nw],[nx,nw]);
Mx = mm{1,1}; Mxw = mm{1,2}; Mw = mm{2,2};

[E,L,G] = care(A,Bw,Mx,Mw,Mxw,Ix);

%%
zw = zeros(nw);

dE = [-E*A-A'*E, (Bw'*E)'; (Bw'*E),  zw] - M_xw;

min(eig(dE))
max(eig(dE))

dE = -E*A-A'*E - Mx+ (Bw'*E)'*Mw^-1*(Bw'*E);
norm(dE)
norm(dE)/norm(E)

%%

gg = 10;

x0 = zeros(nx,1);
E0 = E;
f0 = E0*x0;
g0 = (x0'*E0*x0-0.001);

TVP.reset_id();
iqc = IQCARUS(A,Bw,B,M, 0, gg*E0,gg*f0,gg*g0);

iqc.N_dir = 100;
iqc.init_n_dir();

iqc.Ts = 1e-4;
iqc.Ts_cons = 0.05;
iqc.N_max_add_scaled = N_max_add_parab;

iqc.save_sim = save_sim;

uu = @(t) exp(-t)*ones(nu,1);
tic
iqc.sim_until(Tmax,uu);
t_compute_sim = toc

%%
if show_plots
    P = iqc.Pi{1};
    [E,f,g] = P.curr_Efg();
    eig(E)

    [E,L,K,f,g] = P.vec_to_state(P.X(end,:)');
    Ln = norm(L)
    Kn = norm(K)

    %%
    for t=linspace(0,0.99*iqc.t,100)
        figure(1)
        clf
        hold on
        for i=1:length(iqc.Pi)
            P = iqc.Pi{i};
    %       if P.ti<=t
    %           break
    %       end
    %       P = iqc.Pi{2};
            [E,f,g] = P.Efg_at(t);

            if any(eig(E)<0)
                "break"
                break
            end

            np = [1;zeros(nx-1,1)];

            [tt,xx,xq] = iqc.parab_pts_section_plane(t,300,np,E^-1*f);

            plot(xx*np,xq(:,1),'b')

            [tt,xx,xq] = iqc.pts_section_plane(t,300,np,E^-1*f);
            plot(xx*np,xq,'--k','linewidth',2)

            xc = (E^-1*f)'*np;
            axis([xc-1.5*(xc-min(xx*np)),xc-1.5*(xc-max(xx*np)),0,f'*E^-1*f-g])
        end

        fprintf('%.6f\n',f'*E^-1*f-g);
        pause(0.01)
    end

    %%
    Xpos_sup = [];
    Xpos_inf = [];

    T = linspace(0,0.99*iqc.t,1000);


    ws = zeros(nw,1);
    opt = odeset('RelTol',1e-14,'AbsTol',1e-14);
    [t,y] = ode113(@(t,x) [A*x(1:nx)+B*uu(t)+Bw*ws;[x(1:nx);uu(t);ws]'*M*[x(1:nx);uu(t);ws]],T, [x0;0],opt);

    hh = [];

    for idt=1:length(T)

        t = T(idt);
        z2N = zeros(nx,1);
        xmax = z2N+Inf; xmin = z2N-Inf;
        for i=1:length(iqc.Pi)
            P = iqc.Pi{i};
            if any([t<P.ti,t>P.tf])
                break;
            end
            [E,f,g] = P.Efg_at(t);

            xt = y(idt,1:nx)';
            xqt = y(idt,end)';
            hh(end+1) = xt'*E*xt-2*f'*xt+g+xqt;

            if any(eig(E)<0)
                "break"
                break
            end

            xmax_t = z2N;
            xmin_t = z2N;
            xp = E^-1*f;
            for id_x=1:nx
                np = z2N; np(id_x) = 1;
                xmax_t(id_x) =  sqrt(np'*E^-1*np)+np'*xp;
                xmin_t(id_x) = -sqrt(np'*E^-1*np)+np'*xp;
            end
            xmax = min([xmax,xmax_t],[],2);
            xmin = max([xmin,xmin_t],[],2);
        end
        Xpos_sup(end+1,:) = xmax';
        Xpos_inf(end+1,:) = xmin';
    end

    %%
    for i=1:nx
        figure(i)
        clf
        hold on
        xsup = Xpos_sup(:,i)';
        xinf = Xpos_inf(:,i)';
        xf = [xsup,fliplr(xinf)];
        tf = [T,fliplr(T)];
        fill(tf,xf,'r');
        plot(T,y(:,i))
    end

    xq = y(:,end);
    figure(nx+1)
    plot(t,xq)

    figure(nx+2)
    plot(t,hh)

    %%
    min(xq)
end
%%
t_compute_sim