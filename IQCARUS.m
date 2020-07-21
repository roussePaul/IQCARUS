% implement the simulation of the paraboloids. Each parabolloid is
% propagated using the hamiltonian form of the differential ricatti
% equation

classdef IQCARUS < handle
    properties
        % system definition
        A,B,Bu,M;
        
        % initial paraboloid set
        P0;
        
        Pi = {}; % set of paraboloids
        t % integration time
        
        AA,BB,RR; % hamiltonian dynamic
        
        % parameters
        Ts = 0.01; % sampling of the paraboloids
        N_parab; % maximum number of parabolic sets
        
        Ts_cons = 0.1; % sampling rate of constraint update
        dlambda = 0.1; % steps in the scaled paraboloids
        N_max_add_scaled = 10; % maximum new considered constraints
        
        save_sim = true; % maximum new considered constraints
        
        n_dir; % direction of the contact points to study
        N_dir = 10; % direction of the contact points to study
        
        update_M = []; % update M at T_cons
        
        % system constants
        n; % state space dimension
        m; % disturbance dimension
        p; % input space dimension
    end
    methods
        function obj = IQCARUS(A,B,Bu,M, t0, E0,f0,g0)
            TVP.reset_id();
            obj.A = A;
            obj.B = B;
            obj.Bu = Bu;
            obj.M = M;
            
            obj.n = size(A,1);
            obj.m = size(B,2);
            obj.p = size(Bu,2);
            
            obj.t = t0;
            obj.P0 = {E0,f0,g0};
            
            obj.init_hamil();
            obj.init_Pi();
            obj.init_n_dir();

        end
        
        function rewind(obj, t0)
            obj.t = t0;
            obj.Pi = {};
            obj.init_Pi();
        end
        
        function tvp = new_paraboloid(obj, ti, E0,f0,g0)
            tvp = TVP_Chandrasekhar_LK(obj, ti, E0,f0,g0);
        end
        
        function tvp = add_new_paraboloid(obj, ti, E0,f0,g0)
            tvp = obj.new_paraboloid(ti, E0,f0,g0);
            obj.Pi{end+1} = tvp;
        end
        
        function tvp = add_scaled_paraboloids(obj, ti, P, lambda)
            [E,f,g] = P.Efg_at(ti);
            for l = lambda
                fprintf("[%2.2e] add paraboloid, scaling factor: 1 + %e\n",ti,l-1);
                obj.add_new_paraboloid(ti, l*E,l*f,l*g);
            end
        end
        
        function init_Pi(obj)
            obj.Pi{end+1} = obj.new_paraboloid(0, obj.P0{:});
        end
        
        function init_n_dir(obj)
            n = obj.n;
            if n==1
                obj.n_dir=1;
            else
                n_dir = randn(obj.N_dir,n);
                r = sqrt(sum(n_dir.*n_dir,2));
                obj.n_dir = n_dir./repmat(r,1,n);
            end
        end
        
        % propagate constraintes for dt
        function sim_all_P(obj, tf, u)
            for idP=1:length(obj.Pi)
                P = obj.Pi{idP};
                if  (P.tf>=tf)
                    P.sim([obj.t,tf],u);
                end
            end
        end

        function advance_sim_and_update(obj, tf, u)
            obj.update_Pi(u);
            obj.sim_all_P(tf,u);
            obj.t = tf;
            
            if isempty(obj.update_M)==false
                obj.M = obj.update_M(obj,tf);
            end
            
            if obj.save_sim==false
                obj.reset_trace();
            end
        end
        
        function sim_until(obj, tf, u)
            while obj.t<tf
                obj.advance_sim_and_update(obj.t+obj.Ts_cons,u)
            end
        end
        
        
        function reset_trace(obj)
            for idP=1:length(obj.Pi)
                P = obj.Pi{idP};
                P.reset_trace();
            end
        end
                
                
% parameter definition 
        function [Mx,Mu,Mw,Mwi,Mxw,Mxu,Muw] = split_M(obj)
            M = mat2cell(obj.M,[obj.n,obj.p,obj.m],[obj.n,obj.p,obj.m]);
            Mx = M{1,1};
            Mu = M{2,2};
            Mw = M{3,3}; Mwi = Mw^-1;
            Mxw = M{1,3};
            Muw = M{2,3};
            Mxu = M{1,2};
        end
        
        function [A,B,Bu] = get_LTI(obj)
            A = obj.A; B = obj.B; Bu = obj.Bu;
        end
        
        function [AA,BB,RR] = init_hamil(obj)
            [Mx,Mu,~,Mwi,Mxw,Mxu,Muw] = obj.split_M();
            [A,B,Bu] = obj.get_LTI();
            
            obj.AA = [A-B*Mwi*Mxw'      ,-B*Mwi*B';
                      -(Mx-Mxw*Mwi*Mxw'),-A'+Mxw*Mwi*B'];
            obj.BB = [Bu-B*Mwi*Muw';-(Mxu-Mxw*Mwi*Muw')];
            obj.RR = [Mu-Muw*Mwi*Muw',Muw*Mwi*Mxw'-Mxu', -(Bu'-Muw*Mwi*B')];
        end

        
% return 1 if X belongs to the rechable set
        function is_in_RS(obj,X)
            x = X(:,1:obj.n);
            xq = X(:,end);
        end
        
% return an overapproximation of the maximal energetic level
        function xq = max_xq(obj,t,x)
            xq = Inf*ones(size(x,1),1);
            for idP=1:length(obj.Pi)
                P = obj.Pi{idP};
                if  (P.ti<=t) & (P.tf>=t)
                    [E,f,g] = P.Efg_at(t);
                    if isreal([E,-f;-f',g])==0
                        error('non real values')
                    end
                    xq = min([xq,P.max_xq(t,x)],[],2);
                end
            end
        end
        
        
        function [tmin,tmax,Pmin,Pmax] = get_bounds_section_plane(obj,t,np,xp)
            tmin = [];
            tmax = [];
            Pid = [];
            for idP=1:length(obj.Pi)
                P = obj.Pi{idP};
                if  (P.ti<=t) & (P.tf>=t)
                    p = P.poly_coef_on_plane(t,np,xp);
                    tt = roots(p);
                    if (p(1)>0) & (isreal(tt))
                        tt = sort(tt);
                        tmin = [tmin,tt(1)];
                        Pid  = [Pid,P];
                        tmax = [tmax,tt(2)];
                    end
                end
            end
            [tmin,imin] = max(tmin); Pmin = Pid(imin);
            [tmax,imax] = min(tmax); Pmax = Pid(imax);
        end
        
        function xd = get_center(obj)
            P = obj.Pi{end};
            [E,f,~] = P.Efg_at(Inf);
            xd = E^-1*f;
        end
        
        function [P,L] = unique_P(obj,P,ll)
            [P,~,ic] = unique(P);
            L = [];
            for i=1:size(P)
                L(i) = max(ll(ic==i));
            end
        end
        
        function [X,P,ll] = get_rising_tt_on_boundary(obj,t,u)
            X = [];
            P = [];
            ll = [];
            xd = obj.get_center();
            for nd = obj.n_dir'
                [x1,P1,x2,P2] = get_proj_P(obj,t,nd,xd);
                if isempty(x1)==0
                    dxq1 = P1.eval_dxq(obj.t,x1,u);
                    [E,~,~] = P1.curr_Efg();
                    if (dxq1>0) && (all(eig(E))>0) 
                       X = [X;x1'];
                       P = [P;P1];
                       l1 = P1.eval_lambda(obj.t,x1,u);
                       ll = [ll;l1];
                    end
                end
                if isempty(x2)==0
                    dxq2 = P2.eval_dxq(obj.t,x2,u);
                    [E,~,~] = P2.curr_Efg();
                    if (dxq2>0)  && (all(eig(E))>0) 
                       X = [X;x2'];
                       P = [P;P2];
                       l2 = P2.eval_lambda(obj.t,x2,u);
                       ll = [ll;l2];
                    end
                end
            end
        end
        
        function [x1,P1,x2,P2] = get_proj_P(obj,t,np,xp)
            [t1,t2,P1,P2] = obj.get_bounds_section_plane(t,np,xp);
            if isempty(t1)==0
                x1 = xp+t1*np;
            else
                x1 = [];
            end
            if isempty(t2)==0
                x2 = xp+t2*np;
            else
                x2 = [];
            end
        end
        
        function lambda = get_tabl_lambda(obj,l)
            Nl = (l-1-obj.dlambda)/obj.dlambda;
            if Nl<obj.N_max_add_scaled
                lambda = (1+obj.dlambda):obj.dlambda:l;
            else
                lambda = linspace(1,l,obj.N_max_add_scaled+1);
                lambda = lambda(2:end);
            end
        end
        
        function update_Pi(obj,u)
            [~,P,ll] = obj.get_rising_tt_on_boundary(Inf,u);
            [P,ll] = obj.unique_P(P,ll);
            for i=1:size(P,1)
                lambda = obj.get_tabl_lambda(ll(i));
                obj.add_scaled_paraboloids(obj.t, P(i), lambda);
            end
        end

% -- Plotting functions
% plot the reachable set projection over a given plane
        function [tt,xx,xq,err] = pts_section_plane(obj,t,N,varargin)
            if (obj.n>1) & (nargin<=3) % if the problem is of dimension > 1 and nc is not provided, we can't plot
                error('Not enough input arguments.');
            end
            nargin
            if nargin==5
                np = varargin{1};
                xp = varargin{2};
            end
            if nargin==3
                np = [1];
                xp = [0];
            end
            
            err = 0;

            [tmin,tmax,~,~] = obj.get_bounds_section_plane(t,np,xp);
            if isempty(tmin)
                warning('bounds not found')
                tmin = -1;
                tmax = 1;
                err = 1;
            end
            if isreal([tmin,tmax])==0
                error('non real roots');
            end
            tt = linspace(tmin,tmax,N);
            xx = ([np,xp]*[tt;ones(1,N)])';
            xq = obj.max_xq(t,xx);
        end
        
        function [tt,xx,xq] = parab_pts_section_plane(obj,t,N,varargin)
            if (obj.n>1) & (nargin<=3) % if the problem is of dimension > 1 and nc is not provided, we can't plot
                error('Not enough input arguments.');
            end
            if nargin==5
                np = varargin{1};
                xp = varargin{2};
            end
            if nargin==3
                np = [1];
                xp = [0];
            end
            
            err = 0;

            [tmin,tmax,~,~] = obj.get_bounds_section_plane(t,np,xp);
            if length(tmin)==0
                warning('bounds not found')
                tmin = -1;
                tmax = 1;
                err = 1;
            end
            if isreal([tmin,tmax])==0
                error('non real roots');
            end
            tt = linspace(tmin,tmax,N);
            xx = ([np,xp]*[tt;ones(1,N)])';

            N_parab = size(obj.Pi,2);
            xq = {};
            for i=1:N_parab
                P = obj.Pi{i};
                if (P.ti<=t) & (P.tf>=t)
                    xq{end+1} = P.max_xq(t,xx);
                end
            end
            xq = horzcat(xq{:});
        end
        
        function [conic_cons] = gather_conics_eq_at(obj,t)
            conic_cons = {};
            N_parab = size(obj.Pi,2);
            for i=1:size(obj.Pi,2)
                P = obj.Pi{i};
                if (P.ti<=t) & (P.tf>=t)
                    [E,f,g] = P.Efg_at(t);
                    conic_cons{end+1} = [E,-f;-f',g];
                end
            end
        end
    end
end
