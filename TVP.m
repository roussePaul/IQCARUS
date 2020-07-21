classdef TVP < handle
    properties
        ti,tf,t,tt; % initial and final time
        
        X; % array of the states 
        
        IQC_sys; % associated IQC system
        stop_sim;
    end
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        id;
    end

    methods ( Static, Access = 'public' )
        function reset_id()
            persistent stamp;
            stamp = 0;
        end
            
        function result = increment()
            persistent stamp;
            if isempty( stamp )
                stamp = 0;
            end
            stamp = stamp + uint32(1);
            result = stamp;
        end
    end  

    methods
        function obj = TVP(IQC_sys, ti, E0,f0,g0)
            obj.IQC_sys = IQC_sys;
            
            obj.id = TVP.increment();
            
            obj.ti = ti;
            obj.tf = Inf;
            obj.t = ti;
            obj.tt = [ti];
            
            obj.X = [obj.Efg_to_vec(E0,f0,g0)'];
            
            obj.stop_sim = 0;
        end
        
        function advance_sim(obj,u)
            
        end
        
        function idt = get_idt(obj,t)
            if isinf(t)
                idt = -1;
            elseif t>=obj.t
                idt = -1;
            else
                idt = 1+round((t-obj.ti)/obj.IQC_sys.Ts);
            end
        end
        
        function idt = last_idt(obj)
            idt = size(obj.T,1);
        end

        function [E,f,g] = Efg_at(obj,t)
            idt = obj.get_idt(t);
            if idt>0
                [E,f,g] = obj.vec_to_Efg(obj.X(idt,:)');
            elseif idt==-1
                [E,f,g] = obj.vec_to_Efg(obj.X(end,:)');
            end
        end

        function [E,f,g] = curr_Efg(obj)
            [E,f,g] = obj.vec_to_Efg(obj.X(end,:)');
        end

% return 1 if X belongs to the rechable set
        function is_in_RS(obj,X)
            x = X(:,1:obj.n);
            xq = X(:,end);
        end

% return an overapproximation of the maximal energetic level
        function xq = max_xq(obj,t,x)
            [E,f,g] = obj.Efg_at(t);
            xq = -( sum((x*E).*(x),2) - 2*x*f + g);
        end
        
        function ws = eval_ws(obj,t,x,u)
            [E,f,~] = obj.Efg_at(t);
            [~,~,~,Mwi,Mxw,~,Muw] = obj.IQC_sys.split_M();
            [~,B,~] = obj.IQC_sys.get_LTI();
            ws = -Mwi*(B'*(E*x-f) + Mxw'*x + Muw'*u(t));
        end
        
        function dxq = eval_dxq(obj,t,x,u)
            ws = obj.eval_ws(t,x,u);
            dxq = [x;u(t);ws]'*obj.IQC_sys.M*[x;u(t);ws];
        end
        
        
        function l = eval_lambda(obj,t,x,u)
            [E,f,~] = obj.Efg_at(t);
            [~,~,~,Mwi,Mxw,~,Muw] = obj.IQC_sys.split_M();
            [~,B,~] = obj.IQC_sys.get_LTI();
            ut = u(t);
            Mdxq = [0*x,x;0*ut,ut;-Mwi*B'*(E*x-f), -Mwi*(Mxw'*x + Muw'*ut)]'*obj.IQC_sys.M*[0*x,x;0*ut,ut;-Mwi*B'*(E*x-f), -Mwi*(Mxw'*x + Muw'*ut)];
            p = [Mdxq(1,1),Mdxq(1,2)*2,Mdxq(2,2)];
            l = max(roots(p));
        end
        
        function sim(obj,T,u)
            x0 = obj.X(end,:)';
            Ts = obj.IQC_sys.Ts;
            tt = T(1):Ts:T(2);
            opt = odeset('RelTol',5e-14,'AbsTol',5e-14);
            [tt,xx] = ode113(@(t,x) obj.ff(t,x,u(t)), tt, x0,opt);
            xx = xx(2:end,:);
            obj.X = [obj.X;xx];
            obj.t = tt(end);
            obj.tt = [obj.tt,(tt(2:end))'];
            obj.rescale_last_state();
            if tt(end)~=T(2)
                tt(end)~=T(2)
                tt(end)-T(2)
                xx(end,:)
%                 sdfsgf
                obj.tf = tt(end);
            end
            if obj.tt(end) ~= obj.t
                fprintf("***\n");
                obj.t
                tt
                fprintf("***\n");
            end
        end

% plot utility functions
% function to work on projection of the reachable set
% x = xp + np*t
% [xp + np*t;1] = [np,xp;0,1]*[t;1]
        function p = poly_coef_on_plane(obj,t,np,xp)
            [E,f,g] = obj.Efg_at(t);
            Pp = [np,xp;0,1]'*[E,-f;-f',g]*[np,xp;0,1];
            p = [Pp(1,1),2*Pp(1,2),Pp(2,2)];
        end
        
        function rescale_last_state(obj)
        end
        
        function reset_trace(obj)
            obj.X = [obj.X(end,:)];
            
            obj.ti = obj.tt(end);
            obj.tt = [obj.ti];
        end
    
    end

    methods (Abstract=true)
        [E,f,g] = vec_to_Efg(obj,x);
        x = Efg_to_vec(obj,E,f,g);
        dX = ff(obj,x,u);
    end
end