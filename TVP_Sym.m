classdef TVP_Sym < TVP
    methods
        function obj = TVP_Sym(IQC_sys, ti, E0,f0,g0)
            obj = obj@TVP(IQC_sys, ti, E0,f0,g0);
        end
        
        function [E,f,g] = vec_to_Efg(obj,X)
            [P,~,D,f,g] = obj.vec_to_state(X);
            E = P*D*P';
        end
        
        function rescale_last_state(obj)
            [P,Pi,D,f,g] = obj.vec_to_state(obj.X(end,:)');
            E = P*D*P';
            [Ew,Ev] = eig(E);
            J = (Ev>0) - (Ev<0);
            Ev = sqrt(abs(Ev));
            P = Ew*Ev;
            Pi = P^-1;
            D = J;
            x = obj.state_to_vec(P,Pi,D,f,g);
            obj.X(end,:) = x';
        end
        
        function x = Efg_to_vec(obj,E,f,g)
            [Ew,Ev] = eig(E);
            J = (Ev>0) - (Ev<0);
            Ev = sqrt(abs(Ev));
            P = Ew*Ev;
            Pi = P^-1;
            D = J;
            x = obj.state_to_vec(P,Pi,D,f,g);
        end
        
        function dX = ff(obj,t,X,u)
            n = obj.IQC_sys.n;
            [P,Pi,D,f,~] = obj.vec_to_state(X);
            
            [Mx,Mu,~,Mwi,Mxw,Mxu,Muw] = obj.IQC_sys.split_M();
            [A,B,Bu] = obj.IQC_sys.get_LTI();
            
            G = [B*Mwi*B'       ,Bu-B*Mwi*Muw';
                (Bu-B*Mwi*Muw')', -Mu + Muw*Mwi*Muw'];
            
            
            E = P*D*P';
            
            dP = (-A+B*Mwi*Mxw')'*P;
            dPi = -Pi*(-A+B*Mwi*Mxw')';
            dD = -Pi*(Mx-Mxw*Mwi*Mxw')*Pi' + D*P'*B*Mwi*B'*P*D;
            
            df = -A'*f+(Mxu+E*Bu)*u + (E*B+Mxw)*Mwi*(B'*f-Muw'*u);
            dg = [f;u]'*G*[f;u];
            
            dX = obj.state_to_vec(dP,dPi,dD,df,dg);
        end
        
        
        function [P,Pi,D,f,g] = vec_to_state(obj,X)
            n = obj.IQC_sys.n;
            C = mat2cell(X,[n^2,n^2,n^2,n,1],[1]);
            P  = reshape(C{1},n,n);
            Pi = reshape(C{2},n,n);
            D  = reshape(C{3},n,n);
            f = C{4}; g = C{5};
        end

        function x = state_to_vec(obj,P,Pi,D,f,g)
            x = [P(:);Pi(:);D(:);f;g];
        end
    end
end