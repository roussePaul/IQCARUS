classdef TVP_Ricatti < TVP
    methods
        function obj = TVP_Ricatti(IQC_sys, ti, E0,f0,g0)
            obj = obj@TVP(IQC_sys, ti, E0,f0,g0);
        end
        
        function [E,f,g] = vec_to_Efg(obj,X)
            n = obj.IQC_sys.n;
            C = mat2cell(X,[n^2,n,1],[1]);
            E = reshape(C{1},n,n);
            f = C{2}; g = C{3};
        end
        
        function x = Efg_to_vec(obj,E,f,g)
            x = [E(:);f;g];
        end
        
        function dX = ff(obj,t,X,u)
            n = obj.IQC_sys.n;
            [E,f,~] = obj.vec_to_state(X);

            [Mx,Mu,~,Mwi,Mxw,Mxu,Muw] = obj.IQC_sys.split_M();
            [A,B,Bu] = obj.IQC_sys.get_LTI();

            G = [B*Mwi*B',Bu-B*Mwi*Muw';(Bu-B*Mwi*Muw')', -Mu + Muw*Mwi*Muw'];

            dE = -E*A - A'*E - Mx + (B'*E+Mxw')'*Mwi*(B'*E+Mxw');
            df = -A'*f+(Mxu+E*Bu)*u + (E*B+Mxw)*Mwi*(B'*f-Muw'*u);

            fprintf("min eigs %f\n",min(real(eig(-A'+(E*B+Mxw)*Mwi*(B'*f)))));

            dg = [f;u]'*G*[f;u];
            
            dX = obj.state_to_vec(dE,df,dg);
        end
        
        function [E,f,g] = vec_to_state(obj,X)
            [E,f,g] = obj.vec_to_Efg(X);
        end
        
        function x = state_to_vec(obj,E,f,g)
            x = obj.Efg_to_vec(E,f,g);
        end
    end
end