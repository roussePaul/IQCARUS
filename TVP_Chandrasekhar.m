classdef TVP_Chandrasekhar < TVP
    methods
        function obj = TVP_Chandrasekhar(IQC_sys, ti, E0,f0,g0)
            obj = obj@TVP(IQC_sys, ti, E0,f0,g0);
        end
        
        function [E,f,g] = vec_to_Efg(obj,X)
            [E,~,f,g] = obj.vec_to_state(X);
        end
        
        function x = Efg_to_vec(obj,E,f,g)
            [Mx,~,~,Mwi,Mxw,~,~] = obj.IQC_sys.split_M();
            [A,B,~] = obj.IQC_sys.get_LTI();
            dE0 = -E*A - A'*E - Mx + (B'*E+Mxw')'*Mwi*(B'*E+Mxw');
            L = chol(dE0)';
            x = obj.state_to_vec(E,L,f,g);
        end
        
        function dX = ff(obj,t,X,u)
            n = obj.IQC_sys.n;
            [E,L,f,~] = obj.vec_to_state(X);
            
            [Mx,Mu,~,Mwi,Mxw,Mxu,Muw] = obj.IQC_sys.split_M();
            [A,B,Bu] = obj.IQC_sys.get_LTI();
            
            G = [B*Mwi*B',Bu-B*Mwi*Muw';(Bu-B*Mwi*Muw')', Mu - Muw*Mwi*Muw'];
            
            dE = L*L';
            
            aa = A-B*Mwi*Mxw';
            
            dL = (E*B*Mwi*B'-aa')*L;
            df = -A'*f+(Mxu+E*Bu)*u + (E*B+Mxw)*Mwi*(B'*f-Muw'*u);
            dg = [f;u]'*G*[f;u];
            
            dX = obj.state_to_vec(dE,dL,df,dg);
        end
        
        
        function [E,L,f,g] = vec_to_state(obj,X)
            n = obj.IQC_sys.n;
            C = mat2cell(X,[n^2,n^2,n,1],[1]);
            E  = reshape(C{1},n,n);
            L = reshape(C{2},n,n);
            f = C{3}; g = C{4};
        end
        
        function x = state_to_vec(obj,E,L,f,g)
            x = [E(:);L(:);f;g];
        end
    end
end