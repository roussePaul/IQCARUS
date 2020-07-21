classdef TVP_Hamil < TVP
    methods
        function obj = TVP_Hamil(IQC_sys, ti, E0,f0,g0)
            obj = obj@TVP(IQC_sys, ti, E0,f0,g0);
        end
        
        function [E,f,g] = vec_to_Efg(obj,X)
            n = obj.IQC_sys.n;
            [Pt,Nt,pnc,uHPN,upnc] = vec_to_state(obj,X);
            
            pnc= mat2cell(pnc,[n,n],[1]); pct = pnc{1}; nct = pnc{2};
            Pti = Pt^-1;
            
            E = Nt*Pti;
            f = -1/2*(Pti'*uHPN  + (- Nt*Pti*pct + nct));
            g = -uHPN'*Pti*pct + upnc;
        end

        function x = Efg_to_vec(obj,E,f,g)
            n = obj.IQC_sys.n;
            zn = zeros(n,1);
            P = eye(n); N = E;
            pnc = -[zn;f]; uHPN = -f;
            upnc = g;

            x = obj.state_to_vec(P,N,pnc,uHPN,upnc);
        end
        
        function dX = ff(obj,t,X,u)
            iqc = obj.IQC_sys;
            n = iqc.n; p = iqc.p;
            zuPN = zeros(p,n);
            [P,N,pnc,~,~] = obj.vec_to_state(X);
            
            dPN = iqc.AA*[P;N];
            dpnc = iqc.AA*pnc + iqc.BB*u;
            duHPN = (u'*iqc.RR*[zuPN;P;N])';
            dupnc = u'*iqc.RR*[u;pnc];
            
            dPN = mat2cell(dPN,[n,n],[n]); dP = dPN{1}; dN = dPN{2};
            dX = obj.state_to_vec(dP,dN,dpnc,duHPN,dupnc);
        end
        
        function X = state_to_vec(obj, P,N,pnc,uHPN,upnc)
            X = [P(:);N(:);pnc;uHPN;upnc];
        end
        
        function [P,N,pnc,uHPN,upnc] = vec_to_state(obj,X)
            n = obj.IQC_sys.n;
            C = mat2cell(X,[n^2,n^2,2*n,n,1],[1]);
            P = reshape(C{1},n,n); N = reshape(C{2},n,n);
            pnc = C{3}; uHPN = C{4}; upnc = C{5};
        end
        
    end
end