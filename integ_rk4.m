classdef integ_rk4
    methods(Static)
        function [xpone,upone] = integrate(x,u,e,sptm)
        
            dt = e(2);
            
            k1 = dt*H_dHdu(x,u,e,sptm);
            l1 = -dt*H_dHdx(x,u,e,sptm);
            
            k2 = dt*H_dHdu(x+k1/2,u+l1/2,e,sptm);
            l2 = -dt*H_dHdx(x+k1/2,u+l1/2,e,sptm);
            
            k3 = dt*H_dHdu(x+k2/2,u+l2/2,e,sptm);
            l3 = -dt*H_dHdx(x+k2/2,u+l2/2,e,sptm);
            
            k4 = dt*H_dHdu(x+k3,u+l3,e,sptm);
            l4 = -dt*H_dHdx(x+k3,u+l3,e,sptm);
            
            xpone = x + (k1 + 2*k2 + 2*k3 + k4)/6;
            upone = u + (l1 + 2*l2 + 2*l3 + l4)/6;

        end
    end
end