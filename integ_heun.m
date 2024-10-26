classdef integ_heun
    methods(Static)
        function [xpone,upone] =integrate(x,u,e,sptm)
    
            dt = e(2);
            
            xpone = x + dt*H_dHdu(x,u,e,sptm);
            upone = u - dt*H_dHdx(x,u,e,sptm);
            
            for i = 1:5
                xpone = x + dt*(H_dHdu(x,u,e,sptm) + H_dHdu(xpone,upone,e,sptm))/2;
                upone = u - dt*(H_dHdx(x,u,e,sptm) + H_dHdx(xpone,upone,e,sptm))/2;
            end

        end
    end
end