classdef integ_leapfrog
    methods(Static)
        function [xpone,upone] = integrate(x,u,e)
        
            dt = e(2);
            
            uphalf = u - dt*H_dHdx(x,u,e)/2;
            xpone = x + dt*H_dHdu(x,uphalf,e);
            upone = uphalf - dt*H_dHdx(xpone,uphalf,e)/2;
        end
    end
end