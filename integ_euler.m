classdef integ_euler
    methods(Static)
        function [xpone,upone] = integrate(x,u,e)
        
            dt = e(2);
            
            xpone = x + dt*H_dHdu(x,u,e);
            upone = u - dt*H_dHdx(x,u,e);

        end
    end
end