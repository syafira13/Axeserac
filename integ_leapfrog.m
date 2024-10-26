function [xpone,upone] = leapfrog(x,u,e)

dt = e(2);

uphalf = u - dt*dHdx(x,u,e)/2;
xpone = x + dt*dHdu(x,uphalf,e);
upone = uphalf - dt*dHdx(xpone,uphalf,e)/2;