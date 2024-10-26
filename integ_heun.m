function [xpone,upone] = heunm(x,u,e,sptm)

dt = e(2);

xpone = x + dt*dHdu(x,u,e,sptm);
upone = u - dt*dHdx(x,u,e,sptm);

for i = 1:5
    xpone = x + dt*(dHdu(x,u,e,sptm) + dHdu(xpone,upone,e,sptm))/2;
    upone = u - dt*(dHdx(x,u,e,sptm) + dHdx(xpone,upone,e,sptm))/2;
end