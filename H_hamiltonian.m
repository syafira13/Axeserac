function Hamilton = H_hamiltonian(x,u,e,sptm)

ut = u(1); ur = u(2); uth = u(3); uph = u(4);

eps = e(1);
g = sptm.g_uv(x);

alph = sqrt(...
    -g(1,1) + (g(1,4)^2) / (g(4,4) ) ...
    );
yuu = ((ur^2)/g(2,2)) + ((uth^2)/g(3,3)) + ((uph^2)/g(4,4));
betu = g(1,4)*uph/g(4,4);

Hamilton = ( alph*sqrt(yuu + eps) ) - betu;

