function up = u_cov(x,u,sptm,eps)

guv0 = sptm.g_uv(x);

alph = sqrt(...
    -guv0(1,1) + (guv0(1,4)^2) / (guv0(4,4) ) ...
    );

ur = guv0(2,2) * u(2);
uth = guv0(3,3) * u(3);

uph1 = -alph^2*guv0(4,4)^2*u(4);
uph2 = guv0(4,4)*(guv0(1,4)^2 - alph^2*guv0(4,4));
uph3 =  (-guv0(1,4)^2*((ur^2/guv0(2,2)) +...
        (uth^2/guv0(3,3)) + eps))...
        + alph^2*guv0(4,4)^2*u(4)^2;
uph4 = alph^4*guv0(4,4)^4*u(4)^2;
uph5 = guv0(1,4)^2 - alph^2*guv0(4,4);

uph = (uph1 + sqrt(abs(uph2*uph3 + uph4)))/uph5;

up = [0, ur, uth, uph];

