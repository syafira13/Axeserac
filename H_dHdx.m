function d = H_dHdx(x,u,e,sptm)

dx = e(3);

dxr = [x(1), x(2) + dx, x(3), x(4)];
dxth = [x(1), x(2), x(3) + dx, x(4)];
dxph = [x(1), x(2), x(3), x(4) + dx];

d = [0,...
    (H_hamiltonian(dxr,u,e,sptm) - H_hamiltonian(x,u,e,sptm))/dx,...
    (H_hamiltonian(dxth,u,e,sptm) - H_hamiltonian(x,u,e,sptm))/dx,...
    (H_hamiltonian(dxph,u,e,sptm) - H_hamiltonian(x,u,e,sptm))/dx];