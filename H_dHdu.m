function d = H_dHdu(x,u,e,sptm)

du = e(3);

dur = [u(1), u(2) + du, u(3), u(4)];
duth = [u(1), u(2), u(3) + du, u(4)];
duph = [u(1), u(2), u(3), u(4) + du];

d = [1,...
    (H_hamiltonian(x,dur,e,sptm) - H_hamiltonian(x,u,e,sptm))/du,...
    (H_hamiltonian(x,duth,e,sptm) - H_hamiltonian(x,u,e,sptm))/du,...
    (H_hamiltonian(x,duph,e,sptm) - H_hamiltonian(x,u,e,sptm))/du];