clear all; clc; clf;

%------[Configurations]------%


% Spacetime and particle type %
sptm        = sptm_Ghosh; % Spacetime
ep          = 1; %Particle type; 0 - null, 1 - timelike


% Initial conditions %
r0          = 8;
th0         = pi/2;
phi0        = 0.7484;

ur0         = 0;
uth0        = 0.2/10;
uphi0       = -0.6/10;


% Integration settings %
int_method  = integ_heun; % Integration method
dt          = 0.5; % Time step
dxdu        = 1e-5; % Differentiation interval
t_end       = 5000; % End time


%------[Configurations]------%


e = [ep, dt, dxdu];

x = [0, r0, th0, phi0];
u_con = [dt/abs(dt), ur0, uth0, uphi0];

x_pl = x;
u = u_cov(x,u_con,sptm,eps);

nt = abs(t_end/dt);
H_init = H_hamiltonian(x,u,e,sptm);
dev_H = zeros(1,nt);

tic
for i=1:nt
    [x,u] = int_method.integrate(x,u,e,sptm);
    x_pl(i+1,:) = x;
    H_iter = H_hamiltonian(x,u,e,sptm);
    dev_H(i) = (H_init - H_iter)/H_init;
end
toc

x_plot = x_pl(:,2).*cos(x_pl(:,4)).*sin(x_pl(:,3));
y_plot = x_pl(:,2).*sin(x_pl(:,4)).*sin(x_pl(:,3));
z_plot = x_pl(:,2).*cos(x_pl(:,3));

plim = 15;

subplot(1,2,1)
plot3(x_plot,y_plot,z_plot)
xlim([-plim plim]); ylim([-plim plim]), zlim([-plim plim])
xlabel('x'); ylabel('y'); zlabel('z');
pbaspect([1,1,1])
grid on

subplot(1,2,2)
plot(linspace(0,nt,nt),dev_H)
xlabel('iter'); ylabel('$\Delta H / H_0$','Interpreter','Latex')

