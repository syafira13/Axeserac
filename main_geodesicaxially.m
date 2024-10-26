clear all; clc; clf;

sptm = Ghosh;

eps = 1; 

dt = 1; dxdu = 1e-5;
e = [eps, dt, dxdu];

% x = [0, 8, pi/2, 0]; 
x = [0, 8.4209, pi/2, 0.7484]
x_pl = x;
u_con = [dt/abs(dt), 0, 0, -0.6/10];
% u = u_cov(x,u_con,sptm);
u = [1, 0.0155, 0 , -4.2118]

t_end = 500; nt = abs(t_end/dt);
H_init = H(x,u,e,sptm);
dev_H = zeros(1,nt);


tic
for i=1:nt
    [x,u] = heunm(x,u,e,sptm);
    x_pl(i+1,:) = x;
    H_iter = H(x,u,e,sptm);

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

