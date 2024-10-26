clear all; clc; clf;

% rgbs = linspace(0,1,255);
% rgbs = rgbs';
% maps = [rgbs rgbs.^(1.7) rgbs.^(3)];

%------[Configurations]------%

% Screen and spacetime %
sptm                = sptm_KerrNewman; % Spacetime
Resolution          = [100 100]; % Screen resolution (pixel)
FoV                 = 12*sptm.M; % Field of view
ScreenDistance      = 0.001;

% Initial conditions %
th0                 = pi/2-0.3; % Initial inclination
r0                  = 500; % Observer distance


% Accretion disk properties %
Proret              = 1; % Retrograde (-1) or prograde (1) motion of the accretion disk
doppler_redshift    = true; % Calculate doppler redshift? if false, then only calculate gravitational redshift
all_redshift        = false; % Calculate redshift effects?
thick_accretion     = true; % Optically thick accretion disk?

yg                  = -2; % Intensity \gamma
mu                  = disk_r_isco(3*sptm.M,Proret,sptm); % Intensity \mu
sig                 = sptm.M/4; % Intensity \sigma

lower_R             = mu-0.1; % Accretion disk inner radius
upper_R             = 20*sptm.M; % Accretion disk outer radius


% Integration settings %
err_low             = 1e-10; % Lowest error
err_high            = 1e-7; % Highest error
dxdu                = 1e-6; % Spatial differentiation interval
int_method          = integ_rk4;



%------[Configurations]------%


dt = -1; 
e = [eps, dt, dxdu];
eh = [eps, dt/2, dxdu];

dr0 = 1;
x0 = [0, r0, th0, 0];

param = [yg, mu, sig];

ScreenHeight_y  = ScreenDistance * FoV / r0;
ScreenHeight_z  = ScreenDistance * FoV / r0;
ymage           = linspace(-ScreenHeight_y, ScreenHeight_y, Resolution(1));
zmage           = linspace(-ScreenHeight_z, ScreenHeight_z, Resolution(2));

y_ob = x0(2)*sin(x0(4))*sin(x0(3));
x_ob = x0(2)*cos(x0(4))*sin(x0(3));
z_ob = x0(2)*cos(x0(3));

dz_init         = dr0 * zmage / ScreenDistance;
dy_init         = dr0 * ymage / ScreenDistance;
dph0            = dy_init / (r0*sin(x0(3)));
dth0            = dz_init / r0;

Image           = zeros(Resolution(2),Resolution(1));
Phirec          = zeros(Resolution(2),Resolution(1));
t_end = 10000000; nt = abs(t_end/dt);
tic
for i = 1:Resolution(1)
    parfor j = 1:Resolution(2)
        dt = -1;
        x = x0;

        xcor = x(2)*cos(x(4))*sin(x(3));
        ycor = x(2)*sin(x(4))*sin(x(3));
        zcor = x(2)*cos(x(3));

        u_con = [dt, dr0, -dth0(j), -dph0(i)];
        u = u_cov(x,u_con,sptm,eps);
        H_init = H_hamiltonian(x,u,e,sptm);
        H_e = H_init;
        Inten = 0;

        t_iter = 0;
        iter = 1;

        breaking = false;

        phi_rec = 0;
        t_rec = 0;
        while t_iter<t_end
            
            [x1,u1] = int_method.integrate(x,u,[eps, dt, dxdu],sptm);
            [xh1,uh1] = int_method.integrate(x,u,[eps, dt/2, dxdu],sptm);
            [xh2,uh2] = int_method.integrate(xh1,uh1,[eps, dt/2, dxdu],sptm);

            err = abs((xh2 - x1)/xh2);
            err_dom = max(err);

            if err_dom > err_high
                dt = dt/2;
            elseif max(abs(imag(xh2)),[],'all')>0
                break
            else
                xcor_up = xh2(2)*cos(xh2(4))*sin(x(3));
                ycor_up = xh2(2)*sin(xh2(4))*sin(x(3));
                zcor_up = xh2(2)*cos(xh2(3));

                if zcor_up*zcor < 0
                    r_det = x(2) - (zcor * ( (xh2(2) - x(2)) / (zcor_up - zcor) ) );
                    th_det = pi/2;

                    if and(r_det>lower_R, r_det<upper_R)
                        iplus = disk_intensity_profile(r_det,param);
    
                        gtt = sptm.g_uv_comp([0,r_det,th_det,0],1,1);
                        gtp = sptm.g_uv_comp([0,r_det,th_det,0],1,4);
                        gpp = sptm.g_uv_comp([0,r_det,th_det,0],4,4);
                        omeg = disk_om(r_det,Proret,sptm);

                        uphi = doppler_redshift*omeg/(sqrt(-gtt -( (2*gtp + gpp*omeg)*omeg ) ));
                        ut = (gtp*uphi + sqrt(-gtt + (gtp^2 - gtt*gpp)*uphi^2))/(-gtt);
    
                        g = (all_redshift*1/(ut - uh2(4)*uphi/H_init)) + (not(all_redshift)*1);
                        
                        Inten = Inten + iplus*g^4;
                        breaking = thick_accretion;

                        %-------[For further analysis]-------%
                        % t_rec = t_iter;
                        % H_e = H_hamiltonian(xh2,uh2,e,sptm);
                        % phi_rec = sqrt(xh2(4)^2 + xh2(3)^2);
                        %-------[For further analysis]-------%
                    end
                end
                xcor = xcor_up;
                ycor = ycor_up;
                zcor = zcor_up;

                x = xh2; u = uh2;

                iter = iter+1;
                t_iter = t_iter + abs(dt);
                
                if sptm.g_uv_comp(x,2,2)<=0
                    break
                end
                if x(3)>r0
                    break
                end
                if abs(x(4)) > 40*pi
                    break
                end
                if err_dom > err_low
                    dt = dt/2;
                else
                    dt = dt*2;
                end
            end

            %-------[For further analysis]-------%
            % H_err(j,i) = abs((H_e - H_init)/H_init);
            % Phirec(j,i) = phi_rec;
            % trec(j,i) = t_rec;
            %-------[For further analysis]-------%

            Image(j,i) = Inten;
            if breaking
                break
            end
        end

    end
    imagesc(Image)
    pbaspect([1 1 1])
    pause(0.0001)
    fprintf('Calc finish for pixel %i\n',i)
end
toc

% subplot(2,2,1)
imagesc(Image); pbaspect([1 1 1])

%-------[For further analysis]-------%
% subplot(2,2,2)
% imagesc(abs(real(Phirec))); pbaspect([1 1 1])
% subplot(2,2,3)
% imagesc(H_err); pbaspect([1 1 1])
% subplot(2,2,4)
% imagesc(trec); pbaspect([1 1 1])
%-------[For further analysis]-------%


%-------[Write the data (txt)]-------%

% writematrix(Image,'Image_redshifted_k001horizonless')
% writematrix(Phirec,'Phirec_redshifted_k001horizonless')
% writematrix(trec,'trec_redshifted_k001horizonless')
% writematrix(H_err,'H_err_redshifted_k001horizonless')


