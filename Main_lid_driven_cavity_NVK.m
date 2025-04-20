clc
clear
close all

set(groot,'defaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times New Roman');

addpath("Preevious data/")
%% Initialization
% Domain size and basic informations
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
N = 64;
Nx = N;                         % Number of interior point in x dir at the finest level
Ny = N;                         % Number of interior point in x dir at the finest level
M = 2;
dx1 = (xmax - xmin)/Nx;
dy1 = (ymax - ymin)/Nx;
dt = 1e-3;
t_end = 100;
max_iter_burgers = 1000;
max_iter_PPE = 100;
max_level = floor(log2(N))-2;
Iter_cycle = 15;
Iter_last_cycle = 30;
Iter_smooth = Iter_cycle.*ones(max_level,1);
Iter_smooth(end) = Iter_last_cycle;
tol_PPE = 1e-12;
tol_burgers = 1e-12;
tol_steady_state = 1e-3;
Solver_burgers = "GS";
cmap = [0 0 1;  % Blue (Min)
        1 1 1;  % White (Middle)
        1 0 0]; % Red (Max)
num_color = 256; % Number of colors
customColormap = interp1([1, num_color/2, num_color], cmap, linspace(1, num_color, num_color));

% Viscosity and Reynolds number
u_top = 1;
L = xmax - xmin;
nu = 1e1;
Re = u_top*L/nu;           % Reynolds number

% Mesh grid
x = linspace(xmin-dx1,xmax+dx1,Nx + 3);
y = linspace(ymin-dy1,ymax+dy1,Ny + 3);
[X_edge,Y_edge] = meshgrid(x,y);
[X_cc,Y_cc,X_u_fc,Y_u_fc,X_v_fc,Y_v_fc,Delta_x,Delta_y,dx,dy] = Calculate_coordinate(x,y,Nx,Ny);

% Delta_x: Distance between edge point in x dir
% Delta_y: Distance between edge point in y dir
% dx: Distance between face center point (v_face) in x dir
% dy: Distance between face center point (v_face) in y dir

% Mesh grid for multigrid solver
Nx_level = zeros(max_level,1);
Ny_level = zeros(max_level,1);

for i = 1:max_level
    Nx_level(i) = (Nx)/M^(i-1);
    Ny_level(i) = (Ny)/M^(i-1);
end

dx_level = (xmax - xmin)./(Nx_level - 1);
dy_level = (ymax - ymin)./(Ny_level - 1);


% Initial conditions: cell center
u_cc = zeros(size(X_cc));
v_cc = zeros(size(Y_cc));
p_cc = zeros(size(X_cc));

[u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_top);      % Add BCs

% Get cell center velocity
[U_fc,~] = Get_face_variable(u_cc,Nx,Ny);
[~,V_fc] = Get_face_variable(v_cc,Nx,Ny);

% Velocity at time step: n-1 for 1st step
u_cc_old = u_cc;
v_cc_old = v_cc;
U_fc_old = U_fc;
V_fc_old = V_fc;

%% Time marching
t = 0;
Nt = t_end/dt;

for n = 1:Nt
    t = t + dt; % Update time
    fprintf('Time: %3f s\n',t)

    [u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,p_cc] = NS_solver_NVK(u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,p_cc,dt,Delta_x,Delta_y,dx,dy,Nx,Ny,Re,u_top, ...
                                                                    Solver_burgers,max_iter_burgers,max_iter_PPE,max_level,Nx_level,Ny_level,dx_level,dy_level,Iter_smooth,tol_PPE,tol_burgers);
    res_u = rms(u_cc - u_cc_old,"all");
    res_v = rms(v_cc - v_cc_old,"all");
    fprintf('Residual (u,v): (%e, %e) \n',res_u,res_v)

    if(res_u < tol_steady_state && res_v < tol_steady_state)
        break
    end
    
end

%% Plots
w = calculate_vorticity(u_cc,v_cc,Nx,Ny,dx,dy);
w = -w;

figure(1)
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc(2:end-1,2:end-1),50,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
title(['u: Re = ',num2str(Re),', N = ',num2str(N)],'FontWeight','bold')
xlim([xmin,xmax])
ylim([xmin,xmax])

figure(2)
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),v_cc(2:end-1,2:end-1),50,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
title(['v: Re = ',num2str(Re),', N = ',num2str(N)],'FontWeight','bold')
xlim([xmin,xmax])
ylim([xmin,xmax])

figure(3)
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),w(2:end-1,2:end-1),1000,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
title(['Vorticity: Re = ',num2str(Re),', N = ',num2str(N)],'FontWeight','bold')
xlim([xmin,xmax])
ylim([xmin,xmax])
clim([-6,6])

figure(4);
% set(gcf,'Position',[1200 300 500 400])
% contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),w(2:end-1,2:end-1),100,'LineStyle','none'); hold on
% pcolor(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),w(2:end-1,2:end-1),'FaceColor','interp','EdgeColor','interp'); hold on
% imagesc(X_cc(2,2:end-1),Y_cc(2:end-1,2),w(2:end-1,2:end-1),50); hold on
l = streamslice(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc(2:end-1,2:end-1),v_cc(2:end-1,2:end-1),10); hold on
hold off
% colorbar
% colormap(customColormap)
axis square
set(l,'LineWidth',1)
set(l, 'Color', "k")
% set(l, 'Color', "#0072BD")
title(['Streamline: Re = ',num2str(Re),', N = ',num2str(N)],'FontWeight','bold')
xlim([xmin,xmax])
ylim([xmin,xmax])
% clim([-6,6])

%% Compare velocity

Name_u = "u-Re-" + num2str(Re);
Name_v = "v-Re-" + num2str(Re);
data_u_lit = xlsread("Previous data.xlsx",Name_u);
data_v_lit = xlsread("Previous data.xlsx",Name_v);

Nx_half = floor(Nx/2) + 1;
Ny_half = floor(Ny/2) + 1;

u_half = u_cc(:,Nx_half);
v_half = v_cc(Ny_half,:);
x_half = X_cc(Ny_half,:);
y_half = Y_cc(:,Nx_half);

 % u
 figure(1)
 plot(data_u_lit(:,1),data_u_lit(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Ghia et al. (1982)'); hold on
 plot(u_half,y_half,'color','k','DisplayName','This study'); hold on
 hold off
 xlabel('$ u \ (m/s) $','Interpreter','latex')
 ylabel('$ y \ (m) $','Interpreter','latex')
 grid on
 legend('Location','southeast')
 ylim([0,1])
 title(['Re = ',num2str(Re),': \it{u} \rm{profile}'])
 drawnow

 % v
 figure(2)
 plot(data_v_lit(:,1),data_v_lit(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Ghia et al. (1982)'); hold on
 plot(x_half,v_half,'color','k','DisplayName','This study'); hold on
 hold off
 xlabel('$ x \ (m) $','Interpreter','latex')
 ylabel('$ v \ (m/s) $','Interpreter','latex')
 grid on
 legend('Location','southwest')
 xlim([0,1])
 title(['Re = ',num2str(Re),': \it{v} \rm{profile}'])
 drawnow