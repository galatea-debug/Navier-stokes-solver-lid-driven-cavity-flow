clc
clear
close all

set(groot,'defaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultAxesFontName','Times New Roman');
addpath("Preevious data/")

data_u_lit = xlsread("Previous data.xlsx","u-Re-100");
data_v_lit = xlsread("Previous data.xlsx","v-Re-100");

%% Initialization
% Domain size and basic informations
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;

% Viscosity and Reynolds number
u_top = 1;
L = xmax - xmin;
nu = 1e-2;
Re = u_top*L/nu;           % Reynolds number

% Other parameters
N = [16,32,64,128];
dt = 1e-3;
t_end = 100;
max_iter = 100;
Iter_cycle = 10;
Iter_last_cycle = 25;
tol_PPE = 1e-12;
tol_burgers = 1e-12;
tol_steady_state = 1e-6;
Solver_burgers = "GS";
color_line = ["k","r","b","g"];

for num_N = 1:length(N)  
    % Mesh grid
    Nx = N(num_N);                         % Number of interior point in x dir at the finest level
    Ny = N(num_N);                         % Number of interior point in x dir at the finest level
    M = 2;
    dx1 = (xmax - xmin)/Nx;
    dy1 = (ymax - ymin)/Nx;
    x = linspace(xmin-dx1,xmax+dx1,Nx + 3);
    y = linspace(ymin-dy1,ymax+dy1,Ny + 3);
    [X_edge,Y_edge] = meshgrid(x,y);
    [X_cc,Y_cc,X_u_fc,Y_u_fc,X_v_fc,Y_v_fc,Delta_x,Delta_y,dx,dy] = Calculate_coordinate(x,y,Nx,Ny);
 
    % Delta_x: Distance between edge point in x dir
    % Delta_y: Distance between edge point in y dir
    % dx: Distance between face center point (v_face) in x dir
    % dy: Distance between face center point (v_face) in y dir

    % Mesh grid for multigrid solver
    max_level = floor(log2(N(num_N)))-2;
    Iter_smooth = Iter_cycle.*ones(max_level,1);
    Iter_smooth(end) = Iter_last_cycle;

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
    dp_cc = zeros(size(X_cc));

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

        [u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,p_cc,dp_cc] = NS_solver_VK(u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,p_cc,dp_cc,dt,Delta_x,Delta_y,dx,dy,Nx,Ny,Re,u_top, ...
            Solver_burgers,max_iter,max_level,Nx_level,Ny_level,dx_level,dy_level,Iter_smooth,tol_PPE,tol_burgers);
        res_u = rms(u_cc - u_cc_old,"all");
        res_v = rms(v_cc - v_cc_old,"all");
        fprintf('Residual (u,v): (%e, %e) \n',res_u,res_v)

        if(res_u < tol_steady_state && res_v < tol_steady_state)
            break
        end

    end

%% Plot comparison of u,v
    Nx_half = floor(Nx/2) + 1;
    Ny_half = floor(Ny/2) + 1;
    
    u_half = u_cc(:,Nx_half);
    v_half = v_cc(Ny_half,:);
    x_half = X_cc(Ny_half,:);
    y_half = Y_cc(:,Nx_half);

    % u
    figure(1)
    if (num_N == 1)
        plot(data_u_lit(:,1),data_u_lit(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Ghia et al. (1982)'); hold on
        plot(data_u_lit(:,4),data_u_lit(:,5),'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','DisplayName','Nallasamy & Prasad (1977)'); hold on
        plot(data_u_lit(:,7),data_u_lit(:,8),'LineStyle','none','Marker','diamond','MarkerEdgeColor','b','MarkerFaceColor','b','DisplayName','Agarwal (1981)'); hold on
    end
    plot(u_half,y_half,'color',color_line(num_N),'DisplayName',['This study: \itN = \rm',num2str(N(num_N))]); hold on
    xlabel('$ u \ (m/s) $','Interpreter','latex')
    ylabel('$ y \ (m) $','Interpreter','latex')
    grid on
    legend('Location','southeast')
    ylim([0,1])
    title(['Re = ',num2str(Re),': \it{u} \rm{profile}'])
    drawnow

    % v
    figure(2)
    if (num_N == 1)
        plot(data_v_lit(:,1),data_v_lit(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Ghia et al. (1982)'); hold on
        plot(data_v_lit(:,4),data_v_lit(:,5),'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','DisplayName','Nallasamy & Prasad (1977)'); hold on
    end
    plot(x_half,v_half,'color',color_line(num_N),'DisplayName',['This study: \it N = \rm',num2str(N(num_N))]); hold on
    xlabel('$ x \ (m) $','Interpreter','latex')
    ylabel('$ v \ (m/s) $','Interpreter','latex')
    grid on
    legend('Location','southwest')
    xlim([0,1])
    title(['Re = ',num2str(Re),': \it{v} \rm{profile}'])
    drawnow
end