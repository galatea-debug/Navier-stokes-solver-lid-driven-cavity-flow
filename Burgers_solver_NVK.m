function [u_star,v_star] = Burgers_solver_NVK(u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,dt, ...
                                        Delta_x,Delta_y,dx,dy,Nx,Ny,Re,u_top,solver,tol_burgers,max_iter_burgers)

%% Calculate RHS
fprintf('Burgers solver Begins: \n')
[u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_top);      % Add BCs

% Calculate convection term: C^n, C^n-1
C_u_n = Calculate_convection(U_fc,V_fc,u_cc,Nx,Ny,Delta_x,Delta_y);
C_v_n = Calculate_convection(U_fc,V_fc,v_cc,Nx,Ny,Delta_x,Delta_y);
C_u_n_1 = Calculate_convection(U_fc_old,V_fc_old,u_cc_old,Nx,Ny,Delta_x,Delta_y);
C_v_n_1 = Calculate_convection(U_fc_old,V_fc_old,v_cc_old,Nx,Ny,Delta_x,Delta_y);

C_u = 1.5.*C_u_n - 0.5.*C_u_n_1;
C_v = 1.5.*C_v_n - 0.5.*C_v_n_1;

% Calculate diffusion term: D^n
[Diff_u_x,Diff_u_y] = Calculate_diffusion(u_cc,Nx,Ny,dx,dy);
[Diff_v_x,Diff_v_y] = Calculate_diffusion(v_cc,Nx,Ny,dx,dy);

%% solver for burgers
switch solver
    case "ADI"
        % Assemble RHS for ADI
        % RHS_u = - dt.*C_u + dt/(2*Re).*(Diff_u_x + Diff_u_y);
        % RHS_v = - dt.*C_v + dt/(2*Re).*(Diff_v_x + Diff_v_y);

        RHS_conv_u = -dt * C_u;  % Total convection contribution over full step
        RHS_conv_v = -dt * C_v;

        u_star = ADI_solver(u_cc, Nx, Ny, dt, dx, dy, Re, RHS_conv_u, u_top, true);
        v_star = ADI_solver(v_cc, Nx, Ny, dt, dx, dy, Re, RHS_conv_v, u_top, false);

    case "GS"
        % Assemble RHS for GS
        RHS_u = u_cc - dt.*C_u + dt/(2*Re).*(Diff_u_x + Diff_u_y);
        RHS_v = v_cc - dt.*C_v + dt/(2*Re).*(Diff_v_x + Diff_v_y);

        for k = 1:max_iter_burgers
            u_cc = GS_solver_burgers(u_cc,Nx,Ny,dx,dy,dt,Re,RHS_u);
            v_cc = GS_solver_burgers(v_cc,Nx,Ny,dx,dy,dt,Re,RHS_v);
            [u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_top);

            residual_u = calculate_residual_burgers(u_cc,Nx,Ny,dx,dy,dt,Re,RHS_u);
            residual_v = calculate_residual_burgers(v_cc,Nx,Ny,dx,dy,dt,Re,RHS_v);

            RMS_residual_u = rms(residual_u,"all");
            RMS_residual_v = rms(residual_v,"all");
            fprintf('       Iteration number: %i --- RMS residual (u,v): (%e, %e)\n',k,RMS_residual_u,RMS_residual_v)

            if (RMS_residual_u < tol_burgers && RMS_residual_v < tol_burgers)
                break
            end
        end

        fprintf('Burgers solver Ends \n')
        u_star = u_cc;
        v_star = v_cc;
end

end