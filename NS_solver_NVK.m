function [u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,p_cc] = NS_solver_NVK(u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,p_cc,dt,Delta_x,Delta_y,dx,dy,Nx,Ny,Re,u_top,Solver_burgers,...
    max_iter_burgers,max_iter,max_level,Nx_level,Ny_level,dx_level,dy_level,Iter_smooth,tol_PPE,tol_burgers)

    % Solve for intermediate velocity: u*
    [u_cc_star,v_cc_star] = Burgers_solver_NVK(u_cc,v_cc,U_fc,V_fc,u_cc_old,v_cc_old,U_fc_old,V_fc_old,dt, ...
                                        Delta_x,Delta_y,dx,dy,Nx,Ny,Re,u_top,Solver_burgers,tol_burgers,max_iter_burgers);
    
    [U_fc_star,~] = Get_face_variable(u_cc_star,Nx,Ny);       % Get intermidiate face velocity
    [~,V_fc_star] = Get_face_variable(v_cc_star,Nx,Ny);       % Get intermidiate face velocity

    % PPE
    Div_velocity_star = calculate_divergence(U_fc_star,V_fc_star,Delta_x,Delta_y,Nx,Ny);
    RHS_p = Div_velocity_star./dt;
    
    level = 1;
    Iter = 0;
    Iter_count = zeros(max_iter,1);
    RMS_residual = zeros(max_iter,1);
    
    fprintf('PPE solver Begins: \n')
    [p_cc,~,~] = Multigrid_V_cycle_N(max_iter,Iter,level,max_level,Nx_level,Ny_level, ...
          dx_level,dy_level,RHS_p,p_cc,RMS_residual,Iter_count,Iter_smooth,tol_PPE);        % Multigrid solver for PPE   
    fprintf('PPE solver Ends \n')

    % Update velocity with pressure
    [Grad_p_cc_x,Grad_p_cc_y] = calculate_gradient_p_cc(p_cc,dx,dy,Nx,Ny);
    [Grad_p_fc_x,Grad_p_fc_y] = calculate_gradient_p_fc(p_cc,dx,dy,Nx,Ny);

    u_cc_new = u_cc_star - dt.*Grad_p_cc_x;
    v_cc_new = v_cc_star - dt.*Grad_p_cc_y;
    U_fc_new = U_fc_star - dt.*Grad_p_fc_x;
    V_fc_new = V_fc_star - dt.*Grad_p_fc_y;

    % Update velocity in time
    u_cc_old = u_cc;
    v_cc_old = v_cc;
    U_fc_old = U_fc;
    V_fc_old = V_fc;

    u_cc = u_cc_new;
    v_cc = v_cc_new;
    U_fc = U_fc_new;
    V_fc = V_fc_new;

end