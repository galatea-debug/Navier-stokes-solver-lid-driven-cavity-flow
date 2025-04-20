function [phi,Iter_count,RMS_residual] = Multigrid_V_cycle_N(max_iter,Iter,level,max_level,Nx,Ny, ...
    dx,dy,RHS,phi,RMS_residual,Iter_count,Iter_smooth,tol)

res = cell(1,max_level);
error_restrict = cell(1,max_level);
error_prolong = cell(1,max_level);

for i = 1:max_iter
    Iter = Iter + 1;
    % Pre-smoothing
    for k1 = 1:Iter_smooth(level)
        phi = GS_solver(phi,Nx(level),Ny(level),dx(level),dy(level),RHS);
        % phi = ADI_solver(phi,Nx(level),Ny(level),dx(level),dy(level),f);
        phi = phi - mean(phi,'all');
        phi = Implement_Neumann_BCs(phi,level); % Implement Neumann BCs

    end

    % Get residual
    res{1} = calculate_residual(phi,Nx(level),Ny(level),dx(level),dy(level),RHS);
    error_restrict{1} = zeros(Ny(level)+2,Nx(level)+2);

    % Resitriction
    for level = 2:max_level
        % Restrict residual
        res_coarse = zeros(Ny(level)+2,Nx(level)+2);
        res_coarse = restrict_full_weighting(res{level-1},res_coarse,Nx(level),Ny(level));

        % Calculate error at corase grid
        e_coarse = zeros(Ny(level)+2,Nx(level)+2);

        for k1 = 1:Iter_smooth(level)
            e_coarse = GS_solver(e_coarse,Nx(level),Ny(level),dx(level),dy(level),-res_coarse);
            % e_coarse = ADI_solver(e_coarse,Nx(level),Ny(level),dx(level),dy(level),res_coarse);
            e_coarse = e_coarse - mean(e_coarse,'all');
            e_coarse = Implement_Neumann_BCs(e_coarse,level); % Implement Neumann BCs
        end
        error_restrict{level} = e_coarse;

        res{level} = calculate_residual(e_coarse,Nx(level),Ny(level),dx(level),dy(level),-res_coarse);
    end

    % Prolongation
    error_prolong{level} = e_coarse;
    for level = max_level:-1:2
        e_fine = zeros(Ny(level-1)+2, Nx(level-1)+2);
        e_fine = prolong_bilinear(e_fine, error_prolong{level}, Nx(level), Ny(level));
        e_fine = error_restrict{level-1} + e_fine;

        % Post-smoothing on error
        for k1 = 1:Iter_smooth(level)
            e_fine = GS_solver(e_fine,Nx(level-1),Ny(level-1),dx(level-1),dy(level-1),-res{level-1});
            % e_fine = ADI_solver(e_fine,Nx(level-1),Ny(level-1),dx(level-1),dy(level-1),res{level-1});
            e_fine = e_fine - mean(e_fine,"all");
            e_fine = Implement_Neumann_BCs(e_fine,level); % Implement Neumann BCs
        end
        error_prolong{level - 1} = e_fine;

    end

    % Add error back
    level = level - 1;
    phi = phi + error_prolong{level};

    % Final Post-smoothing on u
    for k1 = 1:Iter_smooth(level)
        phi = GS_solver(phi,Nx(level),Ny(level),dx(level),dy(level),RHS);
        % phi = ADI_solver(phi,Nx(level),Ny(level),dx(level),dy(level),f);
        phi = phi - mean(phi,'all');
        phi = Implement_Neumann_BCs(phi,level); % Implement Neumann BCs
    end

    % phi = phi - mean(phi,'all');
    Residual = calculate_residual(phi,Nx(level),Ny(level),dx(level),dy(level),RHS);

    RMS_residual(i) = rms(Residual,"all");
    Iter_count(i) = Iter;
    fprintf('       Iteration number: %i --- RMS residual: %e\n',Iter,RMS_residual(i))

    if (RMS_residual(i) < tol)
        break
    end
end

Iter_count = Iter_count(1:Iter);
RMS_residual = RMS_residual(1:Iter);

end