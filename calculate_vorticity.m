function w = calculate_vorticity(u,v,Nx,Ny,dx,dy)
    w = zeros(size(u));
    
    % Inner points: Central difference
    for i = 2:Nx+1
        for j = 2:Ny+1
            dv_i = v(j,i+1) - v(j,i-1);
            du_j = u(j+1,i) - u(j-1,i);
            w(j,i) = dv_i/dx/2 - du_j/dy/2;
        end
    end

    % Left BC: 1st order Upwind
    for k = 2:Ny+1
        dv_i = v(k,2) - v(k,1);
        du_j = u(k+1,1) - u(k,1);
        w(k,1) = dv_i/dx - du_j/dy;
    end

    % Right BC: 1st order Upwind
    for k = 2:Ny+1
        dv_i = v(k,end) - v(k,end - 1);
        du_j = u(k+1,i) - u(k,i);
        w(k,end) = dv_i/dx - du_j/dy;
    end

    % Upper BC: 1st order Upwind
    for k = 2:Nx+1
        dv_i = v(2,k) - v(1,k);
        du_j = u(2,k+1) - u(1,k);
        w(1,k) = dv_i/dx - du_j/dy;
    end

    % Lower BC: 1st order Upwind
    for k = 2:Nx+1
        dv_i = v(end,k) - v(end-1,k);
        du_j = u(end,k+1) - u(end-1,k);
        w(end,k) = dv_i/dx - du_j/dy;
    end

end