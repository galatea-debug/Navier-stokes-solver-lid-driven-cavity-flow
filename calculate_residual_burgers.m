function residual = calculate_residual_burgers(phi,Nx,Ny,dx,dy,dt,Re,RHS)

    residual = zeros(Ny+2,Nx+2);
    rx = dt/(2*Re*dx^2);
    ry = dt/(2*Re*dy^2);

    a_P = 1+2*rx+2*ry;
    a_N = ry;
    a_S = ry;
    a_E = rx;
    a_W = rx;
    
    for i = 2:Nx+1
        for j = 2:Ny+1
            residual(j,i) = RHS(j,i) + a_N*phi(j+1,i) + a_S*phi(j-1,i) + a_E*phi(j,i+1) + a_W*phi(j,i-1) - a_P*phi(j,i);
        end
    end

end