function phi = GS_solver_burgers(phi,Nx,Ny,dx,dy,dt,Re,RHS)

    rx = dt/(2*Re*dx^2);
    ry = dt/(2*Re*dy^2);

    a_P = 1+2*rx+2*ry;
    a_N = ry;
    a_S = ry;
    a_E = rx;
    a_W = rx;
    a_P_inv = 1/a_P;

    for i = 2:Nx+1
        for j = 2:Ny+1
            phi(j,i) = a_P_inv*(a_N*phi(j+1,i) + a_S*phi(j-1,i) + a_E*phi(j,i+1) + a_W*phi(j,i-1) + RHS(j,i));
        end
    end

end