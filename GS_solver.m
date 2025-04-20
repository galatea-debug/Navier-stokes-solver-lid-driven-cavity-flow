function phi = GS_solver(phi,Nx,Ny,dx,dy,RHS)

    w = 1.0;
    a_P = 2*(1/dx^2 + 1/dy^2);
    a_E = 1/dx^2;
    a_W = 1/dx^2;
    a_N = 1/dy^2;
    a_S = 1/dy^2;
    a_Pinv=1/a_P;
    for i = 2:Nx+1
        for j = 2:Ny+1
            phit=a_Pinv*(a_W*phi(j,i-1) + a_E*phi(j,i+1) + a_N*phi(j+1,i) + a_S*phi(j-1,i) - RHS(j,i));
            phi(j,i) = (1-w)*phi(j,i) + w*phit;
        end
    end
end