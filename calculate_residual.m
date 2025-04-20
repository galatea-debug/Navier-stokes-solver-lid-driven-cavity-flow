function residual = calculate_residual(phi,Nx,Ny,dx,dy,RHS)
    residual = zeros(size(RHS));
    a_P = 2*(1/dx^2 + 1/dy^2);
    a_E = 1/dx^2;
    a_W = 1/dx^2;
    a_N = 1/dy^2;
    a_S = 1/dy^2;

    for i = 2:Nx+1
        for j = 2:Ny+1
            residual(j,i) = a_E*phi(j,i+1) + a_W*phi(j,i-1) + a_N*phi(j+1,i) + a_S*phi(j-1,i) -a_P*phi(j,i) - RHS(j,i);
        end
    end

end