function [Diff_phi_x,Diff_phi_y] = Calculate_diffusion(phi,Nx,Ny,dx,dy)

Diff_phi_x = zeros(Ny+2,Nx+2);
Diff_phi_y = zeros(Ny+2,Nx+2);

for i = 2:Nx+1
    for j = 2:Ny+1
        Diff_phi_x(j,i) = (phi(j,i+1) - 2*phi(j,i) + phi(j,i-1))/dx^2;
        Diff_phi_y(j,i) = (phi(j+1,i) - 2*phi(j,i) + phi(j-1,i))/dy^2;
    end
end


end