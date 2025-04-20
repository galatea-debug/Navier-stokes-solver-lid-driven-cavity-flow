function [Phi_x,Phi_y] = Get_face_variable(phi,Nx,Ny)

    Phi_x = zeros(Ny+1,Nx+1);
    Phi_y = zeros(Ny+1,Nx+1);
    
    for i = 1:Nx+1
        for j = 1:Ny+1
            Phi_x(j,i) = 0.5*(phi(j,i) + phi(j,i+1));
            Phi_y(j,i) = 0.5*(phi(j,i) + phi(j+1,i));
        end
    end

end