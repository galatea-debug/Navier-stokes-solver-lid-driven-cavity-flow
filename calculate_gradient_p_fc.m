function [Grad_p_fc_x,Grad_p_fc_y] = calculate_gradient_p_fc(p_cc,dx,dy,Nx,Ny)

Grad_p_fc_x = zeros(Ny+1,Nx+1);
Grad_p_fc_y = zeros(Ny+1,Nx+1);

for i = 1:Nx+1
    for j = 1:Ny+1
        Grad_p_fc_x(j,i) = (p_cc(j,i+1) - p_cc(j,i))/dx;
        Grad_p_fc_y(j,i) = (p_cc(j+1,i) - p_cc(j,i))/dy;
    end
end

end