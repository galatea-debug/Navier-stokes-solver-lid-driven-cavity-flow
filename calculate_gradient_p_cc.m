function [Grad_p_cc_x,Grad_p_cc_y] = calculate_gradient_p_cc(p_cc,dx,dy,Nx,Ny)

Grad_p_cc_x = zeros(Ny+2,Nx+2);
Grad_p_cc_y = zeros(Ny+2,Nx+2);

for i = 2:Nx+1
    for j = 2:Ny+1
        Grad_p_cc_x(j,i) = (p_cc(j,i+1) - p_cc(j,i-1))/(2*dx);
        Grad_p_cc_y(j,i) = (p_cc(j+1,i) - p_cc(j-1,i))/(2*dy);
    end
end

end