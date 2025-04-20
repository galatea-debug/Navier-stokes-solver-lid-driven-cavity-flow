function Div_velocity = calculate_divergence(U_fc,V_fc,Delta_x,Delta_y,Nx,Ny)

Div_velocity = zeros(Ny+2,Nx+2);

for i = 2:Nx+1
    for j = 2:Ny+1
        Div_velocity(j,i) = (U_fc(j,i) - U_fc(j,i-1))/Delta_x + ...
                            (V_fc(j,i) - V_fc(j-1,i))/Delta_y;
    end
end

end