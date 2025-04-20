function H = Calculate_convection(U_fc,V_fc,phi,Nx,Ny,Delta_x,Delta_y)

    H = zeros(Ny+2,Nx+2);
    
    for i = 2:Nx+1
        for j = 2:Ny+1
            H(j,i) = (U_fc(j,i)*(phi(j,i)+phi(j,i+1))*0.5 - U_fc(j,i-1)*(phi(j,i)+phi(j,i-1))*0.5)/Delta_x + ...
                     (V_fc(j,i)*(phi(j,i)+phi(j+1,i))*0.5 - V_fc(j-1,i)*(phi(j,i)+phi(j-1,i))*0.5)/Delta_y; %*Delta_yinv instead of divide and so on
        end
    end

end