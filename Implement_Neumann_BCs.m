function phi = Implement_Neumann_BCs(phi,level)

    phi(1, :) = phi(2,:); % u(0, y) boundary
    phi(end, :) = phi(end-1,:); % u(1, y) boundary
    phi(:, 1) =  phi(:,2); % u(x, 0) boundary
    phi(:, end) = phi(:,end-1); % u(x, 1) boundary

end