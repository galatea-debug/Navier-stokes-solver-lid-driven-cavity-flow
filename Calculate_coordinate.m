function [X_cc,Y_cc,X_u_fc,Y_u_fc,X_v_fc,Y_v_fc,Delta_x,Delta_y,dx,dy] = Calculate_coordinate(x,y,Nx,Ny)

    %% Calculate coordinate
    % Initialize matrices
    X_cc = zeros(Ny+2,Nx+2);
    Y_cc = zeros(Ny+2,Nx+2);
    X_u_fc = zeros(Ny+1,Nx+1);
    Y_u_fc = zeros(Ny+1,Nx+1);
    X_v_fc = zeros(Ny+1,Nx+1);
    Y_v_fc = zeros(Ny+1,Nx+1);
    
    % Find coordinate for cell center
    for i = 1:Nx+2
        X_cc(:,i) = (x(i) + x(i+1))/2;
    end
    
    for j = 1:Ny+2
        Y_cc(j,:) = (y(j) + y(j+1))/2;
    end
    
    % Find coordinate for face center
    for i = 2:Nx+2
        X_u_fc(:,i-1) = x(i);
        X_v_fc(:,i-1) = (x(i) + x(i-1))/2;
    end
    
    for j = 2:Ny+2
        Y_u_fc(j-1,:) = (y(j) + y(j-1))/2;
        Y_v_fc(j-1,:) = y(j);
    end
    
    %% Calculate grid size
    % Initialize matrices
    Delta_x = x(2) - x(1);         % Distance between edge point in x dir
    Delta_y = y(2) - y(1);         % Distance between edge point in y dir
    dx = X_cc(1,2) - X_cc(1,1);              % Distance between face center point (v_face) in x dir
    dy = Y_cc(2,1) - Y_cc(1,1);              % Distance between face center point (u_face) in y dir

    % Delta_x = zeros(Ny+2,Nx+2);         % Distance between edge point in x dir
    % Delta_y = zeros(Ny+2,Nx+2);         % Distance between edge point in y dir
    % dx = zeros(Ny+1,Nx+1);              % Distance between face center point (v_face) in x dir
    % dy = zeros(Ny+1,Nx+1);              % Distance between face center point (u_face) in y dir
    
    % % Delta x
    % for i = 1:Nx+2
    %     Delta_x(:,i) = (x(i+1) - x(i))/2;
    % end
    % 
    % % Delta y
    % for j = 1:Nx+2
    %     Delta_y(j,:) = (y(j+1) - y(j))/2;
    % end
    % 
    % % dx
    % for i = 1:Nx+1
    %     dx(:,i) = (X_cc(1,i+1) - X_cc(1,i))/2;
    % end
    % 
    % % dy
    % for j = 1:Nx+1
    %     dy(j,:) = (Y_cc(j+1,1) - Y_cc(j,1))/2;
    % end

end