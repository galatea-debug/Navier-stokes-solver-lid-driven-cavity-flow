function e_fine = prolong_bilinear(e_fine, e_coarse, Nx_coarse, Ny_coarse)
    % Prolongs e_coarse (Ny_coarse x Nx_coarse) to e_fine ((2*Ny_coarse-1) x (2*Nx_coarse-1))
    % using bilinear interpolation. Assumes e_fine is preallocated.

    for i = 1:Ny_coarse+1
        for j = 1:Nx_coarse+1
            i_fine = 2*i-1;
            j_fine = 2*j-1;
            e_fine(i_fine, j_fine) = e_coarse(i, j);                                        % Step 1: Set values at coarse grid points (odd i_fine, odd j_fine)
            e_fine(i_fine, j_fine + 1) = (e_coarse(i, j) + e_coarse(i, j+1)) / 2;           % Step 2: Set midpoints in x-direction (odd i_fine, even j_fine)
            e_fine(i_fine + 1, j_fine) = (e_coarse(i, j) + e_coarse(i+1, j)) / 2;           % Step 3: Set midpoints in y-direction (even i_fine, odd j_fine)
            e_fine(i_fine + 1, j_fine + 1) = (e_coarse(i, j) + e_coarse(i+1, j) + ...
                                      e_coarse(i, j+1) + e_coarse(i+1, j+1)) / 4;           % Step 4: Set center points (even i_fine, even j_fine)
        end
    end

end