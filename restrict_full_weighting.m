function res_coarse = restrict_full_weighting(res_fine,res_coarse,Nx_coarse,Ny_coarse)

for i = 2:Ny_coarse+1
    for j = 2:Nx_coarse+1
        I = 2 * i - 1; % Fine grid i-index
        J = 2 * j - 1; % Fine grid j-index
        
        % Full weighting stencil for 2D (1/16 * [1 2 1; 2 4 2; 1 2 1])
        res_coarse(i, j) = (1/16) * (...
            1 * res_fine(I-1, J-1) + 2 * res_fine(I-1, J) + 1 * res_fine(I-1, J+1) + ...
            2 * res_fine(I,   J-1) + 4 * res_fine(I,   J) + 2 * res_fine(I,   J+1) + ...
            1 * res_fine(I+1, J-1) + 2 * res_fine(I+1, J) + 1 * res_fine(I+1, J+1));

        % res_coarse(i, j) = res_fine(I,J);
    end
end

end