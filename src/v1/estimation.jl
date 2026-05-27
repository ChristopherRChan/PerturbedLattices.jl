"""
    Estimation methods for perturbed lattice models. First case, for Framework 1 with gaussian perturbations and Strauss potential.
"""


function create_int_estimation_grid(B::Matrix{Float64}, Nest::Int64)
    x_grid = range(B[1,1], B[1,2], length=Nest)
    y_grid = range(B[2,1], B[2,2], length=Nest)
    grid = [[x, y] for x in x_grid for y in y_grid]
    return grid
end

function vec_loc_en_W(pl::PerturbedLatticeV1, points_D_n::Vector{Vector{Float64}}, grid_D_n::Vector{Vector{Float64}}, points_W_n_boundary::Vector{Vector{Float64}}, W_n::Matrix{Float64}, Nest::Int64, B::Matrix{Float64})
    n_points_D = length(points_D_n)
    vec_loc_en_W = zeros(n_points_D, Nest^2)
    grid = create_int_estimation_grid(B, Nest)
    config = [points_D_n; points_W_n_boundary]

    # Calculer la matrice d'adjacence pour la configuration initiale
    config_adjacency = compute_adjacency_matrix(config, pl.RS, 2)

    # Pour chaque point dans D_n
    for i in 1:n_points_D
        # Créer une copie de l'adjacency pour ce point
        temp_adjacency = copy(config_adjacency)

        # Pour chaque position testée
        for j in 1:length(grid)
            new_point = grid_D_n[i] .+ grid[j]

            # Mettre à jour la matrice d'adjacence pour la nouvelle position
            refresh_adjacency_matrix(temp_adjacency, config, pl.RS, new_point, i, 2)

            # Calculer l'énergie locale du point i en utilisant local_energy
            loc_en = local_energy(temp_adjacency, i)

            # Stocker l'énergie
            if W_n[1,1] + pl.RS <= new_point[1]  <= W_n[1,2] - pl.RS && W_n[2,1] + pl.RS<= new_point[2] <= W_n[2,2] - pl.RS
                vec_loc_en_W[i,j] = loc_en
            end
        end
    end

    return vec_loc_en_W, config_adjacency
end

function DLR_W(pl::PerturbedLatticeV1, points_D_n::Vector{Vector{Float64}}, grid_D_n::Vector{Vector{Float64}}, points_W_n_boundary::Vector{Vector{Float64}}, beta::Float64, theta::Float64, vec_loc_en_W::Matrix{Float64}, config_adjacency::Matrix{Int}, N_est::Int64, B::Matrix{Float64}, W_n::Matrix{Float64})
    DLR_2 = 0.0
    DLR_1 = 0.0
    N_points = length(points_D_n)
    
    grid_est = create_int_estimation_grid(B, N_est)

    for i in 1:N_points
        num_1 = 0.0
        num_2 = 0.0
        denom = 0.0
        X_i = points_D_n[i] .- grid_D_n[i]
        for j in 1:N_est^2
            x_j = grid_est[j]
            if W_n[1,1] + pl.RS <= grid_D_n[i][1] + x_j[1] <= W_n[1,2] - pl.RS && W_n[2,1] + pl.RS <= grid_D_n[i][2] + x_j[2] <= W_n[2,2] - pl.RS
                exp_term = exp(-theta*(x_j[1]^2 + x_j[2]^2) / 2)
                num_1 += (x_j[1]^2 + x_j[2]^2)/2*exp(-beta*vec_loc_en_W[i,j]) * exp_term 
                num_2 += vec_loc_en_W[i,j] * exp(-beta*vec_loc_en_W[i,j]) * exp_term
                denom += exp(-beta*vec_loc_en_W[i,j]) * exp_term
            end
        end
        DLR_2 += (local_energy(config_adjacency, i) - num_2/denom)
        DLR_1 += ( (X_i[1]^2 + X_i[2]^2)/2- num_1/denom)
    end
    return DLR_1/N_points, DLR_2 / N_points
end


function fit(pl::PerturbedLatticeV1, beta_init::Float64, theta_init::Float64, points_D_n::Vector{Vector{Float64}}, grid_D_n::Vector{Vector{Float64}}, points_W_n_boundary::Vector{Vector{Float64}}, W_n::Matrix{Float64}, N_est::Int64, B::Matrix{Float64})
    vec_loc_en_W, config_adjacency = vec_loc_en_W(pl, points_D_n, grid_D_n, points_W_n_boundary, W_n, N_est, B)

    # Fonction objectif : DLR²
    function objective(params)
        beta, theta = params[1], params[2]

        # Calculer DLR
        dlr_value_1, dlr_value_2 = DLR_W(pl, points_D_n, grid_D_n, points_W_n_boundary, beta, theta, vec_loc_en_W, config_adjacency, N_est, B, W_n)
        dlr_squared = dlr_value_1^2 + dlr_value_2^2
        return dlr_squared
    end

    # Optimisation sans contraintes
    result = optimize(objective, [beta_init, theta_init], NelderMead())

    minimizer = Optim.minimizer(result)
    beta_opt = minimizer[1]
    theta_opt = minimizer[2]

    dlr_squared_opt = Optim.minimum(result)
    dlr_opt = sqrt(dlr_squared_opt)

    return beta_opt, theta_opt, dlr_opt, vec_loc_en_W, config_adjacency
end