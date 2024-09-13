export 
    pack_poly_2d,
    simulation_2d

    
function pack_poly_2d(N, K, D, G, M, P_target, W_factor, seed, plotit)
    mat"""
    addpath('src/matlab_functions/')
    N =double($(N))
    K = double($(K))
    M = double($(M))
    D = double($(D))
    G = double($(G))
    P_target = double($(P_target))
    W_factor = double($(W_factor))
    seed = double($(seed))
    plotit = double($(plotit))
    packing_poly_2d(N, K, D, G, M, P_target, W_factor, seed, plotit)
    """
end

function simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    mat"""
    addpath('src/matlab_functions/')
    K = double($(K))
    M = double($(M))
    Bv = double($(Bv))
    w_D = double($(w_D))
    N =double($(N))
    P = double($(P))
    W = double($(W))
    seed = double($(seed))
    simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    """
end