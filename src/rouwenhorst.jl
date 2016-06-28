function rouwenhorst(N::Int, μ::Real, ρ::Real, σ::Real)
    """
    Use Rouwenhorst's method to produce finite state Markov
    approximation of the AR(1) processes
    y_t = μ*(1-ρ) + ρ y_{t-1} + ε_t,
    where ε_t ~ N (0, σ^2)
    Parameters
    ----------
    N : int
        Number of points in markov process
    ρ : float
        Persistence parameter in AR(1) process
    σ : float
        Standard deviation of random component of AR(1) process
    μ : float, optional(default=0.0)
        Mean of AR(1) process
    Returns
    -------
    y : array(dtype=float, ndim=1)
        1d-Array of nodes in the state space
    Θ : array(dtype=float, ndim=2)
        Matrix transition probabilities for Markov Process
    """
    σ_y = σ / sqrt(1-ρ^2)
    p  = (1+ρ)/2
    Θ = [p 1-p; 1-p p]

    for n = 3:N
        z_vec = zeros(n-1,1)
        z_vec_long = zeros(1, n)
        Θ = p.*[Θ z_vec; z_vec_long] +
            (1-p).*[z_vec Θ; z_vec_long] +
            (1-p).*[z_vec_long; Θ z_vec] +
            p.*[z_vec_long; z_vec Θ]
        Θ[2:end-1,:] ./=  2.0
    end

    ψ = sqrt(N-1) * σ_y
    w = linspace(-ψ, ψ, N)

    w .+= μ   # center process around its mean (wbar / (1 - rho))

    return w, Θ
end