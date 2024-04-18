using .Threads, Serialization, PyCall, #BenchmarkTools

function depickle(max_M,n_dps)
    mpm = pyimport("mpmath")
    mpm.mp.dps = n_dps
    setprecision(BigFloat,n_dps;base=10)
    
    py"""
    import pickle
    def import_pickle_py(max_M,n_dps):
        with open(f'/home/c.delacruz/JuliaYuri/M{max_M}_p{n_dps}/zetazeros_M{max_M}_p{n_dps}_py.txt','rb') as file:
            zz = pickle.load(file)
        return zz
    """
    
    import_pickle = py"import_pickle_py"
    zz_py = import_pickle(max_M,n_dps)
    [parse(BigFloat,mpm.nstr(item,n_dps)) for item in zz_py]
end

function omega(vect::Array{BigFloat})
    n_rows = 2*length(vect)+1
    mat = zeros(BigFloat,(n_rows,n_rows))        # Either BigFloat
    # mat = matrix(RR,ones(n_rows,n_rows))      # or Nemo-Arb
    @threads for j in 1:n_rows
        for i in 1:length(vect)
            tmp = j^(-0.5-vect[i]*im)           # Either BigFloat
            # tmp = j^(-CC(0.5,vect[i]))        # or Nemo-Arb
            @inbounds mat[2*i-1,j]=real(tmp)
            @inbounds mat[2*i,j]=imag(tmp)
        end
    end
    mat[n_rows,n_rows]=1
    transpose(mat)
end

function gauss_elim!(input_mat,mult,addrow)
    @threads for j in axes(input_mat,2) 
        for i in axes(input_mat,1)         
            @inbounds input_mat[i,j] += mult[i] * addrow[j]
        end
    end
end

function delta_coef(omega_mat)
    max_N = size(omega_mat,1)

    # Stage 1: Gauss elimination without pivoting
    st_1 = time()
    for k in 1:max_N
        st = time()
        omega_mat[k+1:end, k] .= -omega_mat[k+1:end, k] / omega_mat[k, k]

        input_mat = @view omega_mat[k+1:end,k+1:end]
        mult = @view omega_mat[k+1:end,k]
        addrow = @view omega_mat[k,k+1:end]
        gauss_elim!(input_mat,mult,addrow)

        println("pGauss - Cycle k=$k: $(time() - st) s")
    end
    println("pGauss - Total: $(time() - st_1) s \n")

    # Stage 2: Product of the applied elementary matrices
    st_1 = time()
    for k in 1:max_N
        st = time()
        omega_mat[k, k] = 1
        for j in k+1:max_N
            @threads for i in j+1:max_N
                @inbounds omega_mat[i, k] += omega_mat[j, k] * omega_mat[i, j] 
            end
        end
        println("Product - Cycle k=$k: $(time() - st) s")
    end
    println("Product - Total: $(time() - st_1) s \n")


    # Stage 3: Computation of δ_{N,n}'s
    st_1 = time()
    delta = Dict{Int, Vector{BigFloat}}()
    @threads for Nn in 1:max_N
        delta[Nn] = omega_mat[Nn, 1:Nn] ./ omega_mat[Nn, 1]
    end
    println("Coeff. Computation: $(time() - st_1) s \n")

    delta
end

function main()
    max_M = 1_300
    n_dps = 10_000
    input_format = "jl"

    setprecision(BigFloat,n_dps;base=10)
    # set_precision!(Balls,Int(ceil(n_dps/log10(2))))             # Nemo-Arb

    # RR = RealField()                              # Nemo-Arb
    # CC = ComplexField()                           # Nemo-Arb

    st = time()
    if input_format == "jl"
        zz = open(deserialize, "/home/c.delacruz/JuliaYuri/M$(max_M)_p$(n_dps)/zetazeros_M$(max_M)_p$(n_dps)_jl.txt")
    elseif input_format == "py"
        zz = depickle(max_M,n_dps)
    end
    println("Importation of zeta zeros: $(time()-st) s. \n")

    omega_mat = omega(zz)
    println("Matrix set-up: $(time()-st) s. \n")

    delta = delta_coef(omega_mat)
    println("Total time: $(time()-st) s")

    open("/home/c.delacruz/JuliaYuri/M$(max_M)_p$(n_dps)/coef_delta_M$(max_M)_p$(n_dps).txt", "w") do f
        serialize(f,delta)
    end
end

main()
