using .Threads, JSON #, BenchmarkTools

# read_zz - reads the imaginary parts of the first (max_z) 
# Riemann zeta zeros from the file (filepath)
function read_zz(filepath,max_z)
    zz = Array{BigFloat}[]
    open(filepath) do f
        for _ in 1:max_z
            line = readline(f)
            numb = parse(BigFloat,line)
            push!(zz,numb)
        end
    end
    zz
end


# addrow! - adds to each row in (input_mat) the row
# (mult * addrow). Parallelized. 
function addrow!(input_mat,mult,addrow)
    @threads for j in axes(input_mat,2) 
        for i in axes(input_mat,1)         
            @inbounds input_mat[i,j] += mult[i] * addrow[j]
        end
    end
end


# omega - sets up the "Vandermonde matrix" Ω_M from the 
# vector (vect) that consists of the imaginary parts of zeta 
# zeros. Here M = 2*length(vec) the output is a matrix of 
# dimensions N x N, where N=2M+1.
function omega(vect::Array{BigFloat})
    n_rows = 2*length(vect)+1
    mat = zeros(BigFloat,(n_rows,n_rows))
    @threads for j in 1:n_rows
        for i in 1:length(vect)
            tmp = j^(-0.5-vect[i]*im)
            @inbounds mat[2*i-1,j]=real(tmp)
            @inbounds mat[2*i,j]=imag(tmp)
        end
    end
    mat[n_rows,n_rows]=1
    transpose(mat)
end


# gauss_elim! - performs Gauss elimination on the matrix (omega_mat)
# while recording the multipliers under the main diagonal.
function gauss_elim!(omega_mat::Matrix{BigFloat})
    max_N = size(omega_mat,1)
    st_1 = time()
    for k in 1:max_N
        st = time()
        omega_mat[k+1:end, k] .= -omega_mat[k+1:end, k] / omega_mat[k, k]

        input_mat = @view omega_mat[k+1:end,k+1:end]
        mult = @view omega_mat[k+1:end,k]
        addrow = @view omega_mat[k,k+1:end]
        addrow!(input_mat,mult,addrow)

        println("pGauss - Cycle k=$k: $(time() - st) s")
    end
    println("pGauss - Total: $(time() - st_1) s \n")
end


# mult_elementary! - multiplies the elementary matrices built from
# the multipliers recorded by (gauss_elim!) under the main diagonal 
# of (omega_mat). After this function, the coefficients δ_{N,n}
# before normalization are contained under the main diagonal. 
function mult_elementary!(omega_mat::Matrix{BigFloat})
    max_N = size(omega_mat,1)
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
end


# normalized_coefficients - From the processed (omega_mat) matrix 
# produces a dictionary (delta) with integer keys N that range from 1 
# to max_N = size of (omega_mat). The entry (delta[N]) is a list of the 
# coefficients δ_{N,n}, normalized so that δ_{N,1} = 1. 
function normalized_coefficients(omega_mat::Matrix{BigFloat})
    st_1 = time()
    delta = Dict{Int, Vector{BigFloat}}()
    @threads for Mm in 1:max_M
        Nn = 2*Mm+1
        delta[Mm] = omega_mat[Nn, 1:Nn] ./ omega_mat[Nn, 1]
    end
    println("Coeff. Computation: $(time() - st_1) s \n")

    delta
end


# delta_coef - produces the coefficients δ_{N,n} from the Vandermonde
# matrix (omega_mat)
function delta_coef(omega_mat)
    gauss_elim!(omega_mat)
    mult_elementary!(omega_mat)
    normalized_coefficients(omega_mat)
end


# cuts - determines the partition of coefficients from M=1:(max_M) into
# chunks of sizes of (chunk_size) assuming the size of one coefficient is 
# (typical_size)
function cuts(max_M,typical_size,chunk_size)
    cutlist = [0]
    new_cut = 0

    while new_cut < max_M
        gsum = cutlist[end]*(cutlist[end]+2)
        disc = 4*(1+chunk_size/typical_size+gsum)
        new_cut = floor(0.5 * (sqrt(disc) - 2))
        push!(cutlist,new_cut)
        println(new_cut)
    end

    if cutlist[end] >= max_M
        cutlist[end] = max_M
    else
        push!(cutlist,max_M)
    end

    cutlist
end


# write_coefs - write the coefficients (delta) into chunk files 
# according to the partition in (cutlist)
function write_coefs(delta,cutlist)
    for i in 1:length(cutlist)-1
        chunk = Dict(i => delta[i] for i in cutlist[i]+1:cutlist[i])
        chunk_json = JSON.json(chunk)
        open("../Data/p{n_precision}/CoefDelta_M$(cutlist[i]+1)-$(cutlist[i+1])_p$(n_dps).txt", "w") do f
            write(f, chunk_json)
        end
    end
end


# Here we go!
function main()
    max_M = 1_300
    max_computed_M = 1_300
    n_dps = 10_000
    chunk_size = Int(3e9)

    filepath = "../Data/p$(n_dps)/NImZetaZero_M$(max_computed_M)_p$(n_dps).txt"
    setprecision(BigFloat,n_dps;base=10)
    
    # Read zeta zeros from file and store in array
    zz = read_zz(filepath,max_z)
    println("Importation of zeta zeros: $(time()-st) s. \n")

    # Set up Vandermonde matrix
    omega_mat = omega(zz)
    println("Matrix set-up: $(time()-st) s. \n")

    # Compute coefficients
    delta = delta_coef(omega_mat)
    println("Total time: $(time()-st) s")

    # Determine partition of files and write files
    typical_size = sizeof(zz[1])
    cutlist = cuts(max_M,typical_size,chunk_size)
    write_coefs(delta,cutlist)
end


main()