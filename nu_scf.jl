##### COMPUTATION OF nu(M,K,s) #####


using FLoops, Nemo, JSON #, BenchmarkTools

# interval - outputs a tuple of numbers (X, Y, Z) from a string 
# (filename) of the form "CoefDelta_MX-Y_pZ.*".
function interval(filename::String)
    str_match = match(r"CoefDelta_M(\d+)-(\d+)_p(\d+)",filename)
    if !isnothing(str_match)
        map(x->parse(Int64,x), str_match.captures)
    end
end


# interval_list - applies the function (interval) to all the files
# in the folder "../RiemannZeta_Data/p$(dps)/", filters the tuples (X,Y,Z) with
# Z = dps, and outputs a list of the corresponding pairs (X,Y).
function interval_list(dps::Int)
    rawlist = Vector{Vector{Int64}}(filter(!isnothing,map(interval,readdir("Data/p$(dps)/"))))
    map(x -> (x[1],x[2]),filter(x -> x[3] == dps, rawlist))
end


# get_chain - for a list of tuples (X,Y) corresponding to the 
# starting and end point of integral intervals, outputs a sequence
# of tuples (X_1,Y_1), ..., (X_k,Y_k) such that 
# 1) m belongs to (X_1,Y_1) and M belongs to (X_k,Y_k)
# 2) Y_i + 1 = X_{i+1} for all i in {1,...,k-1}.
# If no such chain is found, it returns [].
# This function is used to select the CoefDelta files that need 
# to be opened to perform computations in a range with initial point
# m and end point M. 
function get_chain(intervals::Vector{Tuple{Int,Int}}, m::Int, M::Int)
    sorted_intervals = sort(intervals)
    
    i_m = findfirst(x -> (x[1] <= m <= x[2]), sorted_intervals)
    i_M = findlast(x -> (x[1] <= M <= x[2]), sorted_intervals)
    
    if isnothing(i_m) || isnothing(i_M)
        []
    else 
        chain = [sorted_intervals[i_m]]
        i = i_m + 1
        while i <= i_M
            k = findfirst(x -> x[1] == chain[end][2] + 1, sorted_intervals[i:i_M])
            if isnothing(k)
                if chain[end][1] == sorted_intervals[i][1]
                    chain[end] = sorted_intervals[i]
                    i += 1
                else
                    pop!(chain)
                    if chain == []
                        break
                    else 
                        i = findlast(x -> x == chain[end],sorted_intervals) + 1
                        chain[end] = sorted_intervals[i]
                    end
                end
            else
                push!(chain,sorted_intervals[i-1+k])
                if chain[end][1] <= M <= chain[end][2]
                    break
                end
                i += k 
            end
        end
        chain
    end
end


# open_coefdelta - opens the files CoefDelta corresponding to the 
# range (range_M) with a precision of (dps), and outputs a dictionary.
function open_coefdelta(range_M::AbstractRange{Int},dps::Int)
    intervals = interval_list(dps)
    chain = get_chain(intervals,range_M[1],range_M[end])

    data = Dict((Mlo,Mup) => Dict{String,Vector{String}}() for (Mlo,Mup) in chain)
    delta_pre = Dict((Mlo,Mup) => Dict{Int,Vector{BigFloat}}() for (Mlo,Mup) in chain)
    @floop for (Mlo,Mup) in chain
        data[Mlo,Mup] = JSON.parse(read("Data/p$(dps)/CoefDelta_M$(Mlo)-$(Mup)_p$(dps).json",String))
        delta_pre[Mlo,Mup] = Dict(parse(Int,key) => map(v -> parse(BigFloat,v),value) for (key,value) in data[Mlo,Mup])
    end
    merge(values(delta_pre)...)
end


# divis - returns the list of all divisors of an integer (n).
function divis(n::T) where {T<:Integer}
    f = factor(n)
    d = T[one(T)]
    for (p,m) in f
        c = T[p^i for i in 0:m]
        d = d*c'
        d = reshape(d,length(d))
    end
    sort!(d)
end


# power_23adic - returns the s-th power of x. If s is a 
# rational number with denominator whose prime factors are
# only 2 and 3, it performs the operation using integer powers
# and the built-in functions sqrt and cbrt. This implementation
# accelerates significantly computations when x is a BigFloat
# with a high precision.
function power_23adic(x::T,s::Number) where {T<:AbstractFloat}
    pow = x
    if eltype(s) == Rational{Int}
        while denominator(s) % 2 == 0
            pow = sqrt(pow)
            s *= 2
        end
        while denominator(s) % 3 == 0
            pow = cbrt(pow)
            s *= 3
        end
        if denominator(s) == 1
            pow = pow^numerator(s)
            s = 1
        end
    end
    pow^s
end


# nu_single - returns value of the nu function for a vector of delta 
# coefficients (delta_N), with truncation parameter (n_cut) at the 
# value (s).
function nu_single(delta_N::Vector{T},n_cut::Integer,s::Number) where {T<:AbstractFloat}
    mu = zeros(T,n_cut)
    for n in 1:n_cut
        divi = divis(n)
        mu[n] = sum([moebius_mu(Int(n/d)) * (d <= length(delta_N) ? delta_N[d] : 0) for d in divi])
    end
    sum([mu[n]*power_23adic(T(n),-s) for n in axes(mu,1)])
end


# nu - returns dictionary; the key (M,K,s) contains the value
# of nu_single at (M,K,s). The parameters M,K,s are taken within
# the ranges (range_M), (range_K), (range_s), respectively, and
# (delt) containes the coefficients. The computation is parallelized.
function nu(range_M,range_K,range_s,delt)
    I= Iterators.product(range_M,range_K,range_s) 
    nu_dic = Dict(i => BigFloat() for i in I)   
    @floop for (M,K,s) in I
        nu_dic[M,K,s] = 2*K*nu_single(delt[M],K,s)
    end
    nu_dic
end


# scf - returns the simple continued fraction of a floating point (fl)
# up to depth (n_depth), the output is a vector where the first number 
# is the integer part of (fl), and the next numbers correspond to the
# elements of the continued fraction.
function scf(fl::T,n_depth::Integer) where {T<:AbstractFloat}
    q = zeros(Integer,n_depth)
    u = zeros(T,n_depth)
    x = fl

    for k in 1:n_depth
        q[k] = Int(floor(x))
        u[k] = x - q[k]
        u[k] !=0 ? x = 1/u[k] : break
    end
    q
end


# partition_dict - given a dictionary (dic), outputs a list of subdictionaries 
# of (dic) of size in bytes of at most (chunk_size). 
function partition_dict(dic::Dict, typical_size::Int, chunk_size::Int)
    partition = []
    tmp_dic = Dict()
    n_entries_per_chunk = chunk_size ÷ typical_size
    n_chunks = Int(ceil(length(dic)/n_entries_per_chunk))
    count_chunks = 0
    count = 0

    for entry in dic
        push!(tmp_dic,entry)
        count += 1
        if count == n_entries_per_chunk
            push!(partition,tmp_dic)
            tmp_dic = Dict()
            count = 0
            count_chunks += 1
        end
    end

    if n_chunks - count_chunks == 1
        push!(partition,tmp_dic)
    end
    
    partition
end


# write_coefs - write the dictionaries of values of nu contained in the
# list (partition) into chunk files. 
function write_nu(partition,range_M,range_K,range_s,dps)
    function str_range(rang)
        "$(rang[1])-$(step(rang))-$(rang[end])"
    end

    partition = map(JSON.json,partition)
    str_part = length(partition) > 1 ? ["_part_$(i)" for i in 1:length(partition)] : [""]
    for i in 1:length(partition)
        open("Data/p$(dps)/nu_M$(str_range(range_M))_K$(str_range(range_K))_s$(str_range(range_s))_p$(dps)$(str_part[i]).json", "w") do f
            write(f, partition[i])
        end
    end
end


function main()
    # Parameters
    range_M = 1:1500
    range_K = 1:3000
    range_s = 1:1
    n_dps = 10_000
    chunk_size = Int(5e9)

    st = time()
    setprecision(BigFloat,n_dps;base=10)

    # Open CoefDelta files
    delta = open_coefdelta(range_M,n_dps)
    println("Imported delta coefficients: $(time()-st) s.")

    # Computation of nu
    nu_dic = nu(range_M,range_K,range_s,delta)
    println("Computed nu: $(time()-st) s.")

    # Write JSON files
    nu_str = Dict(key => string(value) for (key,value) in nu_dic)
    typical_size = sizeof(nu_str[range_M[end],range_K[end],range_s[end]])
    partition = partition_dict(nu_str,typical_size,chunk_size)
    write_nu(partition,range_M,range_K,range_s,n_dps)
    println("Total time: $(time()-st) s.")
end


main()
