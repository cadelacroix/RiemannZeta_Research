using FLoops, Nemo, JSON, BenchmarkTools

function interval(filename::String)
    str_match = match(r"CoefDelta_M(\d+)-(\d+)_p(\d+)",filename)
    if !isnothing(str_match)
        map(x->parse(Int64,x), str_match.captures)
    end
end


function interval_list(dps::Int)
    rawlist = Vector{Vector{Int64}}(filter(!isnothing,map(interval,readdir("../Data/p$(dps)/"))))
    map(x -> (x[1],x[2]),filter(x -> x[3] == dps, rawlist))
end


function get_chain(intervals::Vector{Tuple{Int,Int}}, m::Int, M::Int)
    sorted_intervals = sort(intervals)
    
    # Find the interval containing m
    i_m = findfirst(x -> (x[1] <= m <= x[2]), sorted_intervals)
    i_M = findlast(x -> (x[1] <= M <= x[2]), sorted_intervals)
    
    # If M wasn't found, return empty array
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


function open_coefdelta(range_M::AbstractRange{Int},dps::Int)
    intervals = interval_list(dps)
    chain = get_chain(intervals,range_M[1],range_M[end])

    data = Dict((Mlo,Mup) => Dict{String,Vector{String}}() for (Mlo,Mup) in chain)
    delta_pre = Dict((Mlo,Mup) => Dict{Int,Vector{BigFloat}}() for (Mlo,Mup) in chain)
    @floop for (Mlo,Mup) in chain
        data[Mlo,Mup] = JSON.parse(read("../Data/p$(dps)/CoefDelta_M$(Mlo)-$(Mup)_p$(dps).json",String))
        delta_pre[Mlo,Mup] = Dict(parse(Int,key) => map(v -> parse(BigFloat,v),value) for (key,value) in data[Mlo,Mup])
    end
    merge(values(delta_pre)...)
end


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


function nu_single(delta_N::Vector{T},n_cut::Integer,s::Number) where {T<:AbstractFloat}
    mu = zeros(T,n_cut)
    for n in 1:n_cut
        divi = divis(n)
        mu[n] = sum([moebius_mu(Int(n/d)) * (d <= length(delta_N) ? delta_N[d] : 0) for d in divi])
    end
    sum([mu[n]*power_23adic(T(n),-s) for n in axes(mu,1)])
end


function nu(range_M,range_K,range_s,delt)
    I= Iterators.product(range_M,range_K,range_s) 
    nu_dic = Dict(i => BigFloat() for i in I)   
    @floop for (M,K,s) in I
        nu_dic[M,K,s] = 2*K*nu_single(delt[M],K,s)
    end
    nu_dic
end


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


function str_range(rang)
    "$(rang[1])-$(step(rang))-$(rang[end])"
end


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
        end
    end

    if n_chunks - count_chunks == 1
        push!(partition,tmp_dic)
    end
    
    partition
end


function write_nu(partition,range_M,range_K,range_s,dps)
    partition = map(JSON.json,partition)
    str_part = length(partition) > 1 ? ["_part_$(i)" for _ in partition] : [""]
    @floop for i in axes(partition,1)
        open("../Data/p$(dps)/nuM$(str_range(range_M))_K$(str_range(range_K))_s$(str_range(range_s))_p$(dps)$(str_part[i]).json", "w") do f
            write(f, partition[i])
        end
    end
end


function main()
    range_M = 1:10
    range_K = 1:100
    range_s = 1:1
    n_dps = 40_000
    chunk_size = Int(3e9)

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