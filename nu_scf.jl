using BenchmarkTools, JSON, Nemo #, Serialization

n_dps = 10_000
setprecision(BigFloat,n_dps;base=10)

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

function nu(delta_N::Vector{T},n_cut::Integer,s::Number) where {T<:AbstractFloat}
    mu = zeros(T,n_cut)
    @Threads.threads for n in 1:n_cut
        divi = divis(n)
        mu[n] = sum([moebius_mu(Int(n/d)) * (d <= length(delta_N) ? delta_N[d] : 0) for d in divi])
    end
    sum([mu[n]*T(n)^(-s) for n in axes(mu,1)])
end

function scf(fl::T,n_depth::Integer) where T<:AbstractFloat
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

function main()
    max_M = 1_300
    n_cut_incr = 100
    s_max = 2
    #n_depth = 20

    n_cut_max = 2*max_M ÷ n_cut_incr * n_cut_incr
    rang_cut = n_cut_incr:n_cut_incr:n_cut_max
    rang_s = range(-s_max,s_max)
    
    ub = [0,800,1100,1300,1500,1700,1800,2000,2150,2300,2400,2500,2600]
    frag = [(ub[i]+1,ub[i+1]) for i in 1:length(ub)-1]

    delta = Dict((Mlo,Mup) => Dict{Integer,Vector{BigFloat}}() for (Mlo,Mup) in frag)
    nu_pre = Dict((Mlo,Mup) => Dict{Tuple{Integer,Integer,Integer},BigFloat}() for (Mlo,Mup) in frag)
    nu_post = Dict{Tuple{Integer,Integer,Integer},BigFloat}()   
    #nu_scf = Dict((Mlo,Mup) => Dict{Integer,Vector{BigInt}}() for (Mlo,Mup) in frag)

    st = time()
    for (Mlo,Mup) in frag
        # READ JSON FILES HERE! (Old variant: delta[(Mlo,Mup)] = open(deserialize, "C:\\Users\\user\\Documents\\c.delacruz\\JuliaYuri\\M$(max_M)_p$(n_dps)\\coef_delta_Mlo$(Mlo)_Mup$(Mup)_p$(n_dps).txt")#"/home/c.delacruz/JuliaYuri/M$(max_M)_p$(n_dps)/coef_delta_Mlo$(Mlo)_Mup$(Mup)_p$(n_dps).txt"))
        file_content = read("/home/c.delacruz/JuliaYuri/M$(max_M)_p$(n_dps)/coef_delta_Mlo$(Mlo)_Mup$(Mup)_p$(n_dps).json",String)
        data = JSON.parse(file_content)
        delta[(Mlo,Mup)] = Dict(parse(Int16,index) => map(v -> parse(BigFloat,v),value) for (index,value) in data)
    end
    println("Completed importing delta coefficients after $(time()-st) s.")

    st2 = time()
    for (Mlo,Mup) in frag
        for (N,delta_vec) in delta[(Mlo,Mup)]
            for K in rang_cut
                for s in rang_s
                    nu_pre[(Mlo,Mup)][(N,K,s)] = 2*K*nu(delta_vec,K,s)
                end
            end
            #nu_scf[(Mlo,Mup)][index] = scf(2*n_cut*nu(value,n_cut,s),n_depth)
            println("N=$(N) took $(time()-st2) s.")
        end
        #open("C:\\Users\\user\\Documents\\c.delacruz\\JuliaYuri\\M$(max_M)_p$(n_dps)\\nu_scf_Mlo$(Mlo)_Mup$(Mup)_cut$(n_cut)_dep$(n_depth)_p$(n_dps).txt", "w") do f #/home/c.delacruz/JuliaYuri/M$(max_M)_p$(n_dps)/nu_scf_Mlo$(Mlo)_Mup$(Mup)_cut$(n_cut)_dep$(n_depth)_p$(n_dps).txt", "w") do f
            #serialize(f,nu_scf[(Mlo,Mup)])
        #end
    end
    println("Completed computation of nu after $(time()-st) s.")

    nu_post = merge(values(nu_pre)...)

    json_nu = JSON.json(nu_post)
    open("nu.json", "w") do f
        write(f, json_nu)
    end
end

main()

# ######### CONVERSION FROM SERIALIZED FILES TO JSON #########
# for (Mlo,Mup) in frag
#     delta = open(deserialize, "C:\\Users\\user\\Documents\\c.delacruz\\JuliaYuri\\M$(max_M)_p$(n_dps)\\coef_delta_Mlo$(Mlo)_Mup$(Mup)_p$(n_dps).txt")#"/home/c.delacruz/JuliaYuri/M$(max_M)_p$(n_dps)/coef_delta_Mlo$(Mlo)_Mup$(Mup)_p$(n_dps).txt")
#     delta_str = Dict(index => map(string,value) for (index,value) in delta)
#     delta_json = json(delta_str)
#     open("coef_delta_Mlo$(Mlo)_Mup$(Mup)_p$(n_dps).json","w") do f
#         write(f,delta_json)
#     end
# end
# ###########################################################