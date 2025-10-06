module ImmuneReceptor

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

# ----------------------------------------------------------------------------------------------- #

using ProgressMeter: @showprogress

using Nucleus

# =============================================================================================== #
# Reading
# =============================================================================================== #

function is_t_gene(::Any)

    false

end

function is_t_gene(st::AbstractString)

    startswith(st, 'T')

end

function is_cdr3(st)

    st[1] == 'C' && st[end] == 'F'

end

function load_cdr3s(csvpath,)
    df = CSV.read(csvpath),DataFrame)
    filtered_df = df[
        (df.chain .∈ ["TRG","TRD","TRA","TRB"]) .&
        .!(df.tcr_gene .∈ ["None", "none","NA", "na", " "]) .&
        .!(df.cdr3 .∈ ["None", "none","NA", "na", " "]) .&
        startswith.(df.cdr3, "C") .&
        endswith.(df.cdr3, "F"),
    :]
    filtered_df_2 = combine(groupby(filtered_df, :cdr3), nrow => :duplicate_count)
    return filtered_df_2
end

# =============================================================================================== #
# Writing
# =============================================================================================== #

function make_trace(an_)

    um = 1000000

    Dict(
        "name" => "Naive",
        "type" => "histogram",
        "histnorm" => "probability",
        "x" => lastindex(an_) <= um ? an_ : rand(an_, um),
    )

end

function writ(ht, st, a1_, a2_)

    Nucleus.Plotly.writ(
        ht,
        (make_trace(a1_), make_trace(a2_)),
        Dict(
            "yaxis" => Dict("title" => Dict("text" => "Probability")),
            "xaxis" => Dict("title" => Dict("text" => st)),
        ),
    )

end

# =============================================================================================== #
# VJ
# =============================================================================================== #

function make_dictionary(an_)

    Dict(an_[nd] => nd for nd in eachindex(an_))

end

function make_vj(s1_, s2_)

    u1_ = unique(s1_)

    u2_ = unique(s2_)

    U = zeros(Int, lastindex(u1_), lastindex(u2_))

    d1 = make_dictionary(u1_)

    d2 = make_dictionary(u2_)

    for nd in eachindex(s1_)

        U[d1[s1_[nd]], d2[s2_[nd]]] += 1

    end

    u1_, u2_, U

end

# =============================================================================================== #
# CDR3
# =============================================================================================== #

function make_hamming_distance(s1, s2)

    sum(s1[nd] != s2[nd] for nd in eachindex(s1))

end

function make_distance(st_)

    u1 = lastindex(st_)

    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)

    po_ = Vector{Int}(undef, u2)

    i1 = 0

    @showprogress for i2 in 1:u1, i3 in (i2+1):u1

        s1 = st_[i2]

        s2 = st_[i3]

        if lastindex(s1) != lastindex(s2)

            continue

        end

        in__[i1+=1] = i2, i3

        po_[i1] = make_hamming_distance(s1, s2)

    end

    resize!(in__, i1), resize!(po_, i1)

end

# =============================================================================================== #
# Motif
# =============================================================================================== #

# count whether motif is present at least once per cdr3 for all cdr3s
function get_motif_counts(
    motif::AbstractVector{<:AbstractString},
    cdrs::AbstractVector{<:AbstractString},
)

    motif_counts = Dict{String,Int}()

    for m in motif

        motif_counts[m] = count(cdr -> occursin(m, cdr), cdrs) # count occurances per cdr3 for each motif

    end

    return motif_counts

end

# makes list of motifs in cdr3s, counts them, filters based on a cutoff, then provides set count
function get_motifs(st_::AbstractVector{<:AbstractString}, min::Int, max::Int)

    mo_ = Dict{String,Int}() # make dictionary

    for s1 in st_

        lastindex(s1) < 7 && continue # only keep cdr3 7 or longer

        s2 = s1[4:(end-3)] # remove first and last 3 aa
        for um in min:max

            um > lastindex(s2) && continue # make sure cdr3 is logner than the motif size

            i1 = 0
            i2 = i1 + um - 1

            while i2 < lastindex(s2)

                m = s2[(i1+=1):(i2+=1)] # get the motif
                mo_[m] = get!(mo_, m, 0) + 1 # count total motif occurances

            end

        end

    end

    mo_ = Dict(m => num for (m, num) in mo_ if num >= 3) # cutoff of 3
    mo2_ = get_motif_counts(collect(keys(mo_)), st_)
    return mo2_

end

# input will be motif_counts dictionary/table
function find_significant_motifs(motifs, cdrs1, cdrs2)

    Random.seed!(1)

    motifs_list = collect(keys(motifs))
    counts_orig = collect(values(motifs))

    counts_sim = Array{Float64}(undef, 1000, length(motifs_list))

    significant_motifs = Dict{String,Float64}()

    for i in 1:1000

        random_cdrs = sample(cdrs2, length(cdrs1); replace=true, ordered=false)
        random_counts = get_motif_counts(motifs_list, random_cdrs)
        counts_sim[i, :] = collect(values(random_counts))

    end

    for (index, m) in enumerate(motifs_list)
        ove = counts_orig[index] / mean(counts_sim[:, index])

        if counts_orig[index] < 2
            continue
        elseif (counts_orig[index] == 2 && ove >= 1000) ||
               (counts_orig[index] == 3 && ove >= 100) ||
               (counts_orig[index] >= 4 && ove >= 10)

            wins = count(x -> x >= counts_orig[index], counts_sim[:, index])
            p_val = (wins + 1) / (1000 + 1)
            if p_val >= 0.05
                significant_motifs[m] = get!(significant_motifs, m, 0.0) + p_val
            end

        end

    end

    return significant_motifs

end


#_________#

# takes list of cdrs and checks for a motif, making pairs
function make_motif_pairs(st_::AbstractVector{<:AbstractString}, motif)
    u1 = lastindex(st_)
    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)
    po_  = Vector{Int}(undef, u2)

    i1 = 0

    for i2 in 1:u1, i3 in (i2+1):u1
        s1, s2 = st_[i2], st_[i3]

        has1 = occursin(motif, s1)
        has2 = occursin(motif, s2)

        if has1 && has2
            i1 += 1
            in__[i1] = (i2, i3)
            po_[i1] = 1
        end
    end

    resize!(in__, i1)
    resize!(po_, i1)

    return in__, po_
end


# TO DO:
# 1) using counts table from above, run simulations on reference dataset + count occurances of all motifs (say sim depth = 1000 for now) ✅
# 1.5) find significant motifs base on these simulations ✅
# 2) filter the motif table from above and then filter the positions table based on the filtered count table (NOT DOING POSITION TRACKING FOR NOW) ❎
# 3) using the filtered significant motifs positions table, then draw local edges for cdr3s with:
#   a) same length
#   b) motif appearance has max overlap difference of x amino acids

#_________#

function get_motif(s1::AbstractString, um)

    s2 = s1[4:(end-3)]

    # TODO: Use Set
    st_ = String[]

    i1 = 0

    i2 = i1 + um - 1

    while i2 < lastindex(s2)

        push!(st_, s2[(i1+=1):(i2+=1)])

    end

    st_

end

function get_motif(st__, u1)

    di = Dict{String,Int}()

    for nd in eachindex(st__), st in st__[nd]

        if !haskey(di, st)

            di[st] = 0

        end

        di[st] += 1

    end

    [st for (st, u2) in di if u1 <= u2]

end


# =============================================================================================== #
# Graphing
# =============================================================================================== #

function add_vertices!(g, cdrs)
    for vertex in 1:nv(g)
        set_prop!(g, vertex, :label, cdrs.cdr3[vertex])
        set_prop!(g, vertex, :v_gene, cdrs.v_gene[vertex])
        set_prop!(g, vertex, :j_gene, cdrs.j_gene[vertex])
        if !(cdrs.d_gene[vertex] in ["None", "none", "NA", "na", " "])
            set_prop!(g, vertex, :d_gene, cdrs.d_gene[vertex])
        end
    end
    return g
end

# make a global edge w/ annotation of hamming disrance
function add_global_edge!(g, u, v, distance::Int64)
    add_edge!(g, u, v)
    set_prop!(g, (u, v), :type, "global")
    set_prop!(g, (u, v), :distance, distance)
end

# add a local edge w/ annotation of motif and the pval (from get_significant_motifs)
function add_local_edge!(g, u, v, motif::String, pval::Int64)
    add_edge!(g, u, v)
    set_prop!(g, (u, v), :type, "local")
    set_prop!(g, (u, v), :motif, motif)
    set_prop!(g, (u, v), :pval, pval)
end

# go through cdrs
function make_edges(cdrs, motifs)
    g = MetaGraph(Graph(length(cdrs)))
    g = add_vertices(g, cdrs)

    # start w/ global distances
    pairs, dists = make_distance(cdrs.cdr3)
    for (index, (cdr1, cdr2)) in enumerate(pairs)
        d = dists[index]
        if dists[index] <= 1
            add_global_edge!(g, cdr1, cdr2, dists[index])
        end
    end

    # TODO: next do local edges
    for (m, pval) in motifs
        pairs2, hasmotif = make_motif_pairs(cdrs.cdr3, m)

        for (index, (cdr1, cdr2)) in enumerate(pairs2)
            if hasmotif[index] == 1
                add_local_edge!(g, cdr1, cdr2, m, pval)
            end

        end

    end

    return g

end

# only make global edges
function make_global_edges(cdrs, motifs)
    g = MetaGraph(Graph(length(cdrs)))
    g = add_vertices(g, cdrs)

    # start w/ global distances
    pairs, dists = make_distance(cdrs.cdr3)
    for (index, (cdr1, cdr2)) in enumerate(pairs)
        d = dists[index]
        if dists[index] <= 1
            add_global_edge!(g, cdr1, cdr2, dists[index])
        end
    end

    return g

end

# only make local edges
function make_local_edges(cdrs, motifs)
    g = MetaGraph(Graph(length(cdrs)))
    g = add_vertices(g, cdrs)

    # local edge making
    for (m, pval) in motifs
        pairs2, hasmotif = make_motif_pairs(cdrs.cdr3, m)

        for (index, (cdr1, cdr2)) in enumerate(pairs2)
            if hasmotif[index] == 1
                add_local_edge!(g, cdr1, cdr2, m, pval)
            end

        end

    end

    return g

end

# to-do return list of communities and generate output file w/ all values
#components = connected_components(g)


# =============================================================================================== #
# Clusters / Significance
# =============================================================================================== #

# length score
# for each group, resample 1000 groups of group size n and measure cdr3 length distribution.

using StatsBase: mode

function score_lengths(g, cdrs2)
    clusters = connected_components(g)

    # TO DO GET PROPERTY LAB FROM EACH VERTEX get_prop(g, vertex_id, :label)

    modes = Vector{Int}()
    props = Vector{Float64}()

    # get most common length for each cluster
    for cluster in clusters
        labels = get_prop(g, v, :label) for v in cluster # fixed this to get labels
        cluster_lengths = [length(s) for s in labels]
        m = mode(cluster_lengths)
        push!(modes, m)
        push!(props, count(==(m), cluster_lengths) / length(cluster))

    end

    # randomly pick n sequences from the null distribution and compare to
    cluster_sizes = [length(cluster) for cluster in clusters]
    p_vals = Float64[]

    for (index, size) in enumerate(cluster_sizes)

        sim_props = Float64[]

        for i in 1:1000
            random_cdrs = sample(cdrs2, size; replace=true, ordered=false)
            cluster_lengths = [length(s) for s in random_cdrs]
            push!(sim_props, (count(==(modes[index]), cluster_lengths)) / size)
        end

        # count number of times sim beats data?
        p_val = ((count(>=(props[index]), sim_props)) + 1) / 1001
        push!(p_vals, p_val)

    end

    return p_vals
    # maybe add the p-vals as a attribute in the graph itself?

end

# v gene score
# hyper geomtric or parameteric when over 200 sequences in one group

# TO DO: lowkey change entire funtion, took wrong approach
function score_vgene(g)
    clusters = connected_components(g)

    # make vector w/ dictionary of vgenes
    counts = Vector{Dict{String,Int}}(undef, length(clusters))
    for (index, cluster) in enumerate(clusters)
        di = Dict{String,Int}()
        for vertex in cluster
            vgene = get_prop(g, vertex, :v_gene)
            di[vgene] = get!(di, vgene, 0) + 1
        end
        counts[index] = di
    end

    props = Vector{Dict{String,Float64}}(undef, length(clusters))
    for (index, di) in enumerate(counts)
        total_vgenes = sum(values(di))

        di2 = Dict{String,Int}()

        for (vgene, count) in di

            prop = count / total_vgenes
            di2[vgene] = get!(di2, vgene, 0) + prop

        end

    end

    # after getting dictionary of vgenes get proportion
    # after counting the vgenes compare to the count in the null

    # redo count for the null distribution


    # make
    # dictionary with every cluster -> list of v-genes for every cdr3
    # for every cluster, count total occurances and determine if over 200
    # if under 200, calculate p-val for every v-gene and report lowest value + gene
    # if over 200, run simulation 1000 times and count, make pvals, etc. (just like length)


    if any(c -> c > 200, [length(cluster) for cluster in clusters])

    else

    end



    # first check if any cluster has > 200 members and if so, use empirical calculation

end

# hla score
#

function score_hla()

    # use sample or donor id from original csv
    # load separare hla typing table
    # figure out best way to determine if a type is enriched - maybe population level HLA distirbution (aka by chance)
    # or maybe make a random cdr cluster w/ random draws and see HLA typing (might be problematic if few donors in sample)

end
