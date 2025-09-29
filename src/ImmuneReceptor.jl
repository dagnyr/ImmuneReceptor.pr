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

# to do - make um a range from min to max motif size + add loop to find motifs in this range

function get_motif_new(st_::AbstractVector{<:AbstractString}, min::Int, max::Int)

    mo_ = Dict(i => Dict{String,Int}() for i in min:max) # make dictionary

    for s1 in st_

        lastindex(s1) < 7 && continue # only keep cdr3 7 or longer

        s2 = s1[4:(end-3)] # remove first and last 3 aa
        for um in min:max

            um > lastindex(s2) && continue # make sure cdr3 is logner than the motif size

            i1 = 0
            i2 = i1 + um - 1

            while i2 < lastindex(s2)

                m = s2[(i1+=1):(i2+=1)] # get the motif
                mo_[um][m] = get(mo_[um], m, 0) + 1 # count motif occurances

            end

        end

    end

    return mo_

end

#function count_motif(listofstrings)
# TO DO:
# - for every motif -> count occurances of each motif through each mo_ motif length pool.
# also note (make dictionary) length of the CDR3 AND index of each occurance of each motif in the string.
# later, use this index + length to helo draw local edges (draw edge if length is same, motif is roughly in same position)
#end



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

end
