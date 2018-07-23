
# Function with equal functionality, but allocating memory
restriction_solve(b::Vector, A::DictionaryOperator, BE::IndexExtensionOperator,
        GR::IndexRestrictionOperator; cutoff=FrameFun.default_cutoff(A), options...) =
    BE*LAPACK.gelsy!(matrix(GR*A*BE),GR*b,cutoff)[1]

function decomposition_solve(b::Vector, A::DictionaryOperator, cart_indices, classified_indices; cutoff=default_cutoff(A), verbose=false, info=false, options...)
    bins = unique(classified_indices)
    # assign sequence numbers to each of the bins.
    seq_nr = assign_sequence_nro(bins)

    primal = basis(src(A))
    omega = grid(dest(A))
    g = supergrid(omega)
    x1 = zeros(primal)
    b_ =  copy(b)
    t = similar(x1)

    for d in 1:1<<length(bins[1])
        # first solve all parts with a low sequence number
        bpart = bins[find(seq_nr.==d)]
        # Each of the parts with an equal sequence number can be solved independently
        verbose && println("$(length(bpart)) parts in $(d)th flow")
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            xx, yy = FrameFun._azselection_restriction_operators(primal, g, omega, mo)
            verbose && println("\t$(i)\t has size ($(size(yy,1)),$(size(xx,1)))")
            op = yy*A*xx'
            a = matrix(op)
            y = LAPACK.gelsy!(a, yy*b_, cutoff)[1]
            # x1 = x1 + xx'*y
            apply!(xx', t, y)
            x1 .+= t
        end
        # Remove the solved part after all parts with equal seq number are dealt with.
        if d!=1<<length(bins[1])
            # b_ = b-A*x1
            apply!(A, b_, x1)
            b_ .= b .- b_
        end
    end
    x1
end

function decomposition_info(b::Vector, A::DictionaryOperator, cart_indices, classified_indices; cutoff=default_cutoff(A), info=false, options...)
    bins = unique(classified_indices)
    seq_nr = assign_sequence_nro(bins)
    primal = basis(src(A))
    omega = grid(dest(A))
    g = supergrid(omega)
    r = zeros(Int, length(bins), 2)
    ii = 1
    for d in 1:1<<length(bins[1])
        bpart = bins[find(seq_nr.==d)]
        println("$(length(bpart)) parts in $(d)th flow")
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            xx, yy = FrameFun._azselection_restriction_operators(primal, g, omega, mo)
            op = yy*A*xx'
            println("\t$(i)\t has size ($(size(yy,1)),$(size(xx,1)))")
            r[ii,1] = size(yy,1)
            r[ii,2] = size(xx,1)
            ii += 1
        end
    end
    r
end

using Plots
function decomposition_plot(A::DictionaryOperator, cart_indices, classified_indices; max_plots=nothing, options...)
    bins = unique(classified_indices)
    seq_nr = assign_sequence_nro(bins)
    # clrs = [:blue,:red,:green,:black]
    clrs = Plots.colormap("blues",1+1<<length(bins[1]))[2:end]
    clrs = clrs[end:-1:1]
    seq_nr = assign_sequence_nro(bins)
    primal = basis(src(A))
    (max_plots==nothing) && (max_plots = 1<<length(bins[1]))
    plot()
    for d in 1:min(max_plots, 1<<length(bins[1]))
        bpart = bins[find(seq_nr.==d)]
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            m = falses(size(primal))
            m[mo] = true
            scatter!(m,c=clrs[d]; options...)
        end
    end
    scatter!()
end

# Function with equal functionality, but allocating memory
restriction_info(b::Vector, A::DictionaryOperator, BE::IndexExtensionOperator,
        GR::IndexRestrictionOperator; cutoff=FrameFun.default_cutoff(A), options...) =
    (println("Selection has size ($(size(GR,1)),$(size(BE,2)))"); [size(GR,1),size(BE,2)]')
