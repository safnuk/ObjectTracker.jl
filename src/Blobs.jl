
struct Blob
    centroid::Tuple{Float64, Float64}
    ratio::Float64
    area::Float64
    source::Int
    θ::Float64
end

function form_blobs(w::Array{T, 3}) where {T}
    (t, n, m) = size(w)
    blobs = Vector{Vector{Blob}}(t)
    for j in 1:t
        blobs[j] = form_blobs(w[j, :, :])
    end
    return blobs
end

function form_blobs(w::Array{T, 2}) where {T}
    blobs = Vector{Blob}(0)
    push!(blobs, Blob((0.0, 0.0), 1.0, 3, 1, 0.0))
    return blobs
end
