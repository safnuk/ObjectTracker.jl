
struct Blob
    centroid::Vector{Float64}
    ratio::Float64
    area::Float64
    source::Int
    θ::Float64
end

function Blob(source::Int, points)
    centroid = vec(mean(points, 1)) 
    area = size(points)[1]
    covariance = cov(points, 1)
    F = eigfact(covariance)
    ratio = (max(F[:values]...) + 0.5) / (min(F[:values]...) + 0.5)
    if F[:values][1] < F[:values][2]
        v = F[:vectors][:, 2]
    else
        v = F[:vectors][:, 1]
    end
    if v[2] == 0.0
        θ = π/2
    else
        θ = atan(v[1] / v[2])
    end
    return Blob(centroid, ratio, area, source, θ)
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
    C = components(w)
    for (label, points) in C
        push!(blobs, Blob(label, points))
    end
    return blobs
end

function components(w)
    X = [5 1; 5 2; 5 3]
    return Dict(1 => X)
end
