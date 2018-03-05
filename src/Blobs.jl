using BackgroundSegmenter
using DataStructures

const AGE_FOR_TRANSIENCE = 2

struct Blob
    centroid::Vector{Float64}
    ratio::Float64
    area::Float64
    source::Int
    θ::Float64
    shape::Array{Float64, 2}
    points::Array{Int, 2}
end

struct Object
    centroid::Vector{Float64}
    ratio::Float64
    area::Float64
    source::Int
    θ::Float64
    points::Array{Int, 2}
    label::Int
    age::Int
end
Object(b::Blob, label=-1, age=0) = Object(b.centroid, b.ratio, b.area, b.source, b.θ, b.points, label, age)

update(object::Object, blob::Blob) = Object(blob, object.label, object.age+1)
is_transient(object::Object) = object.age < AGE_FOR_TRANSIENCE

Blob(centroid, ratio, area, source, θ) = Blob(centroid, ratio, area, source, θ,
                                              Array{Float64}(0, 0),
                                              Array{Int}(0, 0))

function Blob(source::Int, points)
    centroid = vec(mean(points, 1)) 
    area = size(points)[1]
    if area > 1
        covariance = cov(points, 1)
    else
        covariance = [1.0 0.0; 0.0 1.0]
    end
    F = eigfact(covariance)
    ratio = (max(F[:values]...) + 0.5) / (min(F[:values]...) + 0.5)
    if F[:values][1] < F[:values][2]
        v = F[:vectors][:, 2]
        w = F[:vectors][:, 1]
    else
        v = F[:vectors][:, 1]
        w = F[:vectors][:, 2]
    end
    basis = vcat(v', w')
    if v[2] == 0.0
        θ = π/2
    else
        θ = atan(v[1] / v[2])
    end
    shape = calc_shape(points, basis, centroid)
    return Blob(centroid, ratio, area, source, θ, shape, points)
end

function combine(b1::Blob, b2::Blob)
    if b1.area > b2.area
        source = b1.source
    else
        source = b2.source
    end
    points = vcat(b1.points, b2.points)
    return Blob(source, points)
end

function combine(objects)
    transient = true
    source_votes = DefaultDict(Int, Int, 0)
    for object in objects
        source_votes[object.source] += object.area
        if !is_transient(object)
            transient = false
        end
    end
    max_vote = 0
    max_source = 0
    for (source, vote) in source_votes
        if vote > max_vote
            max_vote = vote
            max_source = source
        end
    end
    points = [object.points for object in objects]
    points = vcat(points...)
    blob = Blob(max_source, points)
    return Object(blob, transient)
end

function distance(x::Object, y::Blob)
    sqrt(sum((x.centroid - y.centroid).^2))
end

function cost(x::Object, y::Blob)
    c = 0.0
    c += abs(x.ratio - y.ratio)/x.ratio
    c += abs(x.area - y.area) / x.area
    if x.source != y.source
        c += 1
    end
    c += 1 - (cos(x.θ - y.θ))^2
    return c
end

function composite_cost_comparison(objects, blob::Blob)
    min_cost = Inf
    best_idx = -1
    for (n, object) in enumerate(objects)
        object_cost = cost(object, blob)
        if object_cost < min_cost
            min_cost = object_cost
            best_idx = n
        end
    end
    composite = combine(objects)
    best_match = objects[best_idx]
    c = 0.0
    c += abs(composite.ratio - blob.ratio)/best_match.ratio
    c += abs(composite.area - blob.area) / best_match.area
    if composite.source != blob.source
        c += 1
    end
    c += 1 - (cos(composite.θ - blob.θ))^2
    if c < min_cost
        return -1
    else
        return best_idx
    end
end

function calc_shape(points, basis, centroid)
    return Array{Float64}(0, 0)
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
