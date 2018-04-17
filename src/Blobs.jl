using BackgroundSegmenter
using DataStructures

const AGE_FOR_TRANSIENCE = 2

struct Blob
    x::State
    y::State
    ratio::Float64
    area::Float64
    source::Int
    θ::Float64
    shape1::Array{Float64, 2}
    shape2::Array{Float64, 2}
    points::Array{Int, 2}
end

struct Object
    x::State
    y::State
    ratio::Float64
    area::Float64
    source::Int
    θ::Float64
    shape1::Array{Float64, 2}
    shape2::Array{Float64, 2}
    points::Array{Int, 2}
    label::Int
    age::Int
end
Object(b::Blob, label=-1, age=0) = Object(predict(b.x), predict(b.y), b.ratio, b.area, b.source, b.θ, b.shape1, b.shape2, b.points, label, age)

function update(object::Object, b::Blob, λ=0.9)
    x_velocity = smoothed_velocity_estimate(object.x, b.x, λ)
    y_velocity = smoothed_velocity_estimate(object.y, b.y, λ)
    x_observed = State(b.x.p, x_velocity)
    y_observed = State(b.y.p, y_velocity)
    Object(update(object.x, x_observed), update(object.y, y_observed),
           smooth(b.ratio, object.ratio, λ), smooth(b.area, object.area, λ),
           b.source, smooth(b.θ, object.θ, λ),
           b.shape1, b.shape2, b.points, object.label, object.age+1)
end

smooth(x, y, λ) = λ*x + (1 - λ)*y

is_transient(object::Object) = object.age < AGE_FOR_TRANSIENCE

Blob(centroid, ratio, area, source, θ) = Blob(centroid, ratio, area, source, θ,
                                              Array{Float64}(0, 0),
                                              Array{Float64}(0, 0),
                                              Array{Int}(0, 0))

function Blob(source::Int, points)
    centroid = vec(mean(points, 1)) 
    x = State(centroid[1], 0.0)
    y = State(centroid[2], 0.0)
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
    w .= orientation(v, w) .* w
    basis1 = hcat(v, w)
    basis2 = hcat(-v, -w)
    if v[2] == 0.0
        θ = π/2
    else
        θ = atan(v[1] / v[2])
    end
    shape1 = calc_shape(points, basis1, reshape(centroid, 1, 2))
    shape2 = calc_shape(points, basis2, reshape(centroid, 1, 2))
    return Blob(x, y, ratio, area, source, θ, shape1, shape2, points)
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
    source_votes = DefaultDict(Int, Float64, 0.0)
    for object in objects
        source_votes[object.source] += object.area
        if !is_transient(object)
            transient = false
        end
    end
    max_vote = 0.0
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

function distance(obj::Object, b::Blob)
    sqrt((obj.x.p - b.x.p)^2 + (obj.y.p - b.y.p)^2)
end

function cost(obj::Object, b::Blob)
    scale = (1.0 - pdf(obj.x, b.x.p) * pdf(obj.y, b.y.p)) * distance(obj, b)
    c = 0.0
    c += abs(obj.ratio - b.ratio) / obj.ratio
    c += abs(obj.area - b.area) / obj.area
    if obj.source != b.source
        c += 1
    end
    c += 1 - (cos(obj.θ - b.θ))
    if obj.shape1[1, :]' * b.shape1[1, :] > 0
        c += (1.0 - intersection_over_union(obj.shape1, b.shape1)) / 2
    else
        c += (1.0 - intersection_over_union(obj.shape1, b.shape2)) / 2
    end
    return c * scale
end

function orientation(x::Vector{T}, y::Vector{T}) where {T}
    x1 = zeros(3)
    y1 = zeros(3)
    x1[1:2] = x
    y1[1:2] = y
    z = cross(x1, y1)
    if z[3] > zero(T)
        return one(T)
    else
        return -one(T)
    end
end

function intersection_over_union(x, y)
    if length(x) < length(y)
        return intersection_over_union(y, x)
    end
    x = x[2:end, :]
    y = y[2:end, :]
    (a, _) = size(x)
    (b, _) = size(y)
    c1 = (a + 1.0) / 2.0
    c2 = (b + 1.0) / 2.0
    offset = convert(Int, round(c1 - c2))
    intersection = 0.0
    union = 0.0
    for j in 1:offset
        union += 2x[j, 2]
    end
    for j in (b + offset + 1):a
        union += 2x[j, 2]
    end
    for j in 1:b
        i, u = intersection_and_union(x[j+offset, :], y[j, :])
        intersection += i
        union += u
    end
    if union > 0.0
        return intersection / union
    else
        return 0.5
    end
end

function intersection_and_union(x, y)
    (c1, w1) = x
    (c2, w2) = y
    if c2 < c1
        return intersection_and_union(y, x)
    end
    union = 2.0 * (w1 + w2)
    if (c1 + w1) <= (c2 - w2)
        intersection = 0.0
    elseif (c1 - w1) >= (c2 - w2)
        intersection = 2w1
    elseif (c1 + w1) >= (c2 + w2)
        intersection = 2w2
    else
        intersection = (c1 + w1) - (c2 - w2)
    end
    union -= intersection
    return (intersection, union)
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
    c = cost(composite, blob)
    if c < min_cost
        return -1
    else
        return best_idx
    end
end

function calc_shape(points, basis, centroid)
    deskewed_points = (points .- centroid) * basis
    min_idx = Inf
    max_idx = -Inf
    (n, m) = size(deskewed_points)
    ranges = Dict{Int, Vector{Float64}}()
    for i in 1:n
        idx = convert(Int, round(deskewed_points[i, 1]))
        if idx < min_idx
            min_idx = idx
        end
        if idx > max_idx
            max_idx = idx
        end
        if !haskey(ranges, idx)
            ranges[idx] = [deskewed_points[i, 2]]
        else
            push!(ranges[idx], deskewed_points[i, 2])
        end
    end
    shape = Array{Float64}(max_idx - min_idx + 2, 2)
    shape[1, :] = basis[:, 1]'
    for (k, idx) in enumerate(min_idx:max_idx)
        if haskey(ranges, idx)
            b = max(ranges[idx]...)
            a = min(ranges[idx]...)
            shape[k+1, 1] = (b + a) / 2
            shape[k+1, 2] = (b - a) / 2 
        else
            shape[k+1, 1] = 0.0
            shape[k+1, 2] = 0.0
        end
    end
    return shape
end

function form_blobs(w::Array{T, 3}) where {T}
    (n, m, t) = size(w)
    blobs = Vector{Vector{Blob}}(t)
    for j in 1:t
        blobs[j] = form_blobs(view(w, :, :, j))
    end
    return blobs
end

function form_blobs(w::AbstractArray{T, 2}) where {T}
    blobs = Vector{Blob}(0)
    C = components(w)
    for (label, points) in C
        push!(blobs, Blob(label, points))
    end
    return blobs
end
