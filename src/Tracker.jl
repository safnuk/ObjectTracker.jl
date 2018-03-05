using DataStructures

const MIN_AREA = 10
const MAX_COST = 2.0
const RADIUS = 80

function match!(objects::Vector{Object}, blobs::Vector{Blob}, num_objects::Int;
                radius=RADIUS, min_area=MIN_AREA)
    possible_matches = limit_to_radius(objects, blobs, radius)
    matches = naive_match(objects, blobs, possible_matches)
    finished_merging = false
    remove_duplicate_matches!(matches, objects)
    while !finished_merging
        remove_duplicate_matches!(matches, objects)
        finished_merging = evaluate_object_merges(matches, objects, blobs)
    end
    # parition_merged_blobs!(merges, objects, blobs)
    merge_blob_fragments!(matches, objects, blobs, possible_matches)
    update_with_new_objects!(matches, objects, blobs, num_objects, min_area)
end

function limit_to_radius(objects, blobs, radius)
    possible_matches = zeros(Int, length(objects), length(blobs))
    for (n, object) in enumerate(objects)
        for (m, blob) in enumerate(blobs)
            if distance(object, blob) < radius
                possible_matches[n, m] = 1
            end
        end
    end
    return possible_matches
end

function naive_match(objects, blobs, possible_matches)
    matches = zeros(possible_matches)
    for (n, object) in enumerate(objects)
        min_cost = Inf
        min_blob = -1
        for (m, blob) in enumerate(blobs)
            if possible_matches[n, m] == 0
                continue
            end
            c = cost(object, blob)
            if c < min_cost
                min_cost = c
                min_blob = m
            end
        end
        if min_blob > 0
            matches[n, min_blob] = 1
        end
    end
    return matches
end

function remove_duplicate_matches!(matches, objects)
    (n, m) = size(matches)
    for j in 1:m
        number_matched = 0
        transients = Queue(Int)
        for i in 1:n
            if matches[i, j] != 0
                number_matched += 1
                if is_transient(objects[i])
                    enqueue!(transients, i)
                end
            end
        end
        if number_matched > length(transients)
            for t in transients
                matches[t, j] = 0
            end
        end
    end
    return matches
end

function evaluate_object_merges(matches, objects, blobs)
    finished_checking_merges = true
    (n, m) = size(matches)
    for j in 1:m
        if sum(view(matches, :, j)) > 1
            blob_matches = Vector{Object}()
            match_indices = Vector{Int}()
            for i in 1:n
                if matches[i, j] != 0
                    push!(blob_matches, objects[i])
                    push!(match_indices, i)
                end
            end
            best_match = composite_cost_comparison(blob_matches, blobs[j])
            if best_match == -1  # indicates best match is to composite object
                for i in 1:n
                    if matches[i, j] != 0
                        if best_match < 0
                            best_match = i
                        else
                            matches[i, j] = 0
                        end
                    end
                end
                objects[best_match] = combine(blob_matches)
            else
                best_match = match_indices[best_match]
                finished_checking_merges = false
                for i in 1:n
                    if best_match != i
                        matches[i, j] = 0
                    end
                end
            end
        end
    end
    return finished_checking_merges
end

function merge_blob_fragments!(matches, objects, blobs, possible_matches)
    (n, m) = size(matches)
    unsecure_matches = Vector{Int}()
    secure_matches = Vector{Int}()
    unallocated_blobs = Vector{Int}()
    for j in 1:m
        if sum(view(matches, :, j)) == 0
            push!(unallocated_blobs, j)
        else
            for i in 1:m
                if matches[i, j] != 0 
                    if is_transient(objects[i])
                        push!(unsecure_matches, j)
                    else
                        push!(secure_matches, j)
                    end
                    break
                end
            end
        end
    end
    attempt_merges!(secure_matches, unallocated_blobs, matches, objects, blobs, possible_matches)
    attempt_merges!(secure_matches, unsecure_matches, matches, objects, blobs, possible_matches)
    attempt_merges!(unsecure_matches, unallocated_blobs, matches, objects, blobs, possible_matches)
end

function attempt_merges!(sources, targets, matches, objects, blobs, possible_matches)
    (n, m) = size(matches)
    finished = false
    while !finished
        finished = true
        for idx in sources
            source_blob = blobs[idx]
            object_idx = 0
            for object_idx in 1:n
                if matches[object_idx, idx] != 0
                    break
                end
            end
            object = objects[object_idx]
            for target_idx in length(targets):-1:1
                blob_idx = targets[target_idx]
                if possible_matches[object_idx, blob_idx] == 0
                    continue
                end
                merged = merge_if_match_improves!(object, source_blob, blobs[blob_idx])
                if merged
                    matches[:, blob_idx] = 0
                    blobs[idx] = combine(source_blob, blobs[blob_idx])
                    deleteat!(targets, target_idx)
                    finished = false
                end
            end
        end
    end
end

function merge_if_match_improves!(object, source_blob, target_blob)
    original_cost = cost(object, source_blob)
    merged_blob = combine(source_blob, target_blob)
    merged_cost = cost(object, merged_blob)
    return merged_cost < original_cost
end

function update_with_new_objects!(matches, objects, blobs, num_objects, min_area)
    (n, m) = size(matches)
    for i in n:-1:1
        if sum(view(matches, i, :)) == 0
            deleteat!(objects, i)
        else
            j = 0
            for j in 1:m
                if matches[i, j] != 0
                    break
                end
            end
            objects[i] = update(objects[i], blobs[j])
        end
    end
    for j in 1:m
        if sum(view(matches, :, j)) == 0 && blobs[j].area > min_area
            num_objects += 1
            push!(objects, Object(blobs[j], num_objects))
        end
    end
    return num_objects
end
