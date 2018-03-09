__precompile__()
module ObjectTracker

export
    Blob,
    form_blobs,
    is_transient,
    match!,
    Object,
    predict,
    State,
    update

include("Kalman.jl")
include("Blobs.jl")
include("Tracker.jl")
end # module
