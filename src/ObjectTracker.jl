__precompile__()
module ObjectTracker

export
    Blob,
    form_blobs,
    is_transient,
    match!,
    Object,
    update

include("Blobs.jl")
include("Tracker.jl")
end # module
