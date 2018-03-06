A = zeros(UInt8, 10, 10, 10)
for i in 1:8
    A[i, 5, i:i+2] = 0x01
end
A[9:10, 5, 8:10] = 0x01
objects = Vector{Object}()

blob_series = form_blobs(A)
for blobs in blob_series
    match!(objects, blobs, 1; min_area=2)
    @test length(objects) == 1
end
