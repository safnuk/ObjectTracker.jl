A = zeros(UInt8, 10, 10, 10)
for i in 1:8
    A[i, 5, i:i+2] = 0x01
end
A[9:10, 5, 8:10] = 0x01

blob_series = form_blobs(A)
@test length(blob_series) == 10
for blobs in blob_series
    @test length(blobs) == 1
    blob = blobs[1]
    @test blob.centroid[1] == 5
    @test blob.ratio ≈ 3.0
    @test blob.area == 3
    @test blob.source == 1
    @test blob.θ ≈ 0.0
end

