using ConformalGaussHypergeometric
using Test

@testset verbose = true "ConformalGaussHypergeometric.jl" begin
    @testset "Standard Parameters" begin
        (a,b,c) = (1.1,1.22,1.333)
        @test conformal2F1(a,b,c,.5)         ≈ (2.0054894353193595 + 0.0im)                    rtol=1e-14
        @test conformal2F1(a,b,c,1.25)       ≈ (-4.112910642836204 - 0.0939699849702876im)     rtol=1e-14
        @test conformal2F1(a,b,c,4)          ≈ (-0.34711637240627546 + 0.017257879479033997im) rtol=1e-14
        @test conformal2F1(a,b,c,-1.5)       ≈ (0.39602953929104334 + 0.0im)                   rtol=1e-14
        @test conformal2F1(a,b,c,cispi(1/3)) ≈ (0.4965103388730277 + 0.873467605231226im)      rtol=1e-14
    end

    @testset "Integer Parameters" begin
        @test conformal2F1(1.1,1.2,1.3,1.5+1.5im) ≈ (-0.22147498072051158 + 0.5964515057937457im) rtol=1e-14 skip=true
        @test conformal2F1(-2,-3,-4,.5) ≈ () rtol=1e-14 skip=true
    end
    
    @testset "Error Estimation" begin
        (a,b,c) = (1.1,1.22,1.333)
        @test conformal2F1(a,b,c,.5; esterr=true)[2] < 1e-13
    end
end
