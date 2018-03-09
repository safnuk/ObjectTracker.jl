const TRANSITION = [1.0 1.0; 0.0 1.0]
const POS_NOISE = 1.5
const NOISE = [POS_NOISE 0.0; 0.0 sqrt(2)*POS_NOISE]
const PREDICT_NOISE = [0.1 0.0; 0.0 0.3]

struct State
    p::Float64
    v::Float64
    Σ::Array{Float64, 2}
end
State(p, v) = State(p, v, NOISE)

function predict(s::State, noise=PREDICT_NOISE)
    p, v = TRANSITION * [s.p ; s.v]
    Σ = TRANSITION * s.Σ * TRANSITION' + noise
    State(p, v, Σ)
end

function update(predicted::State, observed::State)
    K = predicted.Σ / (predicted.Σ + observed.Σ)
    # x = [predicted.p ; predicted.v]
    # z = [observed.p; observed.v]
    # p, v = x + K * (z - x)
    p = observed.p
    v = observed.v
    Σ = predicted.Σ - K * predicted.Σ
    State(p, v, Σ)
end

function smoothed_velocity_estimate(old::State, new::State, λ)
    diff = new.p - old.p
    scale = λ
    #scale = 1.0 - 0.5 * (1 - λ) * (1 + 1.0 / sqrt(2π * old.Σ[2, 2]^2))
    return scale * diff + (1 - scale) * old.v
end

function pdf(s::State, p::Float64)
    σ = s.Σ[1, 1]
    exp(-((s.p - p)^2) / (2σ^2)) / (sqrt(2π * σ^2))
end
