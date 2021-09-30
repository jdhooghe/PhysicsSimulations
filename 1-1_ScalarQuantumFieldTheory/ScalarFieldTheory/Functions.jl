using Distributions
using Distributed
struct Model
    dt::Float64
    dx::Float64
    Nt::Int64
    Nx::Int64
    m::Float64
    λ::Float64
end

function CreateGrid(Model, Start)
    if Start == "Cold"
        Grid = zeros(Float64, Model.Nt, Model.Nx)
        Grid[Model.Nt, :] = Grid[Model.Nt-1, :]
    else
        Grid = rand(Uniform(0., 1.), (Model.Nt, Model.Nx))
        Grid[Model.Nt, :] = Grid[Model.Nt-1, :]
    end
    return Grid
end

function V(Model::Model, Field::Array{Float64})
    Quad = Field .* Field
    return 0.5 * Model.m * Quad + 0.25 * Model.λ * Quad .* Quad
end

function T(Model::Model, Field::Array{Float64})::Array{Float64}
    First = zeros(Float64, Model.Nt, Model.Nx)
    Second = zeros(Float64, Model.Nt, Model.Nx)

    for i in 1:Model.Nt-1
        for j in 1:Model.Nx-1
            First[i, j] = (Field[i+1, j] - Field[i, j])^2
            Second[i, j] = (Field[i, j+1] - Field[i, j])^2
        end
    end
    for j in 1:Model.Nx-1
        First[Model.Nt, j] = (Field[1, j] - Field[Model.Nt, j])^2
    end
    for i in 1:Model.Nt-1
        Second[i, Model.Nx] = (Field[i, 1] - Field[i, Model.Nx])^2
    end
    First[Model.Nt, Model.Nx] = First[1, 1]
    Second[Model.Nt, Model.Nx] = Second[1, 1]
    a = 1/Model.dt^2
    b = Model.dx^2
    return 0.5*(First*a - Second*b)
end

function GetEnergyAtSite(Time, Space, Model::Model, Field::Array{Float64})::Float64
    Quad = Field[Time, Space] * Field[Time, Space]
    V = 0.5 * Model.m * Quad + 0.25 * Model.λ * Quad * Quad

    if Time == 1
        Tt = (Field[2, Space] - Field[1, Space])^2 + (Field[1, Space] - Field[Model.Nt, Space])^2
    else
        Tt = (Field[Time+1, Space] - Field[Time, Space])^2 + (Field[Time, Space] - Field[Time-1, Space])^2
    end
    Tt = 0.5 * (Tt / Model.dt^2)

    if Space == 1
        Tx = (Field[Time, 2] - Field[Time, 1])^2 + (Field[Time, 1] - Field[Time, Model.Nx])^2
    else
        Tx = (Field[Time, Space+1] - Field[Time, Space])^2 + (Field[Time, Space] - Field[Time, Space-1])^2
    end
    Tx = 0.5 * (Tx / Model.dx^2)

    return V + Tt + Tx

end

function UpdateField(Time::Int64, Space::Int64, Model::Model, Field::Array{Float64}, NewField::Float64, TotalAction::Float64, dAction::Float64)
    TotalAction += dAction
    Field[Time, Space] = NewField
    if Time == 1
        Field[Model.Nt, Space] = NewField
    end
    if Space == 1
        Field[Time, Model.Nx] = NewField
    end
    Field[Model.Nt, Model.Nx] = Field[1, 1]
    return Field, TotalAction
end

function Iterate(Time::Int64, Space::Int64, Step::Float64, Field::Array{Float64}, Model::Model, TotalAction::Float64)

    InitialEnergy = GetEnergyAtSite(Time, Space, Model, Field)

    NewField = Field

    Change = rand(Uniform(-Step, Step))
    TotalField = NewField[Time, Space] + Change

    NewField[Time, Space] += Change
    NewEnergy = GetEnergyAtSite(Time, Space, Model, NewField)

    dAction = (NewEnergy - InitialEnergy) * Model.dt

    if NewEnergy < InitialEnergy
        return UpdateField(Time, Space, Model, Field, Field[Time, Space] + Change, TotalAction, dAction)

    elseif rand(Uniform(0, 1)) < exp(-dAction)
        return UpdateField(Time, Space, Model, Field, Field[Time, Space] + Change, TotalAction, dAction)
    else
        return Field, TotalAction
    end
end
