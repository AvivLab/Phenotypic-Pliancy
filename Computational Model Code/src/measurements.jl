function genPopMeasurements()
    indMeasures = genMeasures()
    Measurements(indMeasures, Array{Union{Float64,Missing}}(missing, 2, 3),
     Array{Union{Float64,Missing}}(missing, 2, 3),
      zeros(Int64, 3), zeros(Int64, 3), 0., 0., 0.,
       zeros(Float64, 2, 2), zeros(Float64, 2, 2), zeros(Float64, 2, 2),
        zeros(Float64, 2), zeros(Float64, 2), zeros(Float64, 2),
         zeros(Float64, 2, 2), zeros(Float64, 2, 2), zeros(Float64, 2, 2),
          zeros(Float64, 2), zeros(Float64, 2, 2), zeros(Float64, 2),
           zeros(Float64, 2), zeros(Float64, 2), zeros(Float64, 2, 3))
end
