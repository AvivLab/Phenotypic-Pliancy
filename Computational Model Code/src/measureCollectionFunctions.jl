# Functions:
function setDfValueValueVector(df::DataFrame, value::Any, genIndex::Int64, gensSaveData::Vector)
    for i = 1:length(value)
        df[df.Gen .== gensSaveData[genIndex],3+i] = [value[i]]
    end
end

function setDfValueValueMatrix(df::DataFrame, avgValue::Any, stdValue::Any, genIndex::Int64, gensSaveData::Vector)
    numCols = size(df)[2]
    avgCols = [(numCols-(numCols-4)):2:(numCols-1);]
    stdCols = [(numCols-(numCols-4-1)):2:numCols;]
    for i = 1:length(avgValue)
        df[df.Gen .== gensSaveData[genIndex],avgCols[i]] = [avgValue[i]]
        df[df.Gen .== gensSaveData[genIndex],stdCols[i]] = [stdValue[i]]
    end
end

function distancePopAvgs(popMeasurements::Measurements, measurementName::String)
  measurementVals = map(x -> getproperty(popMeasurements.individuals[x], Symbol(measurementName)), 1:N)
  sizeOfMeas = size(measurementVals[1])
  mapNum = sizeOfMeas[1]*sizeOfMeas[2]
  avg = round.(map(y-> mean(skipmissing(map(x -> measurementVals[x][y], 1:N))), 1:mapNum), digits = 3)
  return reshape(avg, sizeOfMeas[1], sizeOfMeas[2])
end

function distancePopStd(popMeasurements::Measurements, measurementName::String)
  measurementVals = map(x -> getproperty(popMeasurements.individuals[x], Symbol(measurementName)), 1:N)
  sizeOfMeas = size(measurementVals[1])
  mapNum = sizeOfMeas[1]*sizeOfMeas[2]
  standdev = round.(map(y-> std(skipmissing(map(x -> measurementVals[x][y], 1:N))), 1:mapNum), digits = 3)
  return reshape(standdev, sizeOfMeas[1], sizeOfMeas[2])
end
