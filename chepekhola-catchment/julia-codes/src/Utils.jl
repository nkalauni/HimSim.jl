module Utils

using DataFrames, Dates, CSV

"Wrapper function to read from csv"
function readfromcsv(filepath)
    return DataFrame(CSV.file(filepath))
end

"Caulate PET using Hargreaves-Samani equation"
function hargreaves(forcings, latitude; tminCol=:tmin, tmaxCol=:tmax, dtCol=:datetime)

    dts = forcings[!, dtCol]
    tmin = forcings[!, tminCol]
    tmax = forcings[!, tmaxCol]
    len = length(tmax)
    Gsc = 0.0820
    latRad = latitude * (pi / 180)

    doy = map(dayofyear, dts)

    tavg = map(mean, zip(tmin, tmax))

    eto = zeros(len)

    for (i, t) in enumerate(doy)
        dr = 1 + 0.33 * cos((2 * pi * t) / 365)
        delta = 0.409 * sin((2 * pi * t) / 365 - 1.39)
        ws = acos(-tan(latRad) * tan(delta))
        Ra = (24 * 60) / pi * Gsc * dr * (ws * sin(latRad) * sin(delta) + cos(latRad) * cos(delta) * sin(ws))
        eto[i] = (0.023 * 0.408 * (tavg[i] + 17.8) * (tmax - tmin)^0.5 * Ra)
    end

    return eto
end

end
