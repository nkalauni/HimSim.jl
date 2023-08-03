export ReadFromCSV, hargreaves

"Wrapper function to read from csv"

function ReadFromCSV(filepath; format="mm/dd/yyyy", dateCol=:datetime)
    # return DataFrame(CSV.File(filepath))
    df = CSV.read(filepath, DataFrame)
    dfmat = DateFormat(format)
    tCol = df[!, dateCol]
    len = length(tCol)
    dts = Array{Date}(undef, len)
    for i in 1:len
        dts[i] = Date(tCol[i], dfmat)
    end
    df[!, dateCol] = dts
    return df
end

"Calculate PET using Hargreaves-Samani equation"
function hargreaves(forcings, latitude; tminCol=:tmin, tmaxCol=:tmax, dtCol=:datetime)

    dts = forcings[!, dtCol]
    tmin = forcings[!, tminCol]
    tmax = forcings[!, tmaxCol]
    len = length(tmax)
    Gsc = 0.0820
    latRad = latitude * (pi / 180)

    doy = map(dayofyear, dts)

    tavg = (tmin + tmax) / 2

    eto = zeros(len)

    for (i, t) in enumerate(doy)
        dr = 1 + 0.33 * cos((2 * pi * t) / 365)
        delta = 0.409 * sin((2 * pi * t) / 365 - 1.39)
        ws = acos(-tan(latRad) * tan(delta))
        Ra = (24 * 60) / pi * Gsc * dr * (ws * sin(latRad) * sin(delta) + cos(latRad) * cos(delta) * sin(ws))
        eto[i] = (0.023 * 0.408 * (tavg[i] + 17.8) * (tmax[i] - tmin[i])^0.5 * Ra)
    end

    return eto
end
