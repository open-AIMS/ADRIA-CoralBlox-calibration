# Writes data/unmapped_disturbance_years.csv: every ReefMonitoring "unknown" ("u") disturbance
# on a target reef that carries NO mortality forcing in the datacube, paired with the survey
# year its impact first shows up in.
#
# NOT currently consumed anywhere - the calibration does not down-weight these observations.
# See https://github.com/open-AIMS/ADRIA-CoralBlox-calibration/issues/26 for the proposed
# down-weighting design; this file only identifies the affected (reef, year) cells.
#
# Storms ("s"), COTS ("c") and "multiple" ("m") are deliberately NOT here: they already carry
# forcing via historical_disturbance_mortality_rates.nc. Bleaching ("b") is not here either -
# DHW forcing covers it, so down-weighting bleaching years would penalise the model for a
# disturbance it does see.

using Dates

# Only "u". Disease ("d") is one entry away if it's ever wanted (a single event GBR-wide over
# 2008-2022, so it changes nothing measurable either way).
UNMAPPED_TYPES = ["u"]

# An event this close to a curated row is the same event - the curated `year` column carries
# rounded/annualised dates for non-cyclone disturbances, so exact matching produces false
# "unmapped" hits (verified: Snake Reef's 2021.94 COTS is curated as 2022.0).
MATCH_TOLERANCE_YEARS = 1.5

disturbance_history_cache_path = joinpath(DATA_DIR, "disturbance_history_raw_cache.csv")
if isfile(disturbance_history_cache_path)
    @info "Found disturbance_history_raw_cache file, skipping its generation."
else
    @info "Generating disturbance_history_raw_cache file."
    rm_reef_names, rm_reef_ids = rm_reef_spec(canonical_gpkg)

    history = DataFrame(
        reef_id=String[], reef_name=String[], type=String[], ddate=Float64[],
        description=String[]
    )
    for uid in target_loc_ids
        matched_names = rm_reef_names[rm_reef_ids.==uid]
        isempty(matched_names) && continue
        reef_name = matched_names[1]

        local disturbances
        try
            disturbances = ReefMonitoring.get_disturbances(reef_name)
        catch
            continue
        end
        (isnothing(disturbances) || isempty(disturbances)) && continue

        for row in eachrow(disturbances[disturbances.sample_type.=="MANTA", :])
            ddate = parse(Float64, row.ddate)
            START_YEAR <= ddate <= END_YEAR || continue
            push!(history, (
                uid, reef_name, row.disturbance, ddate, something(row.description, "")
            ))
        end
    end

    open(disturbance_history_cache_path, "w") do io
        println(io, "# Raw per-reef MANTA disturbance history (ALL types, unfiltered), cached from")
        println(io, "# ReefMonitoring.get_disturbances - read by 07_create_unmapped_disturbance_years.jl.")
        println(io, "# Generated: $(now())")
    end
    CSV.write(disturbance_history_cache_path, history; append=true, writeheader=true)
    @info "disturbance_history_raw_cache file successfully generated."
end

# reef_id forced to String - it's an all-digit id, which CSV.jl would otherwise infer as Int64
# and break the joins against target_loc_ids / RME_UNIQUE_ID.
disturbance_history = CSV.read(
    disturbance_history_cache_path, DataFrame;
    stringtype=String, comment="#", types=Dict(:reef_id => String)
)

curated_manta = CSV.read(
    joinpath(DATA_DIR, "disturbance_years_manta.csv"), DataFrame;
    stringtype=String, comment="#", types=Dict(:reef_id => String)
)

"""
Years in which a reef has a non-missing LTMP cover observation. Read from the same gpkg the
objective function scores against, so an emitted `affected_year` is always a cell that exists
in `LocationDataStore.ltmp_coral_cover`. Duplicate `RME_UNIQUE_ID` rows are treated as
observed if ANY of them reports that year, mirroring `build_calibration_data`'s row merge.
"""
function observed_years_by_reef()::Dict{String,Vector{Int}}
    ltmp = GDF.read(LTMP_REEF_DATA_PATH)
    ltmp = ltmp[(!).(ismissing.(ltmp.RME_UNIQUE_ID)), :]

    year_cols = [yr for yr in START_YEAR:END_YEAR if Symbol(yr) in propertynames(ltmp)]

    out = Dict{String,Vector{Int}}()
    for row in eachrow(ltmp)
        uid = string(row.RME_UNIQUE_ID)
        years = [yr for yr in year_cols if !ismissing(row[Symbol(yr)])]
        out[uid] = sort(union(get(out, uid, Int[]), years))
    end
    return out
end

observed_years = observed_years_by_reef()

results = DataFrame(
    reef_id=String[], reef_name=String[], type=String[], ddate=Float64[],
    affected_year=Int[], description=String[]
)

n_no_survey_after = 0
n_already_forced = 0
for row in eachrow(disturbance_history)
    row.type in UNMAPPED_TYPES || continue

    # Already carries datacube forcing => the model sees this disturbance; leave it alone.
    curated_reef = curated_manta[curated_manta.reef_id.==row.reef_id, :]
    if any(abs.(curated_reef.year .- row.ddate) .<= MATCH_TOLERANCE_YEARS)
        global n_already_forced += 1
        continue
    end

    reef_years = get(observed_years, row.reef_id, Int[])
    # The drop can only surface at the first survey taken AFTER the event.
    year_idx = findfirst(y -> y >= row.ddate, reef_years)
    if isnothing(year_idx)
        global n_no_survey_after += 1
        continue
    end

    # Blank descriptions round-trip through the CSV cache as `missing`, not "".
    push!(results, (
        row.reef_id, row.reef_name, row.type, row.ddate, reef_years[year_idx],
        coalesce(row.description, "")
    ))
end

sort!(results, [:reef_id, :ddate])

println("Unmapped-type events found: ", count(in(UNMAPPED_TYPES), disturbance_history.type))
println("  already covered by a curated disturbance: ", n_already_forced)
println("  no LTMP survey after the event: ", n_no_survey_after)
println("  emitted: ", nrow(results),
    " across ", length(unique(results.reef_id)), " reefs, ",
    length(unique(eachrow(results[:, [:reef_id, :affected_year]]))), " distinct (reef, year) cells")

out_path = joinpath(DATA_DIR, "unmapped_disturbance_years.csv")
open(out_path, "w") do io
    println(io, "# ReefMonitoring \"unknown\"-type disturbances with NO mortality forcing in")
    println(io, "# historical_disturbance_mortality_rates.nc. `affected_year` is the first LTMP survey")
    println(io, "# year at or after `ddate`. Not currently consumed by the calibration - see")
    println(io, "# https://github.com/open-AIMS/ADRIA-CoralBlox-calibration/issues/26.")
    println(io, "# Generated: $(now())")
end
CSV.write(out_path, results; append=true, writeheader=true)
println("Written to: ", out_path)
