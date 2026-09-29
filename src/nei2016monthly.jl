export NEI2016MonthlyEmis, NEI2016ElevatedEmis, NEI2016_ELEVATED_SECTORS

# Hourly scale factors, index 1 = local hour 0.
# DIURNAL_FACTORS: HEMCO GEIA_TOD_FOSSIL (CO, formaldehyde), applied by GEOS-Chem at
#   local time UTC + floor(lon/15) h.
# DIURNAL_FACTORS_NOx: profile of HEMCO's EDGAR_TODNOX file, whose built-in shift is
#   round(lon/15) h.
const DIURNAL_FACTORS = [0.45, 0.45, 0.6, 0.6, 0.6, 0.6, 1.45, 1.45, 1.45, 1.45, 1.4, 1.4, 1.4, 1.4, 1.45, 1.45, 1.45, 1.45, 0.65, 0.65, 0.65, 0.65, 0.45, 0.45]
const DIURNAL_FACTORS_NOx = [0.39598674, 0.31852847, 0.30128068, 0.29590213, 0.33177775, 0.43871498, 0.9094625, 1.5850095, 1.6223788, 1.3429453, 1.2265036, 1.1937649, 1.254314, 1.3282939, 1.331211, 1.4135737, 1.6848333, 1.710925, 1.3491899, 1.0586671, 0.84439224, 0.761263, 0.72693235, 0.5741503]
const DIURNAL_FACTORS_ISOP =[0,0,0,0,0,0, 0.2376,0.7224,1.2048,1.656,2.0496,2.3616,2.5728,2.6616,2.6184,2.4408,2.1288,1.6896,1.1448,0.5136, 0,0,0,0]

# Weekday factors, Monday-first (Dates.dayofweek), from HEMCO GEIA_DOW_NOX and GEIA_DOW_CO.
# HEMCO's CO line repeats the NOx weekday value on Tue-Fri and sums to 6.85; 1.1076 is used
# on all weekdays so the week sums to 7.
const DayofWeekFactors_NOx = [1.0706, 1.0706, 1.0706, 1.0706, 1.0706, 0.863, 0.784]
const DayofWeekFactors_CO = [1.1076, 1.1076, 1.1076, 1.1076, 1.1076, 0.779, 0.683]

# Load and create interpolator for delp_dry_surface
const DELP_DRY_SURFACE_ITP = let
    # Load the delp_dry_surface data
    delp_data = load(joinpath(@__DIR__, "mean_domian_delp_dry_surface.jld2"), "mean_domian_delp_dry_surface")

    # Define the grid coordinates
    domain_lon = collect(-125.0:0.625:-66.875)
    domain_lat = collect(25.0:0.5:49.0)

    # Create 2D interpolator with flat extrapolation (use boundary values for out-of-bounds)
    # Note: delp_data should be (lon, lat) ordered to match (domain_lon, domain_lat)
    itp = interpolate((domain_lon, domain_lat), delp_data, Gridded(Linear()))
    extrapolate(itp, Flat())
end

"""
$(SIGNATURES)

Interpolate the delp_dry_surface field at a given longitude and latitude.
Returns the dry pressure thickness value in Pa.
"""
function delp_dry_surface_itp(lon, lat)
    # Convert from radians to degrees if needed
    lon_deg = rad2deg(lon)
    lat_deg = rad2deg(lat)

    # The interpolator now handles out-of-bounds automatically with Flat() extrapolation
    return DELP_DRY_SURFACE_ITP(lon_deg, lat_deg)
end

# Vertical allocation of the NEI "inline" (elevated point) sectors.
#
# The EPA monthly netCDF files are column totals with a single layer. GEOS-Chem reads
# these sectors from 3-D files that its NEI2016 preprocessing (griddedepa2gc.py, B. Henderson)
# built from the same EPA files by multiplying each sector's column total by one
# representative layer profile (after Simpson et al. 2003, EMEP report 1/2003, Table 4.1)
# on the 72-layer GEOS-FP grid, truncated at layer 11 and renormalized.  The same fractions
# are applied here, so the emissions land in the same model layers as in GEOS-Chem.
# Sectors GEOS-Chem reads as 2-D (surface sources, othpt) get profile "surface".
const NEI2016_PROFILE_NAMES = ["surface", "ptegu", "ptnonipm", "pt_oilgas", "cmv"]
const NEI2016_LAYER_FRACTIONS = let
    raw = Dict(
        "surface" => [1.0],
        "ptegu" => [0.027331187, 0.12260809, 0.192017292, 0.212736741, 0.167549199,
            0.110491495, 0.063653198, 0.036979211, 0.022844579, 0.013909341, 0.008781226],
        "ptnonipm" => [0.751204553, 0.144224191, 0.052822078, 0.023840537, 0.011000364,
            0.005661267, 0.003165207, 0.00197344, 0.001300664, 0.000861865, 0.000621856],
        "pt_oilgas" => [0.918700528, 0.059678397, 0.013323237, 0.003978944, 0.001724366,
            0.001038275, 0.000633395, 0.000387233, 0.000243355, 0.000122994, 6.48583e-05],
        "cmv" => [0.537538889, 0.338982929, 0.123478182]
    )
    [raw[p] ./ sum(raw[p]) for p in NEI2016_PROFILE_NAMES]
end

"""
Default mapping from EPA sector name (as in the file name `2016fh_16j_<sector>_12US1_month_MM.ncf`)
to vertical profile. Sectors not listed here are emitted into the surface layer.
"""
const NEI2016_SECTOR_PROFILES = Dict(
    "emln_ptegu" => "ptegu",
    "emln_ptnonipm" => "ptnonipm",
    "emln_ptnonipm_allinln" => "ptnonipm",
    "emln_pt_oilgas" => "pt_oilgas",
    "emln_pt_oilgas_allinln" => "pt_oilgas",
    "emln_cmv_c3_12" => "cmv",
    "emln_cmv_c1c2_12" => "cmv"
)

"""
$(SIGNATURES)

Fraction of a sector's column emissions that is released into model layer `lev`
(1-based layer index) for vertical profile number `profile_id`
(index into `NEI2016_PROFILE_NAMES`). Zero above the top of the profile.
"""
function nei_layer_fraction(profile_id, lev)
    k = round(Int, lev)
    fracs = NEI2016_LAYER_FRACTIONS[round(Int, profile_id)]
    return (1 <= k <= length(fracs)) ? fracs[k] : 0.0
end

"""
$(SIGNATURES)

Dry pressure thickness (hPa, unitless here) of model layer `lev` at the given location.
Layer 1 comes from the stored mean surface map; higher layers are scaled with the
hybrid-grid coefficients `Ap`/`Bp`, using a mean surface pressure backed out of the
layer-1 thickness.
"""
function delp_dry_itp(lon, lat, lev)
    d1 = delp_dry_surface_itp(lon, lat)
    k = round(Int, lev)
    k <= 1 && return d1
    ps = (d1 - NEI_DAP[1]) / NEI_DBP[1]
    return NEI_DAP[k] + NEI_DBP[k] * ps
end

# Per-layer differences of the GEOS-FP hybrid-grid coefficients, precomputed so the
# layer thickness costs two table reads. Ap is in Pa; stored in hPa to match the
# surface thickness map.
const NEI_DAP = [(Ap(k) - Ap(k + 1)) / 100 for k in 1:72]
const NEI_DBP = [Bp(k) - Bp(k + 1) for k in 1:72]

"""
Emission of one species from all NEI profile groups, as a callable parameter:
`(w::NEILayeredEmission)(t, lon, lat, lev)` returns the flux (kg/m²/s) released into
model layer `lev`, divided by that layer's dry thickness (hPa, unitless) and multiplied
by the species' diurnal and day-of-week factors. Above the highest emitting layer it
returns zero without touching the data, so the cost is only paid where there are
emissions. `fracs[g]` are the layer fractions of group `g`, already folded at the domain
top so that no mass is lost when the domain has fewer layers than a profile.
"""
mutable struct NEILayeredEmission{ITP, CF}
    itps::Vector{ITP}
    fracs::Vector{Vector{Float64}}
    ktop::Int
    clock::CF
    "Model grid of the regridded data: first cell centre and spacing (radians) and size."
    grid::NTuple{6, Float64}
    "Dry layer thickness (hPa) at every model grid cell and emitting layer, shared by all species."
    delp::Array{Float64, 3}
end

# Value of a regridded monthly field at model grid cell (i, j), read directly from the
# interpolator's loaded time slices and blended linearly in time. This is what the
# interpolator itself returns at a grid node, without the general interpolation machinery.
# Returns `nothing` if the loaded slices do not cover `t`; the caller then uses the
# interpolator, which also loads the needed data.
@inline function _nei_grid_value(itp::DataSetInterpolator, t, i, j)
    tc = itp.cache
    tc.initialized || return nothing
    times = tc.times
    D = tc.interp_cache
    n = length(times)
    (ndims(D) == 3 && size(D, 3) == n) || return nothing
    @inbounds for s in 1:(n - 1)
        ta = datetime2unix(times[s])
        tb = datetime2unix(times[s + 1])
        if ta <= t < tb
            wt = (t - ta) / (tb - ta)
            return (1 - wt) * D[i, j, s] + wt * D[i, j, s + 1]
        end
    end
    return nothing
end

function (w::NEILayeredEmission)(t, lon, lat, lev)
    k = round(Int, lev)
    T = typeof(float(lon))
    (k < 1 || k > w.ktop) && return zero(T)
    # Model grid cell of this location, if it is exactly on the grid the data were regridded to.
    x0, dx, nx, y0, dy, ny = w.grid
    fi = (lon - x0) / dx + 1
    fj = (lat - y0) / dy + 1
    ongrid = isfinite(fi) && isfinite(fj)
    i = ongrid ? round(Int, fi) : 0
    j = ongrid ? round(Int, fj) : 0
    ongrid = ongrid && abs(fi - i) < 1e-6 && abs(fj - j) < 1e-6 &&
             1 <= i <= nx && 1 <= j <= ny && k <= size(w.delp, 3)
    s = zero(T)
    @inbounds for g in eachindex(w.itps)
        f = w.fracs[g]
        if k <= length(f) && f[k] != 0
            v = ongrid ? _nei_grid_value(w.itps[g], t, i, j) : nothing
            if v === nothing
                v = interp_unsafe(w.itps[g], t, lon, lat)
            end
            s += f[k] * v
        end
    end
    iszero(s) && return s
    d = ongrid ? (@inbounds w.delp[i, j, k]) : delp_dry_itp(lon, lat, k)
    return s * w.clock(t, lon) / d
end

"""
Model grid description `(x0, dx, nx, y0, dy, ny)` of the interpolator's regridded data,
and whether direct grid reads are possible (x then y as the first two data dimensions).
"""
function nei_model_grid(itp::DataSetInterpolator)
    coords = _model_grid(itp)
    ok = itp.metadata.xdim == 1 && itp.metadata.ydim == 2 && length(coords) == 2
    xs, ys = collect(coords[1]), collect(coords[2])
    step(v) = length(v) > 1 ? (v[end] - v[1]) / (length(v) - 1) : 1.0
    ok = ok && length(xs) > 1 && length(ys) > 1 &&
         maximum(abs.(diff(xs) .- step(xs))) < 1e-9 && maximum(abs.(diff(ys) .- step(ys))) < 1e-9
    grid = ok ? (Float64(xs[1]), Float64(step(xs)), Float64(length(xs)),
                 Float64(ys[1]), Float64(step(ys)), Float64(length(ys))) :
           (0.0, 1.0, 0.0, 0.0, 1.0, 0.0) # nx = ny = 0: never on grid
    return grid, xs, ys
end

is_itp_wrapper(::NEILayeredEmission) = true
wrapped_itps(w::NEILayeredEmission) = w.itps
function lazyload_wrapper!(w::NEILayeredEmission, t)
    for i in eachindex(w.itps)
        w.itps[i] = lazyload!(w.itps[i], t)
    end
    w
end

# Diurnal and day-of-week factors per species (functions of unix time and longitude).
nei_clock_CO(t, lon) = dayofweek_itp_CO(t, lon) * diurnal_itp(t, lon)
nei_clock_NOx(t, lon) = dayofweek_itp_NOx(t, lon) * diurnal_itp_NOx(t, lon)
nei_clock_none(t, lon) = 1.0
function nei_clock(varname)
    varname == "CO" && return nei_clock_CO
    varname == "FORM" && return diurnal_itp
    varname == "ISOP" && return diurnal_itp_ISOP
    varname in ("NO", "NO2") && return nei_clock_NOx
    return nei_clock_none
end

"""
$(SIGNATURES)

Layer fractions of profile `profile_id` for a domain whose top layer is `ktop`: the
fractions above the domain top are added to the top layer, so the column total is kept.
"""
function nei_domain_fractions(profile_id, ktop)
    f = NEI2016_LAYER_FRACTIONS[profile_id]
    length(f) <= ktop && return copy(f)
    g = f[1:ktop]
    g[ktop] += sum(@view f[(ktop + 1):end])
    return g
end

"""
$(SIGNATURES)

Diurnal interpolation function that returns the scale factor for a given time.
Returns different emission scaling factors based on the hour of day.
"""
function diurnal_itp(t, lon)
    ut = Dates.unix2datetime(t)

    # Convert radians to degrees for timezone calculation
    lon_deg = rad2deg(lon)
    dt = floor(lon_deg / 15) # hours; HEMCO local-time rule
    t_local = t + dt * 3600 # in seconds
    ut_local = Dates.unix2datetime(t_local)
    hour_of_day = Dates.hour(ut_local) + 1  # +1 for 1-based indexing

    return DIURNAL_FACTORS[hour_of_day]
end

function diurnal_itp_NOx(t, lon)
    ut = Dates.unix2datetime(t)

    # Convert radians to degrees for timezone calculation
    lon_deg = rad2deg(lon)
    dt = round(lon_deg / 15) # hours; solar time zone as in GEOS-Chem's EDGAR NOx file
    t_local = t + dt * 3600 # in seconds
    ut_local = Dates.unix2datetime(t_local)
    hour_of_day = Dates.hour(ut_local) + 1  # +1 for 1-based indexing

    return DIURNAL_FACTORS_NOx[hour_of_day]
end

function diurnal_itp_ISOP(t, lon)
    ut = Dates.unix2datetime(t)

    # Convert radians to degrees for timezone calculation
    lon_deg = rad2deg(lon)
    dt = round(lon_deg / 15) # hours; same clock as NOx, no GEOS-Chem counterpart
    t_local = t + dt * 3600 # in seconds
    ut_local = Dates.unix2datetime(t_local)
    hour_of_day = Dates.hour(ut_local) + 1  # +1 for 1-based indexing

    return DIURNAL_FACTORS_ISOP[hour_of_day]
end

"""
$(SIGNATURES)

Day of week interpolation function that returns the scale factor for a given time.
Returns different emission scaling factors based on the day of week.
"""
function dayofweek_itp_CO(t, lon)
    ut = Dates.unix2datetime(t)

    # Convert radians to degrees for timezone calculation
    lon_deg = rad2deg(lon)
    dt = floor(lon_deg / 15) # hours; HEMCO local-time rule
    t_local = t + dt * 3600 # in seconds
    ut_local = Dates.unix2datetime(t_local)
    day_of_week = Dates.dayofweek(ut_local)

    return DayofWeekFactors_CO[day_of_week]
end

function dayofweek_itp_NOx(t, lon)
    ut = Dates.unix2datetime(t)

    # Convert radians to degrees for timezone calculation
    lon_deg = rad2deg(lon)
    dt = round(lon_deg / 15) # hours; solar time zone as in GEOS-Chem's EDGAR NOx file
    t_local = t + dt * 3600 # in seconds
    ut_local = Dates.unix2datetime(t_local)
    day_of_week = Dates.dayofweek(ut_local)

    return DayofWeekFactors_NOx[day_of_week]
end

# Register the symbolic function
@register_symbolic diurnal_itp(t, lon)
@register_symbolic diurnal_itp_NOx(t, lon)
@register_symbolic diurnal_itp_ISOP(t, lon)
@register_symbolic dayofweek_itp_CO(t, lon)
@register_symbolic dayofweek_itp_NOx(t, lon)
@register_symbolic delp_dry_surface_itp(lon, lat)

# Dummy function for unit validation. ModelingToolkit will call this function
# with a DynamicQuantities.Quantity to get information about the type and units of the output.
diurnal_itp(t::DynamicQuantities.Quantity, lon) = 1.0
diurnal_itp_NOx(t::DynamicQuantities.Quantity, lon) = 1.0
diurnal_itp_ISOP(t::DynamicQuantities.Quantity, lon) = 1.0
dayofweek_itp_CO(t::DynamicQuantities.Quantity, lon) = 1.0
dayofweek_itp_NOx(t::DynamicQuantities.Quantity, lon) = 1.0
delp_dry_surface_itp(lon::DynamicQuantities.Quantity, lat::DynamicQuantities.Quantity) = 1.0

"""
$(SIGNATURES)

Archived CMAQ emissions data.

Currently, only data for year 2016 is available.
"""
struct NEI2016MonthlyEmisFileSet <: FileSet
    mirror::AbstractString
    sector::Any
    ds::Any
    freq_info::DataFrequencyInfo
    function NEI2016MonthlyEmisFileSet(sector, starttime, endtime)
        NEI2016MonthlyEmisFileSet("https://gaftp.epa.gov/Air/", sector, starttime, endtime)
    end
    function NEI2016MonthlyEmisFileSet(mirror, sector, starttime, endtime)
        floormonth(t) = DateTime(Dates.year(t), Dates.month(t))
        check_times = (floormonth(starttime - Day(16))):Month(1):(endtime + Day(16))
        fs = new(mirror, sector, nothing, DataFrequencyInfo(starttime, Day(1), check_times))
        filepaths = maybedownload.((fs,), check_times)

        start = floormonth(starttime)
        frequency = ((start + Dates.Month(1)) - start) # Only true for the first month.
        centerpoints = [t + Second(t + Month(1) - t) / 2 for t in check_times]
        dfi = DataFrequencyInfo(start, frequency, centerpoints)

        lock(nclock) do
            ds = NCDataset(filepaths, aggdim = "TSTEP")
            new(mirror, sector, ds, dfi)
        end
    end
end

"""
$(SIGNATURES)

File path on the server relative to the host root; also path on local disk relative to `ENV["EARTHSCIDATADIR"]`.
"""
function relpath(fs::NEI2016MonthlyEmisFileSet, t::DateTime)
    @assert Dates.year(t)==2016 "Only 2016 emissions data is available with `NEI2016MonthlyEmis`."
    month = lpad(Dates.month(t), 2, '0')
    return "emismod/2016/v1/gridded/monthly_netCDF/2016fh_16j_$(fs.sector)_12US1_month_$(month).ncf"
end

DataFrequencyInfo(fs::NEI2016MonthlyEmisFileSet) = fs.freq_info

"""
$(SIGNATURES)

Load the NEI data for the given variable name at the given time.
This loads data in kg/s/m^2 units on the NEI source grid for regridding.
"""
function loadslice!(
        data::AbstractArray,
        fs::NEI2016MonthlyEmisFileSet,
        t::DateTime,
        varname
)
    lock(nclock) do
        data = reshape(data, size(data)..., 1)
        var = loadslice!(data, fs, fs.ds, t, varname, "TSTEP")

        # Step 1: Apply unit conversion from the file (typically tons/day to kg/s)
        scale, _ = to_unit(var.attrib["units"])
        if scale != 1
            data .*= scale  # Now data is in kg/s per grid cell
        end

        # The EPA files are labeled "tons/day" but the stored values are actually
        # monthly totals. The unit conversion above divides by one day (86400 s),
        # so applying the value as a constant rate over the whole month over-counts
        # emissions by the number of days in the month (~30x). Normalize to a true
        # daily rate so the per-second rate is correct.
        data ./= Dates.daysinmonth(t)

        # Step 2: Convert from kg/s per grid cell to kg/s/m² for conservative regridding
        # This is the flux density that can be conservatively regridded
        Δx = fs.ds.attrib["XCELL"]  # Cell width in meters
        Δy = fs.ds.attrib["YCELL"]  # Cell height in meters
        data ./= (Δx * Δy)  # Now data is in kg/s/m²
    end
    nothing
end

"""
$(SIGNATURES)

Load the data for the given variable name at the given time.
"""
function loadmetadata(fs::NEI2016MonthlyEmisFileSet, varname)::MetaData
    lock(nclock) do
        timedim = "TSTEP"
        var = fs.ds[varname]
        dims = collect(NCDatasets.dimnames(var))
        @assert timedim ∈ dims "Variable $varname does not have a dimension named '$timedim'."
        time_index = findfirst(isequal(timedim), dims)
        dims = deleteat!(dims, time_index)
        varsize = deleteat!(collect(size(var)), time_index)
        @assert varsize[end]==1 "Only 2D data is supported."
        varsize = varsize[1:(end - 1)] # Last dimension is 1.

        Δx = fs.ds.attrib["XCELL"]
        Δy = fs.ds.attrib["YCELL"]
        _, units = to_unit(var.attrib["units"])
        units /= u"m^2"
        description = var.attrib["var_desc"]

        x₀ = fs.ds.attrib["XORIG"]
        y₀ = fs.ds.attrib["YORIG"]
        Δx = fs.ds.attrib["XCELL"]
        Δy = fs.ds.attrib["YCELL"]
        nx = fs.ds.attrib["NCOLS"]
        ny = fs.ds.attrib["NROWS"]
        xs = x₀ + Δx / 2 .+ Δx .* (0:(nx - 1))
        ys = y₀ + Δy / 2 .+ Δy .* (0:(ny - 1))

        coords = [xs, ys]

        p_alp = fs.ds.attrib["P_ALP"]
        p_bet = fs.ds.attrib["P_BET"]
        #p_gam = fs.ds.attrib["P_GAM"] # Don't think this is used for anything.
        x_cent = fs.ds.attrib["XCENT"]
        y_cent = fs.ds.attrib["YCENT"]
        native_sr = "+proj=lcc +lat_1=$(p_alp) +lat_2=$(p_bet) +lat_0=$(y_cent) +lon_0=$(x_cent) +x_0=0 +y_0=0 +a=6370997.000000 +b=6370997.000000 +to_meter=1"

        xdim = findfirst((x) -> x == "COL", dims)
        ydim = findfirst((x) -> x == "ROW", dims)
        @assert xdim>0 "NEI2016 `COL` dimension not found"
        @assert ydim>0 "NEI2016 `ROW` dimension not found"

        return MetaData(
            coords,
            string(units),
            description,
            dims,
            varsize,
            native_sr,
            xdim,
            ydim,
            -1,
            (false, false, false)
        )
    end
end

function get_geometry(fs::NEI2016MonthlyEmisFileSet, m::MetaData)
    x₀,y₀,Δx,Δy,nx,ny = lock(nclock) do
        x₀ = fs.ds.attrib["XORIG"]
        y₀ = fs.ds.attrib["YORIG"]
        Δx = fs.ds.attrib["XCELL"]
        Δy = fs.ds.attrib["YCELL"]
        nx = fs.ds.attrib["NCOLS"]
        ny = fs.ds.attrib["NROWS"]
        x₀,y₀,Δx,Δy,nx,ny
    end
    # Create edges (nx+1 and ny+1 points) so we get nx*ny cells
    x = range(start=x₀, step=Δx, length=nx+1)
    y = range(start=y₀, step=Δy, length=ny+1)
    # Use column-major (x-fastest) ordering to match vec() on data arrays
    polys = Vector{Vector{NTuple{2, Float64}}}(undef, nx*ny)
    for j in 1:ny, i in 1:nx
        polys[(j-1)*nx + i] = [(x[i], y[j]), (x[i+1], y[j]), (x[i+1], y[j+1]),
            (x[i], y[j+1]), (x[i], y[j])]
    end
    return polys
end

"""
$(SIGNATURES)

Return the variable names associated with this FileSet.
"""
function varnames(fs::NEI2016MonthlyEmisFileSet)
    lock(nclock) do
        return [setdiff(keys(fs.ds), ["TFLAG"; keys(fs.ds.dim)])...]
    end
end

Base.close(fs::NEI2016MonthlyEmisFileSet) = lock(nclock) do; close(fs.ds); end

"""
$(SIGNATURES)

Several NEI2016 sectors on the same grid, summed at load time. A variable that is
missing from one sector's file is treated as zero for that sector.
"""
struct NEI2016MonthlyEmisMultiFileSet <: FileSet
    filesets::Vector{NEI2016MonthlyEmisFileSet}
    function NEI2016MonthlyEmisMultiFileSet(
            sectors::AbstractVector{<:AbstractString}, starttime, endtime)
        isempty(sectors) && error("At least one NEI2016 sector must be given.")
        fss = [NEI2016MonthlyEmisFileSet(String(s), starttime, endtime) for s in sectors]
        gridkeys = ("GDNAM", "XORIG", "YORIG", "XCELL", "YCELL", "NCOLS", "NROWS",
            "P_ALP", "P_BET", "XCENT", "YCENT")
        grid(fs) = lock(nclock) do
            [fs.ds.attrib[k] for k in gridkeys]
        end
        g1 = grid(fss[1])
        for fs in fss[2:end]
            grid(fs) == g1 ||
                error("NEI2016 sector $(fs.sector) is not on the same grid as $(fss[1].sector).")
        end
        new(fss)
    end
end

function _first_with(fs::NEI2016MonthlyEmisMultiFileSet, varname)
    i = findfirst(sub -> varname in varnames(sub), fs.filesets)
    isnothing(i) && error("Variable $varname not found in any of the NEI2016 sectors.")
    fs.filesets[i]
end

mirror(fs::NEI2016MonthlyEmisMultiFileSet) = mirror(first(fs.filesets))
relpath(fs::NEI2016MonthlyEmisMultiFileSet, t::DateTime) = relpath(first(fs.filesets), t)
DataFrequencyInfo(fs::NEI2016MonthlyEmisMultiFileSet) = DataFrequencyInfo(first(fs.filesets))
function loadmetadata(fs::NEI2016MonthlyEmisMultiFileSet, varname)::MetaData
    loadmetadata(_first_with(fs, varname), varname)
end
get_geometry(fs::NEI2016MonthlyEmisMultiFileSet, m::MetaData) = get_geometry(first(fs.filesets), m)
varnames(fs::NEI2016MonthlyEmisMultiFileSet) = unique(vcat(varnames.(fs.filesets)...))

function loadslice!(data::AbstractArray, fs::NEI2016MonthlyEmisMultiFileSet, t::DateTime, varname)
    fill!(data, zero(eltype(data)))
    buf = similar(data)
    for sub in fs.filesets
        varname in varnames(sub) || continue
        loadslice!(buf, sub, t, varname)
        data .+= buf
    end
    nothing
end

Base.close(fs::NEI2016MonthlyEmisMultiFileSet) = foreach(close, fs.filesets)

struct NEI2016MonthlyEmisCoupler
    sys::Any
end

"""
$(SIGNATURES)

A data loader for CMAQ-formatted monthly US National Emissions Inventory data for year 2016,
available from: https://gaftp.epa.gov/Air/emismod/2016/v1/gridded/monthly_netCDF/.
The emissions here are monthly averages, so there is no information about diurnal variation etc.

`sectors` is one sector name or a vector of sector names as they appear in the file names
`2016fh_16j_<sector>_12US1_month_MM.ncf`, e.g. `"mrggrid_withbeis_withrwc"` (all surface
sectors merged) or `"emln_ptegu"` (electricity generating units). Emissions from all listed
sectors are summed. The merged surface file does not contain the "inline" elevated point
sectors (`emln_ptegu`, `emln_ptnonipm`, `emln_pt_oilgas`, `emln_othpt`, `emln_cmv_c3_12`,
`emln_cmv_c1c2_12`, and the fire sectors); add them by listing them here.

Elevated sectors are spread over model layers with the same per-sector layer profiles that
GEOS-Chem uses for its NEI2016 3-D input files (see `NEI2016_SECTOR_PROFILES` and
`NEI2016_LAYER_FRACTIONS`). Sectors not listed in `vertical_profiles` are emitted into the
surface layer. Pass a different `vertical_profiles` dictionary (sector name => profile name,
one of $(NEI2016_PROFILE_NAMES)) to override the defaults.

The emissions are returned as mixing ratios in units of kg/kg/s by converting from the
native flux density (kg/m²/s) using:

    mixing_ratio = flux / (g0_100 * delp_dry)

where g0_100 ≈ 10.197 kg/m² and delp_dry is the dry pressure thickness of the model layer
(physically in hPa, but unitless here) that varies spatially across the domain.

`scale` is a scaling factor to apply to the emissions data. The default value is 1.0.

`stream` specifies whether the data should be streamed in as needed or loaded all at once.

Conservative regridding (via ConservativeRegridding.jl) is used by default to map emissions
from the native NEI Lambert Conformal Conic grid to the simulation domain grid, preserving
total emissions mass.
"""
function NEI2016MonthlyEmis(
        sectors::AbstractVector{<:AbstractString},
        domaininfo::DomainInfo;
        scale = 1.0,
        name = :NEI2016MonthlyEmis,
        stream = true,
        vertical_profiles::AbstractDict = NEI2016_SECTOR_PROFILES
)
    starttime, endtime = get_tspan_datetime(domaininfo)
    for s in sectors
        p = get(vertical_profiles, s, "surface")
        p in NEI2016_PROFILE_NAMES ||
            error("Unknown vertical profile \"$p\" for sector $s; use one of $(NEI2016_PROFILE_NAMES).")
    end
    # Group the sectors by vertical profile. Each group is summed at load time, and all
    # groups share one regridder because every sector is on the same 12US1 grid.
    groups = Tuple{Int, FileSetWithRegridder}[]
    regrid_f = nothing
    for (pid, pname) in enumerate(NEI2016_PROFILE_NAMES)
        secs = [String(s) for s in sectors if get(vertical_profiles, s, "surface") == pname]
        isempty(secs) && continue
        _fs = NEI2016MonthlyEmisMultiFileSet(secs, starttime, endtime)
        if isnothing(regrid_f)
            regrid_f = regridder(_fs, loadmetadata(_fs, first(varnames(_fs))), domaininfo)
        end
        push!(groups, (pid, FileSetWithRegridder(_fs, regrid_f)))
    end
    pvdict = Dict([Symbol(v) => v for v in EarthSciMLBase.pvars(domaininfo)]...)
    @assert :x in keys(pvdict)||:lon in keys(pvdict) "x or lon must be specified in the domaininfo"
    @assert :y in keys(pvdict)||:lat in keys(pvdict) "y or lat must be specified in the domaininfo"
    @assert :lev in keys(pvdict) "lev must be specified in the domaininfo"
    x = :x in keys(pvdict) ? pvdict[:x] : pvdict[:lon]
    y = :y in keys(pvdict) ? pvdict[:y] : pvdict[:lat]
    lev = pvdict[:lev]
    @parameters(Δz=1.0,
        [description = "Couldn't remove Δz without getting errors, so I set it to 1.0 without units"],)
    @parameters t_ref = get_tref(domaininfo) [unit = u"s", description = "Reference time"]
    # Conversion constant: g0_100 = 100 hPa / g0 where g0 = 9.80665 m/s²
    @parameters g0_100 = 100.0 / 9.80665 [unit = u"kg/m^2"]
    eqs = Equation[]
    params = Any[t_ref, g0_100]
    vars = Num[]

    dt = EarthSciMLBase.eltype(domaininfo)
    # Top model layer of the domain; profiles reaching above it are folded into it.
    levidx = findfirst(v -> Symbol(v) === Symbol(lev), EarthSciMLBase.pvars(domaininfo))
    ktop = round(Int, maximum(EarthSciMLBase.grid(domaininfo, (false, false, false))[levidx]))
    allvars = unique(vcat([varnames(fs.fs) for (_, fs) in groups]...))
    # Highest emitting layer of any group in this domain.
    kemit = min(ktop, maximum(length(NEI2016_LAYER_FRACTIONS[pid]) for (pid, _) in groups))
    grid = nothing
    delp_table = Array{Float64, 3}(undef, 0, 0, 0)
    for varname in allvars
        # One interpolator per vertical-profile group, all held by one callable parameter
        # that sums (group column flux) x (group layer fraction) over the groups emitting
        # into layer `lev`, divides by the layer's dry thickness, and applies the species'
        # diurnal and day-of-week factors. Above the emitting layers it returns zero
        # without evaluating anything else.
        itps = DataSetInterpolator[]
        fracs = Vector{Float64}[]
        for (pid, fs) in groups
            varname in varnames(fs.fs) || continue
            push!(itps, DataSetInterpolator{dt}(fs, varname, starttime, endtime, domaininfo;
                stream = stream))
            push!(fracs, nei_domain_fractions(pid, ktop))
        end
        itps = [itps...] # concrete element type
        if grid === nothing
            # Model grid and layer-thickness table, computed once and shared by all species.
            grid, xs, ys = nei_model_grid(first(itps))
            delp_table = grid[3] > 0 ?
                         [delp_dry_itp(xs[i], ys[j], k)
                          for i in eachindex(xs), j in eachindex(ys), k in 1:kemit] :
                         Array{Float64, 3}(undef, 0, 0, 0)
        end
        w = NEILayeredEmission(itps, fracs, min(ktop, maximum(length.(fracs))),
            nei_clock(varname), grid, delp_table)
        n_p = Symbol(varname, "_itp")
        T_w = typeof(w)
        p_w = only(@parameters ($n_p::T_w)(..) = w [
            unit = units(first(itps)),
            description = "Layered NEI emission of $(varname)"
        ])
        push!(params, p_w)
        desc = description(first(itps))
        staggering = first(itps).metadata.staggering

        # Mixing ratio conversion: mixing_ratio = flux / (g0_100 * delp_dry(x, y, lev));
        # the division by the layer thickness happens inside the callable.
        rhs = p_w(t_ref + t, x, y, lev) / Δz * scale / g0_100

        n = Symbol(varname)
        uu = ModelingToolkit.get_unit(rhs)
        lhs = only(@variables $n(t) [
            unit = uu,
            description = desc,
            misc = Dict(:staggering => staggering)
        ])
        push!(eqs, lhs ~ rhs)
        push!(vars, lhs)
    end
    all_params = [x, y, lev, Δz, params...]
    sys = System(
        eqs,
        t,
        vars,
        all_params;
        name = name,
        initial_conditions = _itp_defaults(all_params),
        metadata = Dict(CoupleType => NEI2016MonthlyEmisCoupler,
            SysDiscreteEvent => create_updater_sys_event(name, params, starttime))
    )
    return sys
end

function NEI2016MonthlyEmis(sector::AbstractString, domaininfo::DomainInfo; kwargs...)
    NEI2016MonthlyEmis([sector], domaininfo; kwargs...)
end

"""
The six "inline" elevated point sectors of the 2016 NEI platform that GEOS-Chem reads in
addition to the surface sectors: power plants, other industrial points, oil and gas points
(the inline part only; the low-level part is in the merged surface file), Canada and Mexico
points, and class 1/2 and class 3 commercial marine vessels.
"""
const NEI2016_ELEVATED_SECTORS = ["emln_ptegu", "emln_ptnonipm", "emln_pt_oilgas",
    "emln_othpt", "emln_cmv_c3_12", "emln_cmv_c1c2_12"]

"""
$(SIGNATURES)

The NEI 2016 elevated point sectors ([`NEI2016_ELEVATED_SECTORS`](@ref)) as one emission
system, vertically allocated as in GEOS-Chem. Equivalent to
`NEI2016MonthlyEmis(NEI2016_ELEVATED_SECTORS, domaininfo; kwargs...)`; to combine them with
the merged surface file in a single system, pass both to [`NEI2016MonthlyEmis`](@ref) instead:

    NEI2016MonthlyEmis(["mrggrid_withbeis_withrwc"; NEI2016_ELEVATED_SECTORS], domaininfo)
"""
function NEI2016ElevatedEmis(domaininfo::DomainInfo; name = :NEI2016ElevatedEmis, kwargs...)
    NEI2016MonthlyEmis(NEI2016_ELEVATED_SECTORS, domaininfo; name = name, kwargs...)
end

