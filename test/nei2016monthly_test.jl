@testitem "NEI Setup" begin
    using Dates: DateTime
    using EarthSciMLBase
    using EarthSciData

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 2);
        latrange = deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
        lonrange = deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
        levrange = 1:10
    )
    lon, lat, lev = EarthSciMLBase.pvars(domain)

    ts, te = get_tspan_datetime(domain)
    sample_time = ts

    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    fileset = EarthSciData.NEI2016MonthlyEmisFileSet("mrggrid_withbeis_withrwc", ts, te)
end

@testsnippet NEISetup begin
    using Dates: DateTime
    using EarthSciMLBase
    using EarthSciData

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 2);
        latrange = deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
        lonrange = deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
        levrange = 1:10
    )
    lon, lat, lev = EarthSciMLBase.pvars(domain)

    ts, te = get_tspan_datetime(domain)
    sample_time = ts

    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    fileset = EarthSciData.NEI2016MonthlyEmisFileSet("mrggrid_withbeis_withrwc", ts, te)
end

@testitem "NEI Basics" setup=[NEISetup] begin
    using ModelingToolkit: equations
    eqs = equations(emis)
    @test length(eqs) == 69
end

@testitem "projections" setup=[NEISetup] begin
    import Proj
    fs = EarthSciData.FileSetWithRegridder(fileset, EarthSciData.regridder(fileset,
        EarthSciData.loadmetadata(fileset, first(EarthSciData.varnames(fileset))), domain))
    @testset "correct projection" begin
        itp = EarthSciData.DataSetInterpolator{Float32}(fs, "NOX", ts, te, domain)
        result = interp!(itp, sample_time, deg2rad(-97.0f0), deg2rad(40.0f0))
        @test result > 0.0f0
        @test result < 1.0f-7  # Should be a small positive value
    end

    @testset "Out of domain" begin
        itp = EarthSciData.DataSetInterpolator{Float32}(fs, "NOX", ts, te, domain)
        @test interp!(itp, sample_time, deg2rad(0.0f0), deg2rad(40.0f0)) ≈ 0.0
    end
end

@testitem "polygons" setup=[NEISetup] begin
    itp = EarthSciData.DataSetInterpolator{Float32}(fileset, "NOX", ts, te, domain)
    polys = EarthSciData.get_geometry(fileset, itp.metadata)
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    for poly in polys
        for (x, y) in poly
            global xmin = min(xmin, x)
            global xmax = max(xmax, x)
            global ymin = min(ymin, y)
            global ymax = max(ymax, y)
        end
    end
    @test xmin ≈ -2.556e6
    @test ymin ≈ -1.728e6
    # With nx+1 edges, xmax and ymax include the full last cell
    @test xmax ≈ 2.952e6
    @test ymax ≈ 1.86e6
end

@testitem "monthly frequency" setup=[NEISetup] begin
    using Dates: month
    ts, te = DateTime(2016, 5, 1), DateTime(2016, 6, 1)
    fileset = EarthSciData.NEI2016MonthlyEmisFileSet("mrggrid_withbeis_withrwc", ts, te)
    sample_time = DateTime(2016, 5, 1)
    itp = EarthSciData.DataSetInterpolator{Float32}(fileset, "NOX", ts, te, domain)
    EarthSciData.lazyload!(itp, sample_time)
    ti = EarthSciData.DataFrequencyInfo(itp.fs.fs)
    @test month(itp.cache.times[1]) == 4
    @test month(itp.cache.times[2]) == 5

    sample_time = DateTime(2016, 5, 31)
    EarthSciData.lazyload!(itp, sample_time)
    @test month(itp.cache.times[1]) == 5
    @test month(itp.cache.times[2]) == 6
end

@testitem "run" setup=[NEISetup] begin
    using ModelingToolkit
    using OrdinaryDiffEqTsit5
    sys = mtkcompile(emis)
    prob = ODEProblem(
        sys,
        [lat => deg2rad(40.0), lon => deg2rad(-97.5), lev => 1.0],
        (0.0, 60.0),
    )
    solve(prob, Tsit5())
end

@testitem "run_nei" setup=[NEISetup] begin
    using ModelingToolkit, OrdinaryDiffEqTsit5
    using ModelingToolkit: t, D
    using DynamicQuantities
    domain = DomainInfo(
    DateTime(2016, 5, 16),
    DateTime(2016, 5, 17);
    lonrange=deg2rad(-125):deg2rad(0.625):deg2rad(-66.875),
    latrange=deg2rad(25):deg2rad(0.5):deg2rad(49),
    levrange=1:2,
    u_proto=zeros(Float64, 1, 1, 1, 1))

    emis_nei = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    @constants uc = 1.0 [unit = u"s" description = "unit conversion"]
    @variables ACET(t) = 0.0 [unit = u"1/s"]
    eq = D(ACET) ~ emis_nei.ACET / uc
    sys = compose(System([eq], t, [ACET], [uc]; name = :test_sys), emis_nei)
    sys = mtkcompile(sys)
    # After composition, parameters are namespaced; extract them from the compiled system
    ps = parameters(sys)
    lat_p = only(filter(p -> endswith(string(Symbol(p)), "₊lat"), ps))
    lon_p = only(filter(p -> endswith(string(Symbol(p)), "₊lon"), ps))
    lev_p = only(filter(p -> endswith(string(Symbol(p)), "₊lev"), ps))
    prob = ODEProblem(
        sys,
        [lat_p => deg2rad(40.0), lon_p => deg2rad(-97.5), lev_p => 1.0],
        (0.0, 60.0),
    )
    sol = solve(prob, Tsit5())
    @test sol.u[end][end] != 0.0  # Ensure we get a nonzero result
end

@testitem "diurnal_itp function" setup=[NEISetup] begin
    using Dates

    # Test domain starts at 2016-05-01 00:00:00 UTC
    t_ref_numeric = datetime2unix(DateTime(2016, 5, 1))  # Domain starts at 2016-05-01 00:00:00

    # Test 1: UTC (0° longitude)
    lon_utc = deg2rad(0.0)  # UTC timezone

    # 6 AM UTC
    six_am_utc = t_ref_numeric + 6 * 3600.0  # 7th factor
    @test EarthSciData.diurnal_itp(six_am_utc, lon_utc) == EarthSciData.DIURNAL_FACTORS[7]

    # 6 PM UTC
    six_pm_utc = t_ref_numeric + 18 * 3600.0  # 19th factor
    @test EarthSciData.diurnal_itp(six_pm_utc, lon_utc) == EarthSciData.DIURNAL_FACTORS[19]

    # Test 2: Chicago (UTC-6, longitude ~ -87.6°)
    lon_chicago = deg2rad(-87.6)  # Chicago longitude

    # 6 AM Chicago = 12 PM UTC (6 hours later)
    six_am_chicago = t_ref_numeric + 12 * 3600.0  # 6 AM Chicago = 12 PM UTC, 7th factor
    @test EarthSciData.diurnal_itp(six_am_chicago, lon_chicago) == EarthSciData.DIURNAL_FACTORS[7]

    # 6 PM Chicago = 12 AM UTC next day (6 hours later)
    six_pm_chicago = t_ref_numeric + 24 * 3600.0  # 6 PM Chicago = 12 AM UTC next day, 19th factor
    @test EarthSciData.diurnal_itp(six_pm_chicago, lon_chicago) == EarthSciData.DIURNAL_FACTORS[19]

    # Test that the function wraps around 24 hours correctly
    # 24 hours = 86400 seconds
    next_midnight_utc = t_ref_numeric + 24 * 3600.0  # 24 hours since start
    @test EarthSciData.diurnal_itp(next_midnight_utc, lon_utc) == EarthSciData.DIURNAL_FACTORS[1]

    # Test fractional hours
    half_past_one_utc = t_ref_numeric + 1.5 * 3600.0  # 1.5 hours since start
    @test EarthSciData.diurnal_itp(half_past_one_utc, lon_utc) == EarthSciData.DIURNAL_FACTORS[2]  # Should floor to hour 1
end

@testitem "allocations" setup=[NEISetup] begin
     using AllocCheck
    if !Sys.iswindows() # Allocation tests don't seem to work on windows.
        AllocCheck.@check_allocs checkf(
            itp, t, loc1, loc2) = EarthSciData.interp_unsafe(itp, t, loc1, loc2)

        sample_time = DateTime(2016, 5, 1)
        itp = EarthSciData.DataSetInterpolator{Float32}(fileset, "NOX", ts, te, domain)
        interp!(itp, sample_time, deg2rad(-97.0f0), deg2rad(40.0f0))
        # If there is an error, it should occur in the proj library.
        # https://github.com/JuliaGeo/Proj.jl/issues/104
        try
            checkf(itp, sample_time, deg2rad(-97.0f0), deg2rad(40.0f0))
        catch err
            @warn "Allocation errors:\n$(err.errors)"
            @test length(err.errors) <= 3
            @test all([contains(string(s), "jl_get_pgcstack_static") for s in err.errors])
        end

        itp2 = EarthSciData.DataSetInterpolator{Float64}(fileset, "NOX", ts, te, domain)
        interp!(itp2, sample_time, deg2rad(-97.0), deg2rad(40.0))
        try # If there is an error, it should occur in the proj library.
            checkf(itp2, sample_time, deg2rad(-97.0), deg2rad(40.0))
        catch err
            @warn "Allocation errors:\n$(err.errors)"
            @test length(err.errors) <= 3
            @test all([contains(string(s), "jl_get_pgcstack_static") for s in err.errors])
        end
    end
end

@testitem "Coupling with GEOS-FP" setup=[NEISetup] begin
    using ModelingToolkit
    gfp = GEOSFP("4x5", domain)

    csys = couple(emis, gfp)
    sys = convert(System, csys)
    eqs = observed(sys)

    @test occursin("NEI2016MonthlyEmis₊lat(t) ~ GEOSFP₊lat", string(eqs))
end

@testitem "wrong year" setup=[NEISetup] begin
    using Dates
    sample_time = DateTime(2016, 5, 1)
    itp = EarthSciData.DataSetInterpolator{Float32}(fileset, "NOX", ts, te, domain)
    sample_time = DateTime(2017, 5, 1)
    @test_throws ArgumentError EarthSciData.lazyload!(itp, sample_time)
end

@testitem "delp_dry_surface_itp" setup=[NEISetup] begin
    @test EarthSciData.delp_dry_surface_itp(deg2rad(-94.375), deg2rad(44.5)) ≈ 14.721285536474896
    @test EarthSciData.delp_dry_surface_itp(deg2rad(-88.125), deg2rad(42.0)) ≈ 14.8498301901334
end

@testitem "conservative regridding" setup=[NEISetup] begin
    using ModelingToolkit
    # Test that conservative regridding works through the unified NEI2016MonthlyEmis function
    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 2);
        lonrange=deg2rad(-125):deg2rad(0.625):deg2rad(-66.875),
        latrange=deg2rad(25):deg2rad(0.5):deg2rad(49),
        levrange = 1:10
    )
    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    eqs = equations(emis)
    @test length(eqs) == 69

    @testset "DataSetInterpolator with conservative regridding" begin
        domain = DomainInfo(
            DateTime(2016, 5, 16, 12, 0, 0),
            DateTime(2016, 5, 17, 12, 0, 0);
            lonrange=deg2rad(-125):deg2rad(0.625):deg2rad(-66.875),
            latrange=deg2rad(25):deg2rad(0.5):deg2rad(49),
            levrange = 1:10
        )
        ts, te = get_tspan_datetime(domain)
        fileset = EarthSciData.NEI2016MonthlyEmisFileSet("mrggrid_withbeis_withrwc", ts, te)
        # Use the unified DataSetInterpolator with conservative regridding via regridder()
        fs = EarthSciData.FileSetWithRegridder(fileset, EarthSciData.regridder(fileset,
            EarthSciData.loadmetadata(fileset, "NO"), domain))
        itp = EarthSciData.DataSetInterpolator{Float64}(fs, "NO", ts, te, domain)
        @test itp.varname == "NO"
        @test itp.metadata !== nothing

        # Test that interpolation works
        result = interp!(itp, ts, deg2rad(-88.125), deg2rad(42.0))
        @test result > 0.0  # Should be nonzero for this location
    end

    @testset "Conservative regridding - coarse domain" begin
        domain = DomainInfo(
            DateTime(2016, 5, 1),
            DateTime(2016, 5, 2);
            lonrange=deg2rad(-125):deg2rad(2.5):deg2rad(-66.875),
            latrange=deg2rad(25):deg2rad(2.0):deg2rad(49),
            levrange = 1:10
        )
        ts, te = get_tspan_datetime(domain)
        fileset = EarthSciData.NEI2016MonthlyEmisFileSet("mrggrid_withbeis_withrwc", ts, te)
        fs = EarthSciData.FileSetWithRegridder(fileset, EarthSciData.regridder(fileset,
            EarthSciData.loadmetadata(fileset, "NO"), domain))
        itp = EarthSciData.DataSetInterpolator{Float64}(fs, "NO", ts, te, domain)

        # Test that interpolation works at a specific location
        result = interp!(itp, ts, deg2rad(-87.5), deg2rad(41.0))
        @test result > 0.0
    end
end

@testitem "emission values" setup=[NEISetup] begin
    using ModelingToolkit, DynamicQuantities
    using ModelingToolkit: t, D
    using OrdinaryDiffEqTsit5
    using Dates
    domain = DomainInfo(
        DateTime(2016, 5, 15),
        DateTime(2016, 5, 16);
        lonrange=deg2rad(-88.125):deg2rad(0.625):deg2rad(-86.875),
        latrange=deg2rad(42):deg2rad(0.5):deg2rad(43),
        levrange = 1:2
    )
    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    eqs = equations(emis)
    @test length(eqs) == 69

    ts, te = get_tspan_datetime(domain)

    # Setup output array
    t_end = 1800  # in seconds
    nt = 5
    lon_grid = deg2rad(-88.125):deg2rad(0.625):deg2rad(-86.875)
    lat_grid = deg2rad(42):deg2rad(0.5):deg2rad(43)
    NO_map = Array{Float64}(undef, length(lon_grid), length(lat_grid), nt)
    tspan = (0.0, t_end)

    @constants uc = 1.0 [unit = u"s", description = "unit conversion"]
    @variables NO(t) = 0.0 [unit = u"1/s"]

    saveat = range(tspan[1], tspan[2], length=nt)

    for (i, lon_val) in enumerate(lon_grid)
        for (j, lat_val) in enumerate(lat_grid)
            eq = D(NO) ~ emis.NO / uc
            sys = compose(System([eq], t, [NO], [uc]; name = Symbol("NO_sys_$(i)_$(j)")), emis)
            sys = mtkcompile(sys)
            # After composition, parameters are namespaced; extract them from the compiled system
            ps = parameters(sys)
            lat_p = only(filter(p -> endswith(string(Symbol(p)), "₊lat"), ps))
            lon_p = only(filter(p -> endswith(string(Symbol(p)), "₊lon"), ps))
            lev_p = only(filter(p -> endswith(string(Symbol(p)), "₊lev"), ps))

            prob = ODEProblem(sys,
                [lat_p => lat_val, lon_p => lon_val, lev_p => 1.0],
                tspan)

            sol = solve(prob, Tsit5(), saveat=saveat)
            NO_map[i, j, :] = getindex.(sol.u, 1)
        end
    end

    @test NO_map[1, 1, end] != 0.0  # Ensure we get nonzero emissions
end

@testitem "emission values -- ACET" setup=[NEISetup] begin
    using ModelingToolkit, DynamicQuantities
    using ModelingToolkit: t, D
    using OrdinaryDiffEqTsit5
    using Dates
    domain = DomainInfo(
        DateTime(2016, 5, 15),
        DateTime(2016, 5, 16);
        lonrange=deg2rad(-88.125):deg2rad(0.625):deg2rad(-86.875),
        latrange=deg2rad(42):deg2rad(0.5):deg2rad(43),
        levrange = 1:2
    )
    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)

    ts, te = get_tspan_datetime(domain)

    t_end = 3600  # in seconds
    nt = 5
    lon_grid = deg2rad(-88.125):deg2rad(0.625):deg2rad(-86.875)
    lat_grid = deg2rad(42):deg2rad(0.5):deg2rad(43)
    ACET_map = Array{Float64}(undef, length(lon_grid), length(lat_grid), nt)
    tspan = (0.0, t_end)

    @constants uc = 1.0 [unit = u"s", description = "unit conversion"]
    @variables ACET(t) = 0.0 [unit = u"1/s"]

    saveat = range(tspan[1], tspan[2], length=nt)

    for (i, lon_val) in enumerate(lon_grid)
        for (j, lat_val) in enumerate(lat_grid)
            eq = D(ACET) ~ emis.ACET / uc
            sys = compose(System([eq], t, [ACET], [uc]; name = Symbol("ACET_sys_$(i)_$(j)")), emis)
            sys = mtkcompile(sys)
            # After composition, parameters are namespaced; extract them from the compiled system
            ps = parameters(sys)
            lat_p = only(filter(p -> endswith(string(Symbol(p)), "₊lat"), ps))
            lon_p = only(filter(p -> endswith(string(Symbol(p)), "₊lon"), ps))
            lev_p = only(filter(p -> endswith(string(Symbol(p)), "₊lev"), ps))

            prob = ODEProblem(sys,
                [lat_p => lat_val, lon_p => lon_val, lev_p => 1.0],
                tspan)

            sol = solve(prob, Tsit5(), saveat=saveat)
            ACET_map[i, j, :] = getindex.(sol.u, 1)
        end
    end

    @test ACET_map[1, 1, end] != 0.0  # Ensure we get nonzero emissions
end

@testitem "local-time clock conventions" begin
    using EarthSciData
    using Dates: DateTime, datetime2unix

    # Dallas (96.8°W), 12:00 UTC: NOx clock round(lon/15) = -6 h, CO clock floor(lon/15) = -7 h.
    t = datetime2unix(DateTime(2016, 7, 15, 12))
    dallas = deg2rad(-96.8)
    @test EarthSciData.diurnal_itp_NOx(t, dallas) == EarthSciData.DIURNAL_FACTORS_NOx[6 + 1]
    @test EarthSciData.diurnal_itp(t, dallas) == EarthSciData.DIURNAL_FACTORS[5 + 1]

    # Half-hour meridians on the 0.625° grid; EDGAR file offsets are -6, -6, -8 (round half to even).
    for (lon_deg, offset) in ((-97.5, -6), (-82.5, -6), (-112.5, -8))
        @test EarthSciData.diurnal_itp_NOx(t, deg2rad(lon_deg)) ==
              EarthSciData.DIURNAL_FACTORS_NOx[12 + offset + 1]
    end

    # 2016-07-17 06:30 UTC is Sunday; at 97.5°W it is Sunday on the NOx clock, Saturday on the CO clock.
    ts = datetime2unix(DateTime(2016, 7, 17, 6, 30))
    @test EarthSciData.dayofweek_itp_NOx(ts, deg2rad(-97.5)) == EarthSciData.DayofWeekFactors_NOx[7]
    @test EarthSciData.dayofweek_itp_CO(ts, deg2rad(-97.5)) == EarthSciData.DayofWeekFactors_CO[6]

    @test EarthSciData.DayofWeekFactors_NOx == [1.0706, 1.0706, 1.0706, 1.0706, 1.0706, 0.863, 0.784]
    @test EarthSciData.DayofWeekFactors_CO == [1.1076, 1.1076, 1.1076, 1.1076, 1.1076, 0.779, 0.683]
    @test sum(EarthSciData.DayofWeekFactors_NOx) ≈ 7
    @test sum(EarthSciData.DayofWeekFactors_CO) ≈ 7
end

@testitem "layer fractions" begin
    using EarthSciData
    names = EarthSciData.NEI2016_PROFILE_NAMES
    for (i, p) in enumerate(names)
        @test sum(EarthSciData.NEI2016_LAYER_FRACTIONS[i]) ≈ 1
        @test sum(EarthSciData.nei_layer_fraction(i, k) for k in 1:15) ≈ 1
    end
    f(p, k) = EarthSciData.nei_layer_fraction(findfirst(==(p), names), k)
    @test f("surface", 1) == 1
    @test f("surface", 2) == 0
    # Layer fractions of GEOS-Chem's NEI2016 3-D files (2016fh_16j_<sector>_0pt1degree_3D_month_05).
    @test f("ptegu", 1) ≈ 0.0279 atol = 1e-4
    @test sum(f("ptegu", k) for k in 4:11) ≈ 0.651 atol = 1e-3
    @test f("ptnonipm", 1) ≈ 0.7537 atol = 1e-4
    @test f("pt_oilgas", 1) ≈ 0.9188 atol = 1e-4
    @test f("cmv", 1) ≈ 0.5375 atol = 1e-4
    @test f("cmv", 3) > 0
    @test f("cmv", 4) == 0
    @test f("ptegu", 12) == 0
    @test EarthSciData.NEI2016_SECTOR_PROFILES["emln_ptegu"] == "ptegu"
    @test EarthSciData.NEI2016_SECTOR_PROFILES["emln_cmv_c3_12"] == "cmv"
    @test !haskey(EarthSciData.NEI2016_SECTOR_PROFILES, "emln_othpt")
end

@testitem "delp_dry_itp" begin
    using EarthSciData
    lon, lat = deg2rad(-94.375), deg2rad(44.5)
    d1 = EarthSciData.delp_dry_surface_itp(lon, lat)
    @test EarthSciData.delp_dry_itp(lon, lat, 1) == d1
    # The lowest GEOS-FP layers are all about 15 hPa thick.
    for k in 2:11
        @test 0.9 < EarthSciData.delp_dry_itp(lon, lat, k) / d1 < 1.1
    end
    # Consistent with the hybrid grid: thickness = ΔAp + ΔBp * ps (Ap is in Pa, delp in hPa).
    ΔAp(k) = (EarthSciData.Ap(k) - EarthSciData.Ap(k + 1)) / 100
    ΔBp(k) = EarthSciData.Bp(k) - EarthSciData.Bp(k + 1)
    ps = (d1 - ΔAp(1)) / ΔBp(1)
    @test EarthSciData.delp_dry_itp(lon, lat, 5) ≈ ΔAp(5) + ΔBp(5) * ps
end

@testitem "multi-sector sum" setup=[NEISetup] begin
    # Listing the same sector twice must give exactly twice the emissions.
    sector = "mrggrid_withbeis_withrwc"
    fs2 = EarthSciData.NEI2016MonthlyEmisMultiFileSet([sector, sector], ts, te)
    @test EarthSciData.varnames(fs2) == EarthSciData.varnames(fileset)
    @test EarthSciData.verify_fileset_interface(EarthSciData.NEI2016MonthlyEmisMultiFileSet)
    rg = EarthSciData.regridder(fs2, EarthSciData.loadmetadata(fs2, "NOX"), domain)
    itp1 = EarthSciData.DataSetInterpolator{Float64}(
        EarthSciData.FileSetWithRegridder(fileset, rg), "NOX", ts, te, domain)
    itp2 = EarthSciData.DataSetInterpolator{Float64}(
        EarthSciData.FileSetWithRegridder(fs2, rg), "NOX", ts, te, domain)
    v1 = interp!(itp1, sample_time, deg2rad(-97.0), deg2rad(40.0))
    v2 = interp!(itp2, sample_time, deg2rad(-97.0), deg2rad(40.0))
    @test v1 > 0
    @test v2 ≈ 2v1
    @test_throws ErrorException EarthSciData.NEI2016MonthlyEmisMultiFileSet(String[], ts, te)
    @test_throws ErrorException NEI2016MonthlyEmis([sector], domain;
        vertical_profiles = Dict(sector => "stack"))
end

@testsnippet NEIRate begin
    using ModelingToolkit, DynamicQuantities, OrdinaryDiffEqTsit5
    using ModelingToolkit: t, D
    using Dates: DateTime
    using EarthSciMLBase, EarthSciData

    # NO emission rate (1/s) of an emission system at a location and model layer.
    function no_rate(emis, lon_val, lat_val, lev_val)
        @constants uc = 1.0 [unit = u"s"]
        @variables NO(t) = 0.0 [unit = u"1/s"]
        sys = compose(System([D(NO) ~ emis.NO / uc], t, [NO], [uc]; name = :no_sys), emis)
        sys = mtkcompile(sys)
        ps = parameters(sys)
        getp(suffix) = only(filter(p -> endswith(string(Symbol(p)), suffix), ps))
        prob = ODEProblem(sys,
            [getp("₊lat") => lat_val, getp("₊lon") => lon_val, getp("₊lev") => lev_val],
            (0.0, 60.0))
        sol = solve(prob, Tsit5())
        sol.u[end][1] / 60.0
    end
end

@testitem "vertical allocation" setup=[NEIRate] begin
    domain = DomainInfo(
        DateTime(2016, 5, 15),
        DateTime(2016, 5, 16);
        lonrange = deg2rad(-88.125):deg2rad(0.625):deg2rad(-86.875),
        latrange = deg2rad(42):deg2rad(0.5):deg2rad(43),
        levrange = 1:11
    )
    lon_val, lat_val = deg2rad(-88.125), deg2rad(42.0)
    sector = "mrggrid_withbeis_withrwc"

    surface = NEI2016MonthlyEmis(sector, domain)
    # The same file treated as an elevated sector: the column total is spread over layers.
    elevated = NEI2016MonthlyEmis([sector], domain;
        vertical_profiles = Dict(sector => "ptegu"), name = :elevated)
    @test length(equations(elevated)) == 69

    r1 = no_rate(surface, lon_val, lat_val, 1.0)
    @test r1 > 0
    @test no_rate(surface, lon_val, lat_val, 2.0) == 0

    ptegu = findfirst(==("ptegu"), EarthSciData.NEI2016_PROFILE_NAMES)
    f(k) = EarthSciData.nei_layer_fraction(ptegu, k)
    d(k) = EarthSciData.delp_dry_itp(lon_val, lat_val, k)
    @test no_rate(elevated, lon_val, lat_val, 1.0) ≈ r1 * f(1) rtol = 1e-6
    @test no_rate(elevated, lon_val, lat_val, 4.0) ≈ r1 * f(4) * d(1) / d(4) rtol = 1e-6
    @test no_rate(elevated, lon_val, lat_val, 12.0) == 0
end

@testitem "surface regression" begin
    # The default single-sector call must reproduce the surface emissions of the code before
    # multi-sector support and vertical allocation were added (commit af45420): emission
    # values (1/s) at 2016-05-15 00:00 UTC in layer 1, evaluated directly from the emission
    # system with that version. Evaluating at one instant avoids the solver noise of an ODE
    # integration.
    using ModelingToolkit, Dates, EarthSciMLBase, EarthSciData
    using SymbolicIndexingInterface: getsym, setp
    domain = DomainInfo(
        DateTime(2016, 5, 15),
        DateTime(2016, 5, 16);
        lonrange = deg2rad(-88.125):deg2rad(0.625):deg2rad(-86.875),
        latrange = deg2rad(42):deg2rad(0.5):deg2rad(43),
        levrange = 1:2
    )
    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    @test length(equations(emis)) == 69
    @test !any(occursin("₊", string(p)) for p in parameters(emis)) # no namespaced groups

    sys = mtkcompile(emis)
    prob = ODEProblem(sys,
        [sys.lon => deg2rad(-88.125), sys.lat => deg2rad(42.0), sys.lev => 1.0], (0.0, 60.0))
    set_lon, set_lat, set_lev = setp(prob, sys.lon), setp(prob, sys.lat), setp(prob, sys.lev)
    function value(sp, lo, la, lev)
        set_lon(prob, deg2rad(lo)); set_lat(prob, deg2rad(la)); set_lev(prob, lev)
        getsym(prob, getproperty(sys, Symbol(sp)))(prob)
    end

    expected = Dict(
        ("NO", -88.125, 42.0) => 3.696630216752866e-12,
        ("NO2", -88.125, 42.0) => 4.3730420959836544e-13,
        ("CO", -88.125, 42.0) => 9.478783058412906e-12,
        ("FORM", -88.125, 42.0) => 4.037229705850393e-14,
        ("ISOP", -88.125, 42.0) => 5.184836779214539e-13,
        ("ACET", -88.125, 42.0) => 1.3005873181398365e-13,
        ("NO", -87.5, 42.5) => 1.3329110156415248e-13,
        ("NO2", -87.5, 42.5) => 1.42753478357713e-14,
        ("CO", -87.5, 42.5) => 3.628395844549008e-13,
        ("FORM", -87.5, 42.5) => 1.7564883575452577e-15,
        ("ISOP", -87.5, 42.5) => 9.160316615712789e-15,
        ("ACET", -87.5, 42.5) => 7.801309939082226e-15,
        ("NO", -86.875, 43.0) => 2.714099586734287e-15,
        ("NO2", -86.875, 43.0) => 2.774412889837166e-16,
        ("CO", -86.875, 43.0) => 1.2543254094733087e-14,
        ("FORM", -86.875, 43.0) => 2.803387435072575e-17,
        ("ISOP", -86.875, 43.0) => 1.1051762504955001e-17,
        ("ACET", -86.875, 43.0) => 3.841633770766277e-18
    )
    # Differences are judged against each species' largest value: at grid cells the loader
    # now reads the regridded field directly, while the old code's interpolator carried
    # rounding noise of about 1e-15 of the neighbouring values.
    speciesmax = Dict(sp => maximum(v for ((s, _, _), v) in expected if s == sp)
                      for ((sp, _, _), _) in expected)
    for ((sp, lo, la), v) in expected
        @test isapprox(value(sp, lo, la, 1.0), v; rtol = 1e-12, atol = 1e-12 * speciesmax[sp])
    end
    @test value("NO", -88.125, 42.0, 2.0) == 0.0
end

@testitem "direct grid reads" setup=[NEIRate] begin
    # At model grid cells the emission callable reads the loaded monthly slices directly;
    # the values must equal what the interpolator returns there. Off the grid it falls
    # back to the interpolator.
    using Dates: datetime2unix
    domain = DomainInfo(
        DateTime(2016, 5, 15),
        DateTime(2016, 5, 16);
        lonrange = deg2rad(-90.625):deg2rad(0.625):deg2rad(-89.375),
        latrange = deg2rad(28):deg2rad(0.5):deg2rad(29),
        levrange = 1:4
    )
    emis = NEI2016MonthlyEmis(["mrggrid_withbeis_withrwc", "emln_cmv_c3_12"], domain)
    no_p = only(filter(p -> string(Symbol(p)) == "NO_itp", parameters(emis)))
    w = ModelingToolkit.getdefault(no_p)
    tref = datetime2unix(DateTime(2016, 5, 15))
    xs, ys = EarthSciMLBase.grid(domain, (false, false, false))[1:2]
    @test w.grid[3] == length(xs) && w.grid[6] == length(ys)
    @test size(w.delp) == (length(xs), length(ys), 3)
    w(tref, xs[1], ys[1], 1.0) # load data
    for g in eachindex(w.itps)
        fieldmax = 0.0
        maxdiff = 0.0
        for (i, x) in enumerate(xs), (j, y) in enumerate(ys), dt in (0.0, 1234.567, 86399.9)
            a = EarthSciData._nei_grid_value(w.itps[g], tref + dt, i, j)
            b = EarthSciData.interp_unsafe(w.itps[g], tref + dt, x, y)
            fieldmax = max(fieldmax, abs(b))
            maxdiff = max(maxdiff, abs(a - b))
        end
        @test fieldmax > 0
        @test maxdiff <= 1e-12 * fieldmax
    end
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys), k in 1:3
        @test w.delp[i, j, k] ≈ EarthSciData.delp_dry_itp(x, y, k)
    end
    # Off-grid points use the interpolator: a point very close to a cell gives nearly the same value.
    x, y = xs[2], ys[2]
    on = w(tref, x, y, 1.0)
    off = w(tref, x + 1e-9, y, 1.0)
    @test on > 0
    @test off ≈ on rtol = 1e-5
    @test w(tref, x, y, 5.0) == 0
end

@testitem "domain-top fold" setup=[NEIRate] begin
    # A power-plant profile reaches layer 11; on a 3-layer domain the fractions of layers
    # 3-11 go into layer 3 so the column total is kept, and nothing is emitted above.
    domain = DomainInfo(
        DateTime(2016, 5, 15),
        DateTime(2016, 5, 16);
        lonrange = deg2rad(-88.125):deg2rad(0.625):deg2rad(-86.875),
        latrange = deg2rad(42):deg2rad(0.5):deg2rad(43),
        levrange = 1:3
    )
    lon_val, lat_val = deg2rad(-88.125), deg2rad(42.0)
    sector = "mrggrid_withbeis_withrwc"
    surface = NEI2016MonthlyEmis(sector, domain)
    elevated = NEI2016MonthlyEmis([sector], domain;
        vertical_profiles = Dict(sector => "ptegu"), name = :elevated)

    ptegu = findfirst(==("ptegu"), EarthSciData.NEI2016_PROFILE_NAMES)
    f = EarthSciData.NEI2016_LAYER_FRACTIONS[ptegu]
    folded = EarthSciData.nei_domain_fractions(ptegu, 3)
    @test length(folded) == 3
    @test sum(folded) ≈ 1
    @test folded[3] ≈ sum(f[3:end])
    @test EarthSciData.nei_domain_fractions(ptegu, 30) == f

    d(k) = EarthSciData.delp_dry_itp(lon_val, lat_val, k)
    r1 = no_rate(surface, lon_val, lat_val, 1.0)
    @test no_rate(elevated, lon_val, lat_val, 1.0) ≈ r1 * f[1] rtol = 1e-6
    @test no_rate(elevated, lon_val, lat_val, 3.0) ≈ r1 * sum(f[3:end]) * d(1) / d(3) rtol = 1e-6
    @test no_rate(elevated, lon_val, lat_val, 4.0) == 0
    # Column total over the domain's layers equals the surface-only column.
    col = sum(no_rate(elevated, lon_val, lat_val, Float64(k)) * d(k) for k in 1:3)
    @test col ≈ r1 * d(1) rtol = 1e-6
end

@testitem "multi-sector groups" setup=[NEIRate] begin
    # Surface file plus ocean-going ships (Gulf of Mexico off Louisiana).
    domain = DomainInfo(
        DateTime(2016, 5, 15),
        DateTime(2016, 5, 16);
        lonrange = deg2rad(-90.625):deg2rad(0.625):deg2rad(-89.375),
        latrange = deg2rad(28):deg2rad(0.5):deg2rad(29),
        levrange = 1:3
    )
    lon_val, lat_val = deg2rad(-90.0), deg2rad(28.5)
    sectors = ["mrggrid_withbeis_withrwc", "emln_cmv_c3_12"]

    combined = NEI2016MonthlyEmis(sectors, domain; name = :combined)
    surface = NEI2016MonthlyEmis(sectors[1], domain; name = :surface_only)
    ships = NEI2016MonthlyEmis(sectors[2], domain; name = :ships_only)

    @test length(equations(combined)) == 69 # union of the variables of both files
    # One callable parameter per species, holding one interpolator per profile group.
    no_p = only(filter(p -> string(Symbol(p)) == "NO_itp", parameters(combined)))
    w = ModelingToolkit.getdefault(no_p)
    @test w isa EarthSciData.NEILayeredEmission
    @test length(w.itps) == 2
    @test w.ktop == 3

    s1 = no_rate(ships, lon_val, lat_val, 1.0)
    s2 = no_rate(ships, lon_val, lat_val, 2.0)
    @test s1 > 0
    cmv = findfirst(==("cmv"), EarthSciData.NEI2016_PROFILE_NAMES)
    f(k) = EarthSciData.nei_layer_fraction(cmv, k)
    d(k) = EarthSciData.delp_dry_itp(lon_val, lat_val, k)
    @test s2 ≈ s1 * f(2) / f(1) * d(1) / d(2) rtol = 1e-6
    @test no_rate(ships, lon_val, lat_val, 4.0) == 0

    @test no_rate(combined, lon_val, lat_val, 1.0) ≈
          no_rate(surface, lon_val, lat_val, 1.0) + s1 rtol = 1e-6
    @test no_rate(combined, lon_val, lat_val, 2.0) ≈ s2 rtol = 1e-6
end
