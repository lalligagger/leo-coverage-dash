import numpy as np
from skyfield.api import load, wgs84, Distance
from skyfield.toposlib import ITRSPosition
from skyfield.framelib import itrs
from scipy.spatial.transform import Rotation

import folium
import branca
import shapely
import plotly.express as px
from shapely.geometry import Polygon
import geopandas as gpd
import pandas as pd
import json

# TODO: Create main function that generates either a FOV (Swath) map or a revisit (gridded raster) map.
# i.e., python coverage.py -t "swath" for swath type from cmd line, looks at settings.json for config.


def gen_sats(sat_nos=[39084, 49260]):
    """
    Skyfield satellite lookup from Celestrack, based on catalog ID.
    Landsat 8 & 9 defaults.
    """
    sats = []
    for n in sat_nos:
        url = "https://celestrak.com/satcat/tle.php?CATNR={}".format(n)
        tle_filename = "tle-CATNR-{}.txt".format(n)
        sat = load.tle_file(url, filename=tle_filename)
        sats.append(sat)

    return sats


def gen_times(start_yr=2021, start_mo=11, start_day=20, days=1, step_min=1):
    """
    Generate skyfield timespan over desired range.
    """
    ts = load.timescale()
    times = ts.utc(start_yr, start_mo, start_day, 0, range(0, 60 * 24 * days, step_min))

    return times


def gen_instrument(
    name="instrument", fl=178, pitch=0.025, h_pix=1850, v_pix=1800, mm=True
):
    """
    Takes in instrument parameters and calculates the azimuth offset to generate azimuth angles to top corners, and the half-diagonal FOV in angle space.
    For v2, we use the horizontal and vertical FOVs instead, which are divided by 2 and applied in the test notebook as (az, el) offsets in LVLH...
    The complete function needs to be integrated somewhere below.

    Defaults are TIRS.
    """
    hfov_deg = np.degrees(h_pix * pitch / fl)
    vfov_deg = np.degrees(v_pix * pitch / fl)
    instrument = {
        "name": name,
        "fl": fl,
        "pitch": pitch,
        "h_pix": h_pix,
        "v_pix": v_pix,
        "mm": mm,
        "hfov_deg": hfov_deg,  # update to atan
        "vfov_deg": vfov_deg,  # update to atan
        "half_diag_deg": np.degrees(
            (pitch / fl) * np.sqrt(h_pix**2 + v_pix**2) / 2
        ),
        "az1": np.degrees(np.arctan2(h_pix, v_pix)),
        "az2": 360 - np.degrees(np.arctan2(h_pix, v_pix)),
        "corners": {
            "c1": {"X": -hfov_deg / 2, "Y": vfov_deg / 2},
            "c2": {"X": hfov_deg / 2, "Y": vfov_deg / 2},
            "c3": {"X": hfov_deg / 2, "Y": -vfov_deg / 2},
            "c4": {"X": -hfov_deg / 2, "Y": -vfov_deg / 2},
        },
    }
    return instrument


def los_to_earth(position, pointing):
    """Find the intersection of a pointing vector with the Earth
    Finds the intersection of a pointing vector u and starting point s with the WGS-84 geoid
    Source: https://stephenhartzell.medium.com/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6

    Args:
        position (np.array): length 3 array defining the starting point location(s) in meters
        pointing (np.array): length 3 array defining the pointing vector(s) (must be a unit vector)
    Returns:
        np.array: length 3 defining the point(s) of intersection with the surface of the Earth in meters
    """

    a = 6378.137
    b = 6378.137
    c = 6356.752314245
    x = position[0]
    y = position[1]
    z = position[2]
    u = pointing[0]
    v = pointing[1]
    w = pointing[2]

    value = (
        -(a**2) * b**2 * w * z - a**2 * c**2 * v * y - b**2 * c**2 * u * x
    )
    radical = (
        a**2 * b**2 * w**2
        + a**2 * c**2 * v**2
        - a**2 * v**2 * z**2
        + 2 * a**2 * v * w * y * z
        - a**2 * w**2 * y**2
        + b**2 * c**2 * u**2
        - b**2 * u**2 * z**2
        + 2 * b**2 * u * w * x * z
        - b**2 * w**2 * x**2
        - c**2 * u**2 * y**2
        + 2 * c**2 * u * v * x * y
        - c**2 * v**2 * x**2
    )
    magnitude = (
        a**2 * b**2 * w**2 + a**2 * c**2 * v**2 + b**2 * c**2 * u**2
    )

    if radical < 0:
        raise ValueError("The Line-of-Sight vector does not point toward the Earth")
    d = (value - a * b * c * np.sqrt(radical)) / magnitude

    if d < 0:
        raise ValueError("The Line-of-Sight vector does not point toward the Earth")

    return np.array(
        [
            x + d * u,
            y + d * v,
            z + d * w,
        ]
    )


def get_los(sat, time):
    geo = sat.at(time)
    xyz_dist_rates = geo.frame_xyz_and_velocity(itrs)
    xyz_dist = xyz_dist_rates[0]
    pointing = -xyz_dist.km / xyz_dist.length().km

    xyz_vel = xyz_dist_rates[1]
    # bearing = xyz_vel.km_per_s / np.linalg.norm(xyz_vel.km_per_s)

    los_xyz = los_to_earth(xyz_dist.km, pointing)  # input is meters

    los = Distance(km=los_xyz)
    los_itrs = ITRSPosition(los)
    los_itrs.at(time).frame_xyz(itrs).km

    los_lat, los_lon = wgs84.latlon_of(los_itrs.at(time))
    d = np.sqrt(np.sum(np.square(xyz_dist.km - los_xyz)))

    return los_lat, los_lon, d


def get_lvlh_pointing(sat, time):
# def get_lvlh_pointing(sat_v_time):
    geo = sat.at(time)
    xyz_dist_rates = geo.frame_xyz_and_velocity(itrs)
    xyz_dist = xyz_dist_rates[0]
    pointing = -xyz_dist.km / xyz_dist.length().km

    xyz_vel = xyz_dist_rates[1]
    bearing = xyz_vel.km_per_s / np.linalg.norm(xyz_vel.km_per_s)

    neg_orb_normal = -np.cross(pointing, bearing)
    x_axis = np.cross(neg_orb_normal, pointing)
    lvlh = {"X": x_axis, "Y": neg_orb_normal, "Z": pointing}

    return lvlh, pointing


def get_inst_fov(sat, time, inst):
# def get_inst_fov(sat_v_time, inst):
    lvlh, pointing = get_lvlh_pointing(sat, time)
    # lvlh, pointing = get_lvlh_pointing(sat_v_time)
    xyz_dist_rates = sat.at(time).frame_xyz_and_velocity(itrs)
    # sat_v_time.frame_xyz...
    xyz_dist = xyz_dist_rates[0]
    z_rate = xyz_dist_rates[1]

    # Empty dict to populate with lat/ lons, mapped from cs_dict in angle space
    cs_lla_dict = {
        "c1": {"lat": None, "lon": None},
        "c2": {"lat": None, "lon": None},
        "c3": {"lat": None, "lon": None},
        "c4": {"lat": None, "lon": None},
    }

    # For each corner in FOV...
    for c in inst["corners"]:
        # Generate X and Y rotation vectors
        rot_X_deg = inst["corners"][c]["X"]
        rot_X_rad = np.radians(rot_X_deg)
        rot_X_ax = lvlh["X"]

        rot_Y_deg = inst["corners"][c]["Y"]
        rot_Y_rad = np.radians(rot_Y_deg)
        rot_Y_ax = lvlh["Y"]

        # Rotations with scipy:
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html
        # Rotate about X
        rot_X_vec = rot_X_rad * rot_X_ax
        rot_X = Rotation.from_rotvec(rot_X_vec)

        # Rotate about Y for final LOS
        rot_Y_vec = rot_Y_rad * rot_Y_ax
        rot_Y = Rotation.from_rotvec(rot_Y_vec)
        rot = rot_X * rot_Y
        los_XY = rot.apply(pointing)

        # Get Earth intercept of LOS, create ITRS position object
        los_xyz = los_to_earth(xyz_dist.km, los_XY)
        los_itrs = ITRSPosition(Distance(km=los_xyz))

        # Calculate intercept lat/ lon from ITRS frame
        los_lat, los_lon = wgs84.latlon_of(los_itrs.at(time))
        cs_lla_dict[c]["lat"] = los_lat.degrees
        cs_lla_dict[c]["lon"] = los_lon.degrees

    return cs_lla_dict


def forecast_fovs(sat, times, inst):
# def forecast_fovs(sat_v_time, inst):
    # Create temporary function that can be vectorized
    def gen_fov_poly(time):
        # gets satellite geoms at time steps
        # Get the ITRS position of the satellite as origin of LVLH frame.

        xyz_dist_rates = sat.at(time).frame_xyz_and_velocity(itrs)
        # xyz_dist_rates = sat_v_time.frame_xyz_and_velocity(itrs)

        # xyz_dist = xyz_dist_rates[0]
        z_rate = xyz_dist_rates[1]

        # TODO: Make this filter a selection (ascending/ descending)
        # Descending only filter:
        if z_rate.km_per_s[2] < 0:
            cs_lla_dict = get_inst_fov(sat, time, inst)

            # Add lat, lon offset for each corner of FOV
            return Polygon(
                [
                    (cs_lla_dict["c1"]["lon"], cs_lla_dict["c1"]["lat"]),
                    (cs_lla_dict["c2"]["lon"], cs_lla_dict["c2"]["lat"]),
                    (cs_lla_dict["c3"]["lon"], cs_lla_dict["c3"]["lat"]),
                    (cs_lla_dict["c4"]["lon"], cs_lla_dict["c4"]["lat"]),
                    (cs_lla_dict["c1"]["lon"], cs_lla_dict["c1"]["lat"]),
                ]
            )

    vfunc = np.vectorize(gen_fov_poly)
    polys = vfunc(times)

    poly_df = gpd.GeoDataFrame(data=polys, columns=["geometry"], crs="EPSG:4326")
    poly_df["satellite"] = sat.name
    poly_df["id"] = np.abs(sat.target)
    poly_df["time"] = times.utc_strftime()

    return poly_df


def create_grid(bounds, xcell_size, ycell_size):
    (xmin, ymin, xmax, ymax) = bounds

    # Create grid of points with regular spacing in degrees
    # projection of the grid
    crs = "EPSG:4326"

    xcells = np.arange(xmin, xmax + xcell_size, xcell_size)
    ycells = np.arange(ymin, ymax + ycell_size, ycell_size)
    grid_shape = (len(xcells), len(ycells))

    # create the grid points in a loop
    grid_points = []
    for x0 in xcells:
        for y0 in ycells:
            grid_points.append(shapely.geometry.Point(x0, y0))

    grid = gpd.GeoDataFrame(grid_points, columns=["geometry"], crs=crs)

    return grid, grid_shape


# TODO: Delete these default settings in favor of .json
settings = {
    "sat_nos": [39084, 49260],
    "time_range": dict(start_yr=2021, start_mo=12, start_day=10, days=1, step_min=1),
    "inst_params": dict(
        name="tirs", fl=178, pitch=0.025, h_pix=1850, v_pix=4000, mm=True
    ),
    "cell": 5,
    # "continent": "North America",
    "continent": "Brazil",
}


def gen_coverage_plot(continent, maptype, settings):
    tles = gen_sats(json.loads(settings["sat_ids"]).values())
    # tles = gen_sats(sat_nos=settings["sat_nos"]) # old way

    # Select cell size for coverage map
    xcell_size = ycell_size = settings["cell"]
    # Select AOI from gpd naturalearth dataset (filter by .name for country, .continent for continent)
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    # aoi = world[world.continent == continent].geometry
    aoi = world[world.name == continent].geometry

    ## TODO: Support country selection
    # aoi =  world[world.name == "Brazil"].geometry
    # aoi = world[world.name == "United States of America"]#.geometry

    # Batch FOV generation over N satellites

    gdfs = []
    colors = px.colors.qualitative.Plotly
    color_idx = 0

    for tle in tles:
        sat = tle[0]
        poly_df = forecast_fovs(
            sat,
            gen_times(**settings["time_range"]),
            # gen_instrument(settings["inst_params"]),
            gen_instrument(**settings["inst_params"]),
        )
        poly_df["color"] = colors[color_idx]
        gdfs.append(poly_df)
        color_idx += 1

    poly_df = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs="epsg:4326")
    poly_df["lonspan"] = poly_df.bounds["maxx"] - poly_df.bounds["minx"]

    # Filter shapes crossing anti-meridian
    plot_df = poly_df[poly_df["lonspan"] < 20].copy()
    print("original plot_df length:" + str(len(plot_df)))
    # plot_df.to_json("df_dump.json")
    # Filter by AOI
    xmin, ymin, xmax, ymax = aoi.total_bounds
    plot_df = plot_df.cx[xmin:xmax, ymin:ymax]
    print("filtered plot_df length:" + str(len(plot_df)))

    if maptype == "Swath Viewer":

        m = plot_df.explore(color="color", tooltip=["satellite", "time"])

        m.save("./cache/map.html")
        # return m

    if maptype == "Revisit Viewer":
        ## Coverage data analysis for single satellite/ batch of satellites
        # 1) Create a grid of equally spaced points
        grid, grid_shape = create_grid(aoi.total_bounds, xcell_size, ycell_size)

        # 2) Add "n_visits" column to grid using sjoin/ dissolve
        shapes = gpd.GeoDataFrame(plot_df.geometry)
        merged = gpd.sjoin(shapes, grid, how="left", predicate="intersects")

        # init merged, will be replaced with nan or positive int where n_visits > 0
        merged["n_visits"] = 0
        dissolve = merged.dissolve(
            by="index_right", aggfunc="count"
        )  # no difference in count vs. sum here?
        grid.loc[dissolve.index, "n_visits"] = dissolve.n_visits.values

        # 3) Form 2D array of n_visits based on grid shape
        img = np.rot90(grid.n_visits.values.reshape(grid_shape))

        # 4) Create colormap and apply to img
        colormap = branca.colormap.step.viridis.scale(1, grid.n_visits.max())

        def colorfunc(x):
            if np.isnan(x):
                return (0, 0, 0, 0)
            else:
                return colormap.rgba_bytes_tuple(x)

        # Apply cmap to img array and rearrange for RGBA
        cmap = np.vectorize(colorfunc)
        rgba_img = np.array(cmap(img))
        rgba_img = np.moveaxis(rgba_img, 0, 2)

        m = folium.Map()
        m.fit_bounds([[ymin, xmin], [ymax, xmax]])
        m.add_child(
            folium.raster_layers.ImageOverlay(
                rgba_img,
                opacity=0.4,
                mercator_project=True,  # crs="EPSG:4326",
                bounds=[[ymin, xmin], [ymax, xmax]],
            )
        )
        colormap.add_to(m)

        m.save("./cache/map.html")
        # return m


def gen_orbit_plot(settings):
    # TODO: Add configurable time step size (using "step_min")
    tr = settings["time_range"]
    ts = load.timescale()

    time_steps = ts.utc(
        tr["start_yr"],
        tr["start_mo"],
        tr["start_day"],
        0,
        # added 10 min sampling to orbit
        range(0, 60 * 24 * tr["days"], 10),
    )

    satellites = gen_sats(settings["sat_nos"])

    gdf_all = gpd.GeoDataFrame()
    for sat in satellites:

        # Gets satellite ground track at time steps
        geos = sat[0].at(time_steps)
        lat, lon = wgs84.latlon_of(geos)
        altitude = wgs84.height_of(geos).km

        df = pd.DataFrame({"lat": lat.degrees, "lon": lon.degrees})
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat))
        gdf["altitude_km"] = altitude
        gdf["sat_name"] = sat[0].name
        gdf.index = time_steps
        gdf_all = gpd.GeoDataFrame(pd.concat([gdf_all, gdf]))

    fig = px.line_geo(
        gdf_all,
        lat="lat",
        lon="lon",
        color="sat_name",
        projection="orthographic",
        # width=400
    )
    return fig
