import dash
from dash.dependencies import Input, Output
from dash.long_callback import DiskcacheLongCallbackManager
import dash_bootstrap_components as dbc
from dash import dcc, html, dash_table
from datetime import date, datetime
import geopandas as gpd
import pandas as pd
import pprint
import plotly.express as px
import json
import diskcache
from sgp4.api import Satrec, WGS72
from skyfield.api import EarthSatellite, load, wgs84

import coverage as cov

# From: https://towardsdatascience.com/long-callbacks-in-dash-web-apps-72fd8de25937
cache = diskcache.Cache("./cache")
lcm = DiskcacheLongCallbackManager(cache)

# Standard Dash/ Flask setup
app = dash.Dash(
    __name__, long_callback_manager=lcm, external_stylesheets=[dbc.themes.BOOTSTRAP]
)

server = app.server

# Load ephemeris, timescale, and world
eph = load("de421.bsp")
ts = load.timescale()
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

# Load settings from coverage model (TODO: from file)
settings = cov.settings

# TODO: Add check for pre-existing .html?
# cov.gen_coverage_plot()

## Can use stations/ satcat text files for more sats to add
# stations_url = "http://celestrak.com/NORAD/elements/stations.txt"
# stations_url = "./satcat.txt"
# satellites = load.tle_file(stations_url)

satellites = cov.gen_sats(settings["sat_nos"])
df = pd.DataFrame(
    {
        "name": [sat[0].name for sat in satellites],
        "satnum": [sat[0].model.satnum for sat in satellites],
        "epoch": [sat[0].epoch.utc_jpl() for sat in satellites],
    }
)

app.layout = html.Div(
    [
        dcc.Store(id="settings", data=json.dumps(settings)),
        dbc.Row(
            [
                html.H1(children="LEO Satellite Orbits & Coverage"),
                html.Div(
                    children="""
            Plots the 2D (lat, lon) ground tracks for multiple low-Earth orbit satellites, using Celestrack TLE's and `skyfield` for propagation with SPG4.
            For best performance, limit date range to 1-3 days.
            """
                ),
                dcc.DatePickerRange(
                    id="date-range",
                    month_format="MMMM Y",
                    end_date_placeholder_text="MMMM Y",
                    start_date=date(2022, 4, 20),
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dash_table.DataTable(
                                id="datatable-row-ids",
                                columns=[
                                    {"name": i, "id": i, "deletable": True}
                                    for i in df.columns
                                    # omit the id column
                                    if i != "id"
                                ],
                                data=df.to_dict("records"),
                                editable=True,
                                filter_action="native",
                                sort_action="native",
                                sort_mode="multi",
                                row_selectable="multi",
                                row_deletable=True,
                                selected_rows=[0, 1],
                                page_action="native",
                                page_current=0,
                                page_size=10,
                            ),
                        ],
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        dcc.Graph(
                            id="orbit-graph",
                            # figure=blank_map (TODO)
                        ),
                    ),
                    width=6,
                ),
            ],
            justify="left",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div("Continent/ Country selection"),
                        dcc.Dropdown(
                            # world.continent.unique(),
                            # "North America",
                            world.name.unique(),
                            "Brazil",
                            id="continent-select",
                        ),
                        dcc.RadioItems(
                            ["Swath Viewer", "Revisit Viewer"],
                            "Swath Viewer",
                            id="map-type",
                        ),
                    ]
                ),
                dbc.Col(
                    [
                        html.Div(["This is my map"]),
                        html.Iframe(
                            id="map-frame",
                            srcDoc=open("./cache/map.html", "r").read(),
                            width="80%",
                            height="400",
                        ),
                    ]
                ),
            ]
        ),
    ]
)


@app.callback(
    Output("settings", "data"),
    Input("date-range", "start_date"),
    Input("date-range", "end_date"),
    Input("datatable-row-ids", "selected_rows"),
    prevent_initial_call=False,
)
def update_settings(startdate, enddate, satrows):
    # Update date range
    format = "%Y-%m-%d"
    start_dt = datetime.strptime(startdate, format)
    if enddate:
        end_dt = datetime.strptime(enddate, format)
        ndays = end_dt.day - start_dt.day
    else:
        ndays = 1
    settings["time_range"] = dict(
        start_yr=start_dt.year,
        start_mo=start_dt.month,
        start_day=start_dt.day,
        days=ndays,
        step_min=1,
    )

    # Update satellites
    settings["sat_ids"] = df.iloc[satrows].satnum.to_json()
    with open("settings.json", "w") as outfile:
        json.dump(settings, outfile, indent=4)

    return settings


# @app.long_callback(
@app.callback(
    Output("orbit-graph", "figure"),
    Input("settings", "data"),
    prevent_initial_call=False,
)
def update_orbit_plot(settings):

    return cov.gen_orbit_plot(settings)


@app.long_callback(
    Output("map-frame", "srcDoc"),
    Input("continent-select", "value"),
    Input("map-type", "value"),
    Input("settings", "data"),
    prevent_initial_call=True,
)
def update_map(continent, maptype, settings):
    m = cov.gen_coverage_plot(continent, maptype, settings)
    return open("./cache/map.html", "r").read()


if __name__ == "__main__":
    app.run_server(debug=True)
