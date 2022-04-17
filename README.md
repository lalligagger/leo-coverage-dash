# leo-coverage-dash
Dash app for low-Earth orbit (LEO) satellite coverage

Install repo2docker CLI tool https://repo2docker.readthedocs.io/en/latest/getting-started/index.html

Pull this repo, build container but don't run:

`repo2docker --image-name dash_app https://github.com/lalligagger/leo-coverage-dash.git`

Run app with port forwarding:

`docker run -p 8050:8050 dash_app python app.py`
