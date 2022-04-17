# leo-coverage-dash
Dash app for low-Earth orbit (LEO) satellite coverage

Install repo2docker CLI tool https://repo2docker.readthedocs.io/en/latest/getting-started/index.html

Pull this repo, open port forwarding for Dash, and run app:
`repo2docker -p 8050:8050 https://github.com/lalligagger/leo-coverage-dash.git python app.py`

That line should work on its own but I have had issues with port forwarding. If the above generates a link but it cannot be opened, try:
`docker run -p 8050:8050 <docker ID> python app.py`
