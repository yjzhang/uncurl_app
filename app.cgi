#!/cse/web/homes/yjzhang/uncurl/bin/python

from wsgiref.handlers import CGIHandler
from app import app

activate_this = '/cse/web/homes/yjzhang/uncurl/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))

app.debug=True
CGIHandler().run(app)
