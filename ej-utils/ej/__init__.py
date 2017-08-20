from flask import Flask, request, Response
from flask import render_template, send_from_directory, url_for
from flask.ext.sqlalchemy import SQLAlchemy

app = Flask(__name__)

app.config.from_object('ej.settings')
db = SQLAlchemy(app)
session = db.session
app.url_map.strict_slashes = False

import ej.models
import ej.controllers