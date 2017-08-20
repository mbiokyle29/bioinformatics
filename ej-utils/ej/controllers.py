from flask import Flask, request, Response
from flask import render_template, url_for, redirect, send_from_directory, jsonify
from flask import send_file, make_response, abort, g
from ej import app, db, session

# default route
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/desc/<string:gene_name>')
def gene_desc(gene_name):
	print("OK")