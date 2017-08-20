from ej import db, session

class GeneDesc(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    uni_name = db.Column(db.String)
    description = db.Column(db.String)

    def __init__(self,id,uni_name,description):
    	self.id = id
    	self.uni_name = uni_name
    	self.description = description