import os, sys, json
from sqlalchemy import Column, ForeignKey, Integer, String, Float, Enum, TypeDecorator
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base()

class Json(TypeDecorator):
	"""Useful to turn seriable objects like lists into JSON objects"""
	impl = String
	def process_bind_param(self, value, dialect):
		return json.dumps(value)
	def process_result_value(self, value, dialect):
		return json.loads(value)

class Pdb(Base):
	__tablename__ = 'pdb'
	id = Column(Integer, primary_key=True) # SQL id 
	pdb_name = Column(String(250), nullable=False) 

class Interhelix_Interactions(Base):
	__tablename__ = 'interhelix_base_interactions'
	id = Column(Integer, primary_key=True)
	hbonds = Column(Json) # HO-atoms
	kihs = Column(Json)  # Knobs-Into-Holes 
	# Foreign key 
	pdb_id = Column(Integer, ForeignKey('pdb.id'))
	pdb = relationship(Pdb)

class Chain2Complex_BB_distance(Base):
	__tablename__ = 'com2com_distance'
	id = Column(Integer, primary_key=True)
	bb_distances = Column(Json)
	# Foreign key 
	pdb_id = Column(Integer, ForeignKey('pdb.id'))
	pdb = relationship(Pdb)

class HOLE_Output(Base):
	# HOLE pore conductance estimates (Gpred), different correction factors 
	__tablename__ = 'hole'
	id = Column(Integer, primary_key=True)
	HOLE_data = Column(Json)		# Non-corrected conductance
	# Foreign key 
	pdb_id = Column(Integer, ForeignKey('pdb.id'))
	pdb = relationship(Pdb)

if __name__ == "__main__":
	outdb = sys.argv[1] # Output name of database (.db)
	engine = create_engine('sqlite:///'+outdb)
	Base.metadata.create_all(engine)
