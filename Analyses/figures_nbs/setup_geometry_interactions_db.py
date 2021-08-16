#!/usr/bin/python

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

class Tags(Base):
    __tablename__ = 'tag'
    id = Column(Integer, primary_key=True) # SQL id 
    mutant = Column(String(250), nullable=False)
    group = Column(String(250), nullable=False)
    pdb_name = Column(String(250), nullable=False)

class Interhelix_Interactions(Base):
    __tablename__ = 'interhelix_base_interactions'
    id = Column(Integer, primary_key=True)
    hbonds = Column(Json) # HO-atoms
    kihs = Column(Json)  # Knobs-Into-Holes 
    # Foreign key 
    tag_id = Column(Integer, ForeignKey('tag.id'))
    tag = relationship(Tags)

class RigidBody(Base):
    __tablename__ = 'rigid_body'
    id = Column(Integer, primary_key=True)
    # COM metrics
    com_bb = Column(Json)
    com_aa = Column(Json)
    com_primitives = Column(Json)
    com_assembly = Column(Json)
    # Euler angles from chain primitives
    euler_angles = Column(Json)
    # Foreign key 
    tag_id = Column(Integer, ForeignKey('tag.id'))
    tag = relationship(Tags)

class RadialProfiles(Base):
    __tablename__ = 'radial_profiles'
    id = Column(Integer, primary_key=True)
    primitive = Column(Json)
    punctual_bb = Column(Json)
    punctual_aa = Column(Json)
    vdw_bb_lower = Column(Json)
    vdw_bb_upper = Column(Json)
    vdw_aa_lower = Column(Json)
    vdw_aa_upper = Column(Json)
    # Foreign key 
    tag_id = Column(Integer, ForeignKey('tag.id'))
    tag = relationship(Tags)

class Miscellaneous(Base):
    __tablename__ = 'miscellaneous'
    id = Column(Integer, primary_key=True)
    vdw_minima = Column(Json)
    # Foreign key 
    tag_id = Column(Integer, ForeignKey('tag.id'))
    tag = relationship(Tags)

if __name__ == "__main__":
    outdb = sys.argv[1] # Output name of database (.db)
    engine = create_engine('sqlite:///'+outdb)
    Base.metadata.create_all(engine)