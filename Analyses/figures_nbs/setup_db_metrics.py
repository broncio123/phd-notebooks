import os
import sys
import json
from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base()

class Tags(Base):
    __tablename__ = 'tag'
    id = Column(Integer, primary_key=True) # SQL id 
    mutant = Column(String(250), nullable=False)
    group = Column(String(250), nullable=False)
    pdb_name = Column(String(250), nullable=False)
    frame = Column(String(250), nullable=False)

class Pore_Dimensions(Base):
    # Pore dimensions using HOLE and Isambard
    __tablename__ = 'pore_dimensions'
    id = Column(Integer, primary_key=True)
    #######################################
    pore_Rmin = Column(Float) # From HOLE
    pore_length = Column(Float) # From Isambard
    #######################################
    # Foreign key 
    tag_id = Column(Integer, ForeignKey('tag.id'))
    tag = relationship(Tags)
    
class Radii_of_Gyration(Base):
    # Pore dimensions using HOLE and Isambard
    __tablename__ = 'radii_of_gyration'
    id = Column(Integer, primary_key=True)
    #######################################
    Rg_x = Column(Float)
    Rg_y = Column(Float)
    Rg_z = Column(Float)
    #######################################
    # Foreign key 
    tag_id = Column(Integer, ForeignKey('tag.id'))
    tag = relationship(Tags)

if __name__ == "__main__":
    outdb = sys.argv[1] # Output name of database (.db)
    engine = create_engine('sqlite:///'+outdb)
    Base.metadata.create_all(engine)

