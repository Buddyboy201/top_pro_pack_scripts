
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String


def init_db(engine):
    meta = MetaData()
    cliques_table = Table(
        'cliques', meta,
        Column("id", Integer, primary_key=True),
        Column("size", Integer),
        Column("clique", String),
        Column("resid", String),
        Column("oldresid", String),
        Column("layerinfo", String),
        Column("pdbname", String)
    )
    meta.create_all(engine)

    conn = engine.connect()
    return conn, meta, engine

