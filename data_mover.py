from io import StringIO
from h3 import h3
import pandas as pd
import geopandas as gpd
from sqlalchemy.engine.url import URL
from sqlalchemy import create_engine
import gdal
import csv
# import matplotlib
import time
import geojson


def df_to_postgresql(df, table, db_conn_params, if_exists='replace', sep='\t', encoding='utf8'):

    db_url = URL(drivername='postgresql+psycopg2', host=db_conn_params['host'], database=db_conn_params['db_name'],
                 username=db_conn_params['name'], port=db_conn_params['port'], password=['pass'])

    db_engine = create_engine(db_url)
    # Create Table
    df[:0].to_sql(table, db_engine, if_exists=if_exists, index=False)
    print("{} table prepared".format(table))
    # Prepare data
    output = StringIO()
    df.to_csv(output, sep=sep, index=False, header=False, escapechar="\\", encoding=encoding, quoting=csv.QUOTE_NONE)
    output.seek(0)
    # Insert data
    connection = db_engine.raw_connection()
    print("Opened database successfully")
    insert_start_time = time.time()
    cursor = connection.cursor()
    cursor.copy_from(output, table, sep=sep, null='')
    connection.commit()
    print("insert finish %s seconds ---" % (time.time() - insert_start_time))
    connection.close()
    return


def postgresql_to_df(db_table_name, db_conn_params):

    db_url = URL(drivername='postgresql+psycopg2', host=db_conn_params['host'], database=db_conn_params['db_name'],
                 username=db_conn_params['name'], port=db_conn_params['port'], password=['pass'])
    db_engine = create_engine(db_url)
    # Load the data
    df = pd.read_sql_table(db_table_name, db_engine)
    return df



