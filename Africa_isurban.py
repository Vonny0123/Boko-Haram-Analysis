# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:18:46 2020

@author: ewand
"""

import pandas as pd
import geopandas as gpd
import tqdm
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import multiprocessing
import numpy as np
from pdl import pdl
import os
import sys
import time
import threading


def containment_tests(data, shapes, long_name='longitude', lat_name='latitude'):
  data = pd.DataFrame(data)
  points = gpd.GeoDataFrame(data.loc[:,[long_name,lat_name]], geometry=gpd.points_from_xy(data.loc[:,long_name], data.loc[:,lat_name])) #create a series of point objects representing location of events
  polys = shapes.geometry #This is a series of polygons
  containment_checker = polys.geometry.buffer(0).contains
  tqdm.tqdm.pandas(position=0, leave=True)
  r = points.geometry.progress_apply(containment_checker)
  return r.any(axis=1)

def multi_process_containment_tests(data, shapes, long_name='longitude', lat_name='latitude', cores=int(np.round(multiprocessing.cpu_count()/2))):
  data = pd.DataFrame(data)
  
  points = gpd.GeoDataFrame(data.loc[:,[long_name,lat_name]], geometry=gpd.points_from_xy(data.loc[:,long_name], data.loc[:,lat_name])) #create a series of point objects representing location of events
  polys = shapes.geometry #This is a series of polygons
  containment_checker = polys.geometry.buffer(0).contains
  with ProgressBar():  
    r = dd.from_pandas(points.geometry, npartitions=cores).map_partitions(lambda dframe: pd.Series(np.any(dframe.apply(containment_checker), axis=1)), meta=pd.Series(dtype=bool)).compute(scheduler='processes')  
  return r

class SpinnerThread(threading.Thread):

    def __init__(self):
        super().__init__(target=self._spin)
        self._stopevent = threading.Event()

    def stop(self):
        self._stopevent.set()

    def _spin(self):

        while not self._stopevent.isSet():
            for t in '|/-\\':
                sys.stdout.write(t)
                sys.stdout.flush()
                time.sleep(0.1)
                sys.stdout.write('\b')


if __name__ == '__main__':
  def yes_no(question):
    yes = set(['yes','y'])
    no = set(['no','n'])
     
    while True:
        choice = input(question).lower()
        if choice in yes:
           return True
        elif choice in no:
           return False
        else:
           print('Please respond with y/n.')
  
  def get_directory(question, error_message):
    dir_ = input(question)
    if os.path.exists(dir_):
      return dir_
    else:
      print(error_message)
      dir_ = get_directory(question, error_message)
      return dir_
          
  data_dir = get_directory('Please enter the directory containing your dataset, \neg. "C://Users/name/Desktop/Data/:"\n',
                           "That path doesn't exist, please enter the correct path.")

  data_filename = input('Please enter the name of your data file (csv), \neg. "data.csv":\n')
  multiprocess = yes_no('To speed up computation for large datasets\nmultiprocessing can be utilised. Should this be done? (y/n)\n')
  
  if multiprocess:
    cores = input('How many cores should be used? If left blank half will be used.\n')
    if cores == '':
      cores = int(np.round(multiprocessing.cpu_count()/2))
    else: cores = int(cores)
  
  long_name = input('Name of longitude column: (If "longitude", may leave blank)\n')
  if long_name == '':
    long_name = 'longitude'
    
  lat_name = input('Name of latitude column: (If "latitude", may leave blank)\n')
  if lat_name == '':
    lat_name = 'latitude'
    
  def load_data():
    global africapolis
    try: 
      africapolis = gpd.read_file(os.path.join(data_dir, 'Africapolis.shp'))
    except Exception:
      africapolis_url = 'http://www.africapolis.org/download/Africapolis_2015_shp.zip'
      pdl.download(africapolis_url, data_dir=data_dir, keep_download=False, overwrite_download=True, verbose=True)
      africapolis = gpd.read_file(os.path.join(data_dir,'Africapolis.shp'))

  print('Loading Data ')
  task = threading.Thread(target=load_data)
  task.start()

  spinner_thread = SpinnerThread()
  spinner_thread.start()

  task.join()
  data = pd.read_csv(os.path.join(data_dir, data_filename))
  spinner_thread.stop()
  
  
  print('Starting processing...')
  if multiprocess:
    isurban = multi_process_containment_tests(data=data, 
                                              shapes=africapolis,
                                              long_name=long_name,
                                              lat_name=lat_name,
                                              cores=cores)
  else:
    isurban = containment_tests(data=data, 
                                shapes=africapolis,
                                long_name=long_name,
                                lat_name=lat_name)
  
  data['is_urban'] = isurban
  print('Saving Data...')
  data.to_csv(os.path.join(data_dir, 'data_isurban.csv'))
  print('Done!')
