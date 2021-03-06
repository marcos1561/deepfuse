'''
  Baixa as imagens no legacy
'''

import pandas as pd
import urllib.request
import urllib.error
import os
from math import ceil

#TODO: What to do with the label column?

def fetch(csv, output_dir, id_df_path, survey, verbose=False):
  """
  Download the galaxys from the file given in "csv" and store then in the "output_dir".
  
  Parameters:
  ---------- 
  csv: str
    csv file with the coordenates of the images to download.
  
  output_dir: str 
    The directory to save de images.
  
  id_df_path: str
    The path to .csv file that contains all galaxys already downloaded properly.
    This file will be used to verify if a given galaxy needs to be downloaded.
  
  survey: str
    The survey where the images are downloaded.

  Return:
    error_occur: bool
      True if some error occurred during download, false otherwise.

  Exemple of usage::

    if __name__ ==  '__main__':
      train = '/content/drive/My Drive/Deepfuse/data/train.csv'
      train_dir = '/content/drive/My Drive/Deepfuse/data/train-ps'
      fetch(train, train_dir, '/content/drive/My Drive/Deepfuse/data/id-train-ps.csv', survey="ls-dr9")
      val = '/content/drive/My Drive/Deepfuse/data/val.csv'
      val_dir = '/content/drive/My Drive/Deepfuse/data/val-ps'
      fetch(val, val_dir, '/content/drive/My Drive/Deepfuse/data/id-val-ps.csv', survey="ls-dr9")
      test = '/content/drive/My Drive/Deepfuse/data/test.csv'
      test_dir = '/content/drive/My Drive/Deepfuse/data/test-ps'
      fetch(test, test_dir, '/content/drive/My Drive/Deepfuse/data/id-test-ps.csv', survey="ls-dr9")
  """
  error_occur = False

  os.makedirs(output_dir, exist_ok=True)
  df = pd.read_csv(csv, index_col=0)

  if os.path.exists(id_df_path): 
    id_df = pd.read_csv(id_df_path, index_col=0)
  else:
    id_df = pd.DataFrame({'name':[], 'label':[]})

  counter = 1
  max_counter = df.shape[0]
  count_step = ceil(1/10 * max_counter)

  for i, row in df.iterrows():
    counter += 1
    if counter % count_step == 0:
      print(f"Progress:{round((counter-1)/max_counter*100, 2)} %") 

    candidate = f'{i}.fits'
    ra = row['ra']
    dec = row['dec']  
    ps = row['ps']      
    label = row['label']

    # check if candidate is in the folder
    if os.path.exists(os.path.join(output_dir, candidate)): 
      # if candidate already exists, skip iteration
      try:
          if candidate in id_df.loc[i, 'name']:
              continue

      # if candidate is in the folder but not in id_df
      # it means it was not downloaded properly, so
      # we redownload it
      except Exception as e:
          os.remove(os.path.join(output_dir, candidate))
          if verbose:
            print(f'\nredownloading id={i} ra={ra}, dec={dec}, ps={ps}')
          legacy_survey = f'https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer={survey}&pixscale={ps}'
          try:
            urllib.request.urlretrieve(legacy_survey, os.path.join(output_dir, candidate))
          except urllib.error.HTTPError as e:
            print(e)
            error_occur = True
            continue
          except Exception as e:
            print(e)
            error_occur = True
            continue 

          id_df.loc[i, 'name'] = candidate
          id_df.loc[i, 'label'] = label

          id_df['label'] = id_df['label'].astype(int)
          id_df.to_csv(id_df_path)
          if verbose:
            print('csv saved')
            
    # download new candidate
    else:           
      legacy_survey = f'https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer={survey}&pixscale={ps}'
      if verbose:
        print(f'\ndownloading id={i} ra={ra}, dec={dec}, ps={ps}')
        print(legacy_survey)
      try:
        urllib.request.urlretrieve(legacy_survey, os.path.join(output_dir, candidate))
      except urllib.error.HTTPError as e:
        error_occur = True
        print(e)
        continue
      except Exception as e:
        print(e)
        error_occur = True
        continue

      id_df.loc[i, 'name'] = candidate
      id_df.loc[i, 'label'] = label

      id_df['label'] = id_df['label'].astype(int)
      id_df.to_csv(id_df_path)
      if verbose:
        print('csv saved')

  return error_occur

##### TEST #####
# cwd = os.getcwd()
# train = cwd + "\\data\\teste.csv"
# train_dir = cwd + "\\data\\train-ps"
# id_train = cwd + "\\data\\id-train-ps.csv"
# fetch(train, train_dir, id_train, survey="ls-dr9")
