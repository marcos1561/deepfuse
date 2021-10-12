import pandas as pd
import urllib.request
import urllib.error
import os

#TODO: What to do with the label column?

def fetch(csv, output_dir, id_df_path, survey): # same ps as DECam
  """
  Download the galaxys from the file given in "csv" and store then in the "output_dir".
  
  Parameters:
    csv: .csv with the coordenates of the images to download.
    output_dir: The directory to save de images.
    id_df_path: is the path to .csv file that contains all galaxys already downloaded properly.
  This file will be used to verify if a given galaxy needs to be downloaded.
    survey: The survey where the images are downloaded.

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
  os.makedirs(output_dir, exist_ok=True)
  df = pd.read_csv(csv, index_col=0)

  if os.path.exists(id_df_path): 
    id_df = pd.read_csv(id_df_path, index_col=0)
  else:
    id_df = pd.DataFrame({'name':[], 'label':[]})

  for i, row in df.iterrows():
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
          print(f'\nredownloading ra={ra}, dec={dec}, ps={ps}')
          legacy_survey = f'https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer={survey}&pixscale={ps}'
          try:
            urllib.request.urlretrieve(legacy_survey, os.path.join(output_dir, candidate))
          except urllib.error.HTTPError as e:
            print(e)
            continue

          id_df.loc[i, 'name'] = candidate
          id_df.loc[i, 'label'] = label

          id_df['label'] = id_df['label'].astype(int)
          id_df.to_csv(id_df_path)
          print('csv saved')
            
    # download new candidate
    else:           
      print(f'\ndownloading ra={ra}, dec={dec}, ps={ps}')
      legacy_survey = f'https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer={survey}&pixscale={ps}'
      print(legacy_survey)
      try:
        urllib.request.urlretrieve(legacy_survey, os.path.join(output_dir, candidate))
      except urllib.error.HTTPError as e:
        print(e)
        continue
      except Exception as e:
        continue

      id_df.loc[i, 'name'] = candidate
      id_df.loc[i, 'label'] = label

      id_df['label'] = id_df['label'].astype(int)
      id_df.to_csv(id_df_path)
      print('csv saved')

##### TEST #####
train = "./data/teste.csv"
train_dir = "./data/train-ps"
id_train = "./data/id-train-ps.csv"
fetch(train, train_dir, id_train, survey="ls-dr9")
