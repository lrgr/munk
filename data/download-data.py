# Script to download ppi-gi-data stored on the filestore
# Requires Python 3.2+
import os
import sys
import yaml
from urllib.request import urlretrieve

# This script dowloads all files from urls to specified file paths from supplied
# yaml file, where keys are file paths to save files to and values are 
# corresponding urls

# Usage:
#     python download-data <yaml-file>

def download_and_save(fp, url):
    print('DOWNLOADING from:', url)
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    urlretrieve(url, fp)
    print('SAVING to:', fp, '\n')

def main():
    with open(sys.argv[1], 'r') as IN:
        d = yaml.load(IN)
        for fp, url in d.items():
            download_and_save(fp, url)

if __name__ == '__main__':
    main()