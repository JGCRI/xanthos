import os
import requests



class InstallSupplement:

    def __init__(self):

        # full path to the xanthos root directory where the example dir will be stored
        self.root_dir = os.path.split(os.path.dirname(__file__))[0]

        self.save_zip = os.path.join(self.root_dir, 'example.zip')

        self.url = ''

    def fetch_zenodo(self):
        """Retrieve data supplement from Zenodo and save as zipped directory."""

        r = requests.get(self.url)

        with open(self.save_zip, 'wb') as out:
            out.write(r.content)

    def unzip_zenodo(self):
        """Unzip data supplement and destroy zipped directory."""
        #TODO: add in Python 2 and 3 support to read compressed data to memory and then to uncompressed file using
        #TODO:  StringIO and BytesIO respectively;
        #TODO: see (https://stackoverflow.com/questions/5710867/downloading-and-unzipping-a-zip-file-without-writing-to-disk)
        pass
