import os
import requests
import sys
import zipfile

# get Python major version number
pyversion = sys.version_info.major

if pyversion <= 2:
    from StringIO import StringIO as BytesIO
else:
    from io import BytesIO as BytesIO


class InstallSupplement:

    # URL for DOI minted example data hosted on Zenodo
    URL = 'https://zenodo.org/record/2578287/files/example.zip?download=1'

    def __init__(self):

        # full path to the Xanthos root directory where the example dir will be stored
        self.root_dir = os.path.split(os.path.dirname(__file__))[0]

    def fetch_zenodo(self):
        """Retrieve data supplement from Zenodo and save as zipped directory."""

        # retrieve content from URL
        r = requests.get(InstallSupplement.URL)

        with zipfile.ZipFile(BytesIO(r.content)) as zipped:

            # extract each file in the zipped dir to the project
            for f in zipped.namelist():
                print("Unzipped: {}".format(os.path.join(self.root_dir, f)))
                zipped.extract(f, self.root_dir)
