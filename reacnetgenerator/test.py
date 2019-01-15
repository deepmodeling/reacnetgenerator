'''Test ReacNetGen'''


import os
import unittest

import reacnetgenerator
import requests


def download_file(url, local_filename):
    # from https://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
    r = requests.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
    return local_filename


class TestReacNetGen(unittest.TestCase):
    def test_reacnetgen(self):
        # download bonds.reaxc
        file_url = "https://dev.tencent.com/api/share/download/18504f9b-4995-40fd-8322-c56571dc9049"
        folder = "test"
        filename = "bonds.reaxc"
        bondfilename = os.path.join(folder, filename)

        if not os.path.exists(folder):
            os.makedirs(folder)
        print(f"Downloading {filename} ...")
        download_file(file_url, bondfilename)

        r = reacnetgenerator.ReacNetGenerator(inputfilename=bondfilename, atomname=[
                                              'H', 'O'], inputfiletype='lammpsbondfile', runHMM=True)
        r.runanddraw()

        print("Here are reactions:")
        with open(r.reactionfilename) as f:
            for line in f:
                print(line)
        self.assertTrue(isinstance(os.path.exists(r.resultfilename)))


if __name__ == '__main__':
    unittest.main()
