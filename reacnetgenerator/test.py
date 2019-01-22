'''Test ReacNetGen'''


import os
import unittest
import hashlib

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


def checkmd5(filename):
    if not os.path.isfile(filename):
        return
    myhash = hashlib.md5()
    with open(filename, 'rb') as f:
        while True:
            b = f.read(8096)
            if not b:
                break
        myhash.update(b)
    return myhash.hexdigest()


class TestReacNetGen(unittest.TestCase):
    def test_reacnetgen(self):
        # download bonds.reaxc
        file_url = "https://drive.google.com/uc?authuser=0&id=1CJ22BZTh2Bg3MynHyk_CVZl0rcpSQzRn&export=download"
        file_md5 = "d41d8cd98f00b204e9800998ecf8427e"

        folder = "test"
        filename = "bonds.reaxc"
        bondfilename = os.path.join(folder, filename)

        while not os.path.isfile(bondfilename) or checkmd5(bondfilename) != file_md5:
            if not os.path.exists(folder):
                os.makedirs(folder)
            print(f"Downloading {filename} ...")
            download_file(file_url, bondfilename)

        r=reacnetgenerator.ReacNetGenerator(inputfilename = bondfilename, atomname = [
                                              'H', 'O'], inputfiletype = 'lammpsbondfile', runHMM = True)
        r.runanddraw()

        print("Here are reactions:")
        with open(r.reactionfilename) as f:
            for line in f:
                print(line.strip())
        self.assertTrue(os.path.exists(r.resultfilename))


if __name__ == '__main__':
    unittest.main()
