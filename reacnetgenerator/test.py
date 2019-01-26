'''Test ReacNetGen'''


import hashlib
import json
import os
import unittest
import math

import pkg_resources

import reacnetgenerator
import requests
from tqdm import tqdm


def _download_file(url, local_filename):
    # from https://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
    r = requests.get(url, stream=True)
    total_size = int(r.headers.get('content-length', 0))
    block_size = 1024
    with open(local_filename, 'wb') as f:
        for chunk in tqdm(r.iter_content(chunk_size=1024), total=math.ceil(total_size//block_size), unit='KB', unit_scale=True, desc=f"Downloading {local_filename}..."):
            if chunk:
                f.write(chunk)
    return local_filename


def _checkmd5(filename):
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
    '''Test ReacNetGenerator'''

    def test_reacnetgen(self):
        ''' Test main process of ReacNetGen'''
        testparm = json.load(
            pkg_resources.resource_stream(__name__, 'test.json'))
        pathfilename = os.path.join(testparm['folder'], testparm['filename'])
        # download if not exists
        while not os.path.isfile(pathfilename) or _checkmd5(pathfilename) != testparm['md5']:
            try:
                os.makedirs(testparm['folder'])
            except OSError:
                pass
            print(f"Downloading  ...")
            _download_file(testparm['url'], pathfilename)

        r = reacnetgenerator.ReacNetGenerator(
            inputfilename=pathfilename, atomname=testparm['atomname'], inputfiletype=testparm['inputfiletype'], runHMM=testparm['hmm'])
        r.runanddraw()

        print("Here are reactions:")
        with open(r.reactionfilename) as f:
            for line in f:
                print(line.strip())
        self.assertTrue(os.path.exists(r.resultfilename))


if __name__ == '__main__':
    unittest.main()
