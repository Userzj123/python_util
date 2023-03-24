
#!/usr/bin/python

import sys
from pathlib import Path
import glob
import contextlib
from PIL import Image


def main():
    args = sys.argv
    if len(args) != 3:
        print('png2gif.py [filename] [varname]')
        sys.exit()
    else:
        filename = args[1]
        varname = args[2]
    
    ROOT_DIR = Path(filename).parent


    # filepaths
    fp_in = ROOT_DIR + "/%s_*.png" % varname
    fp_out = ROOT_DIR + "/%s.gif" % varname

    # use exit stack to automatically close opened images
    with contextlib.ExitStack() as stack:

        # lazily load images
        imgs = (stack.enter_context(Image.open(f))
                for f in sorted(glob.glob(fp_in)))

        # extract  first image from iterator
        img = next(imgs)

        # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
        img.save(fp=fp_out, format='GIF', append_images=imgs,
                save_all=True, duration=200, loop=0)