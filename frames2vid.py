#!/bin/python3

import glob
import sys
import re
import os


GLOB_PATH = "./frames/frame*.png"
RE_PATTERN = r"frame(\d+)"

def ParseInput():
    if len(sys.argv) == 1:
        output_file = "out.mp4"
    elif len(sys.argv) == 2:
        output_file = sys.argv[1]
    else:
        print(f"usage: {sys.argv[1]} output_file")
        sys.exit()
    return output_file

def PrepareFramesFolder():
    def NextPowerOf10(n):
        out = 1
        while out < n:
            out *= 10
        return out

    frame_files = glob.glob(GLOB_PATH)
    next_power_of_10 = NextPowerOf10(len(frame_files))

    for filename in frame_files:
        idx = int(re.search(RE_PATTERN, filename).group(1))
        prefix = str(next_power_of_10 + idx)[1:]

        os.system(f"mv {filename} ./frames/{prefix}-frame.png")

if __name__ == "__main__":
    # Discover the name of the output file
    output_file = ParseInput()

    # Reorganize the frames folder for ffmpeg
    PrepareFramesFolder()

    cmd = f"""
        ffmpeg -framerate 10 -pattern_type glob -i '*.png' \
        -c:v libx264 -pix_fmt yuv420p {output_file}
        """
    
    os.chdir("./frames")
    os.system(cmd)

