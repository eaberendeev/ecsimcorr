#!/bin/bash
cd ../Anime
mkdir Video

ffmpeg -r 4 -pattern_type glob -i './Res*.png' -y -c:v libx264 -r 25 -pix_fmt yuv420p -s 1024x720 ./Video/Video1.mp4
