#!/bash/bin -l

convert -delay 1 -resize 50% $(ls -tr /mnt/e/john/vec/merged/*.jpg)  /mnt/e/john/vec/vec_1.gif
