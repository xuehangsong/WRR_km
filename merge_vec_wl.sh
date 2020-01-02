#!/bash/bin -l

for itime in $(seq 0 743)
do
    echo $(printf "%04d" $itime)
    echo "/mnt/e/john/vec/stage/"$itime".png"
    echo "/mnt/e/john/vec/arrow/vec."$(printf "%04d" $itime)".png"
    convert -append "/mnt/e/john/vec/stage/"$itime".png" "/mnt/e/john/vec/arrow/vec."$(printf "%04d" $itime)".png"  "/mnt/e/john/vec/merged/"$itime".jpg"
done
