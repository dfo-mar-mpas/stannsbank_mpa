batches=30
for i in {1..30}
do
    # echo "loop $i"
    Rscript code/power_crab_function.R $batches $i &
done
