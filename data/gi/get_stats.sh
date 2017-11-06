echo "Total number of interations"
wc -l $1 | awk '{print $1-1}'
echo "Total number of SLs"
grep -P "\tSL" $1 | wc -l

echo "Total number of Non-Sls"
grep Non-SL $1 | wc -l

echo "Total number of Unknowns"
grep Inco $1 | wc -l
