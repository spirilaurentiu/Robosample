cp tools/flexor*py bin/
cp tools/process_flex.py bin/
cp tools/batstat.py bin/
rm bin/robosample*.so
cp buildpy/robosample.*.so bin/
echo "Installation complete."
