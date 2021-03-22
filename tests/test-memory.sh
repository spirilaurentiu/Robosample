#!/bin/bash

date_formatted=$(date '+%Y-%m-%d %H:%M:%S')
date_formatted="${date_formatted// /_}"
date_formatted=valgrind-report-$date_formatted.txt

valgrind --leak-check=full --show-leak-kinds=all ./tests/Robosample inp > $date_formatted 2>&1