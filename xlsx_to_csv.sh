#!/bin/sh
for i in *.xlsx; do xlsx2csv $i > $(echo $i | sed \"s/\.xlsx/\.csv/\") ; done
