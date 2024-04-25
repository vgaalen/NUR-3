#!/bin/bash

echo "Run handin 3"

#echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Creating Plotting Directory!"
  mkdir plots
fi

#echo "Creating the output directory if it does not exist"
if [ ! -d "output" ]; then
  echo "Creating Output Directory!"
  mkdir output
fi

#echo "Creating the data directory if it does not exist"
if [ ! -d "data" ]; then
  echo "Creating Data Directory!"
  mkdir data
fi

echo "Download the datasets if they are not present"
if [ ! -e "data/satgals_m11.txt" ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m11.txt -P data/
fi
if [ ! -e "data/satgals_m12.txt" ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m12.txt -P data/
fi
if [ ! -e "data/satgals_m13.txt" ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m13.txt -P data/
fi
if [ ! -e "data/satgals_m14.txt" ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m14.txt -P data/
fi
if [ ! -e "data/satgals_m15.txt" ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m15.txt -P data/
fi

if [ ! -e "output/ex1a.txt" ] || [ ! -e "output/ex1b.txt" ] || [ ! -e "output/ex1c.txt" ] || [ ! -e "output/ex1d.txt" ]; then
  if [ -e "output/ex1a.txt" ] || [ -e "output/ex1b.txt" ] || [ -e "output/ex1c.txt" ] || [ -e "output/ex1d.txt" ]; then
    rm output/ex1*.txt
  fi
  echo "Run the script for 1"
  python3 ex1.py
fi

echo "Generating the pdf"
pdflatex NUR-3.tex > latex_output1.txt
#bibtex template.aux > bibtex_output.txt
pdflatex NUR-3.tex > latex_output2.txt
#pdflatex template.tex > latex_output3.txt


