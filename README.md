To use the kd-tree:

Download the example dataset from here: https://rcdata.nau.edu/gowanlock_lab/datasets/real/iono_57min_5.16Mpts_2D.txt

Download the files in the kd-tree directory

Run Make

Run the program, the parameters are as follows: <input dataset> <input queryset> <# points dataset> <# points queryset> <epsilon/search distance> <number of data dimensions> 

For example:
./main iono_57min_5.16Mpts_2D.txt iono_57min_5.16Mpts_2D.txt 5159737 5159737 0.3 2

The output of this file is summarized in stats.txt, which will say:
Dataset/Queryset: iono_57min_5.16Mpts_2D.txt, iono_57min_5.16Mpts_2D.txt, NPOINTS/NQPOINTS: 5159737, 5159737, 0.3, 3890036889, Times (s): /Total/Total-Index/Index: 9.7816, 7.23186, 2.54974, NUMCPUTHREADS/IDIM/DTYPE: 32, 2, double, SHORTCIRCUIT: 1, RESULTSETSTORAGECPU: 2, COUNTONLY: 0



