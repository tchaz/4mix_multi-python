
4mix_multi in python 3.

Credits: (lifted from 4mix_multi.R:) CC BY-NC-SA original code by Eurogenes DESEUK1, modified by AJ14, 
         python translation & adaptation of the readme by tchaz 

Quick instructions for 4mix_multi.py:

Cliff notes: syntax for get_mix() is very similar to the original 4mix_multi.R - but no need to specify an output filename,
if you ran that and have python 3(.4) on your machine you should be good to go.


1. Make sure you are in the directory with the 4mix_multi.py file (or edit the file to change the default working directory).
   Make sure you have write priveleges for that directory (ie that you can create a file there) in order to be able to write the results file.
   Run the script 4mix_multi.py under python 3. (Tested under python 3.4.)
   

2. To run the script, there are seven arguments (the last optional). You write:
   "getMix(input file, target file, pop1, pop2, pop3, pop4, output filename prefix)"
   
   input file  = input file name (file CSV format)
   target file = target file name (file in CSV format)
   
   pop1 = The 1st population you want to use to model the target
   pop2 = The 2nd population you want to use to model the target
   pop3 = The 3rd population you want to use to model the target
   pop4 = The 4th population you want to use to model the target

   optional: output filename prefix (output file will in CSV format) 
             the script tries to make a sensible filename containing the target populations
             - you can add things to that name here 
             (or edit the script if, for example, you don't want the target pop names)


   So for example, if the input file is K8avg.csv, the target file is target.csv (with various European populations as the target groups), 
   the 1st population is HungaryGamba_EN, 
   the 2nd population is HungaryGamba_HG, 
   the 3rd population is Karelia_HG, 
   and the 4th population is Corded_Ware_LN, you write:

   getMix('K8avg.csv','target.csv','PPNB','La_Brana-1','HungaryGamba_HG','Karelia_HG')


3. After it finishes running (perhaps half a second for a target list of length ~200), 
   you'll get the closest solution along with the Euclidean distance between each target and the solution. 
   In this example, using the file target.csv we get, in the results_PPNB-La_Brana-1-HungaryGamba_HG-Karelia_HG.csv file 
   : 
  
   Population,PPNB,La_Brana-1,HungaryGamba_HG,Karelia_HG,D statistic
   Ashkenazi,69,6,0,25,0.026
   Basque_French,47,0,39,14,0.0078
   Belorussian,33,0,26,41,0.0254
   Bosnian,46,0,20,34,0.0233
   Bulgarian,52,0,16,32,0.0163
   Central_Sicilian,68,0,11,21,0.025
   ...
   
  The input file, and target file can be directly created from spreadsheets 
  saved in CSV format, and the output file which is in CSV format can be 
  direct imported into a spreadsheet.


4. If you run get_mix() without passing it arguments it will run 

   get_mix('K8avg.csv', 'target.csv', 'Bashkir', 'PPNB', 'Central_Greek', 'Corded_Ware_LN', 'results_') 
   
   - an arbitary default used during testing.

* NOTE ON THE INPUT FILE - The input file must be comma delimited,
  have one header row with the admixture component names, and have
  only one column of labels with the population and/or individual
  IDs.

* NOTE ON THE TARGET FILE - The format of the target file must be
  comma delimited, have admixture population headers, and have the
  target name to the left of the admixture values on the second
  and subsequent lines.


