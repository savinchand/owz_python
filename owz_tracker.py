import os, fnmatch
from calendar import monthrange
import thrsh as Threshold
import tracker as Tracker
import sys, getopt

indir = "/var/climatenas/ERA5/OWZ_Tracker/data"
outdir = ""
var1="rh"
var2="mrsh"
var3="temp"
var4="uwnd"
var5="vwnd"


#
#  Index number of level in netcdf file
#  Change this index number according to the data
#  e.g. ERA5 data have levels in assending order 200, 500,700,850,950
#  So I have set the index values accoring to the data
#
#
#  The script will rearrange the data in following format during preprocessing stage:
#  950, 850, 700, 500, 200
LVL200 = 1
LVL500 = 2
LVL700 = 3
LVL850 = 4
LVL950 = 5


topog_NH="topog/topogNH.nc"
topog_SH="topog/topogSH.nc"

#os.system("ls")

#**
# -  convert variable name to the name expected by threshold algorithm
# - extract levels data (950,850,700, 500, 200) for specified date for 00 and 12
# module assumes that levels are in sequence 200, 500, 700, 850, 950
# input :
#    varName - Name of the the variable
#    timestamp    - yyymmddhh format
#    timeIndex  - index value of 0 hour
#    inputFile - name of input file
# output: Temporary output file will generated with extrated data
##
def extractLevelData(varName, timestamp,timeIndex,inputFile):
    print("ncks -O -F -v "+ varName + " -d time,"+str(timeIndex)+" -d lvl,"+str(LVL200)+" "+indir+"/"+inputFile+" "+varName+"200_"+str(timestamp)+"_tmp.nc")
    os.system("ncks -O -F -v "+ varName + " -d time,"+str(timeIndex)+" -d lvl,"+str(LVL200)+" "+indir+"/"+inputFile+" "+varName+"200_"+str(timestamp)+"_tmp.nc")
    os.system("ncks -O -F -v "+ varName + " -d time,"+str(timeIndex)+" -d lvl,"+str(LVL500)+" "+indir+"/"+inputFile+" "+varName+"500_"+str(timestamp)+"_tmp.nc")
    os.system("ncks -O -F -v "+ varName + " -d time,"+str(timeIndex)+" -d lvl,"+str(LVL700)+" "+indir+"/"+inputFile+" "+varName+"700_"+str(timestamp)+"_tmp.nc")
    os.system("ncks -O -F -v "+ varName + " -d time,"+str(timeIndex)+" -d lvl,"+str(LVL850)+" "+indir+"/"+inputFile+" "+varName+"850_"+str(timestamp)+"_tmp.nc")
    os.system("ncks -O -F -v "+ varName + " -d time,"+str(timeIndex)+" -d lvl,"+str(LVL950)+" "+indir+"/"+inputFile+" "+varName+"950_"+str(timestamp)+"_tmp.nc")

    #Merge levels in sequence 950 to 200 for timestamp
    os.system("cdo merge "+varName+"950_"+timestamp+"_tmp.nc "+varName+"850_"+timestamp+"_tmp.nc "+varName+"700_"+timestamp+"_tmp.nc "+varName+"500_"+timestamp+"_tmp.nc "+varName+"200_"+timestamp+"_tmp.nc "+varName+"_"+timestamp+".nc")
    os.system("rm *_tmp.nc")



#
#   Append all variables into a single file
#
def appendFiles(timestamp):
    print("appending Files " + timestamp)
    #Append all variables of 00:00 hr into a single file
    os.system("ncks -A "+var2+"_"+timestamp+".nc "+ var1+"_"+timestamp+".nc")
    os.system("ncks -A "+var3+"_"+timestamp+".nc "+ var1+"_"+timestamp+".nc")
    os.system("ncks -A "+var4+"_"+timestamp+".nc "+ var1+"_"+timestamp+".nc")
    os.system("ncks -A "+var5+"_"+timestamp+".nc "+ var1+"_"+timestamp+".nc")

    os.system("rm "+var2+"_"+timestamp+".nc")
    os.system("rm "+var3+"_"+timestamp+".nc")
    os.system("rm "+var4+"_"+timestamp+".nc")
    os.system("rm "+var5+"_"+timestamp+".nc")



#
#  Split NH and SH data and save in separate files for further processing
#
def split_NH_SH_data(timestamp):
    print("split_NH_SH_data Files " + timestamp)

    #NH
    os.system("ncks -F -d lon,21,351 -d lat,86,151 "+var1+"_"+timestamp+".nc "+ var1+"_"+timestamp+"NH_tmp.nc")

    #SH
    os.system("ncks -F -d lon,21,351 -d lat,31,96 "+var1+"_"+timestamp+".nc "+ var1+"_"+timestamp+"SH_tmp.nc")

    os.system("rm "+var1+"_"+timestamp+".nc")



#
#   Append togography data to NH and SH files
#
def appendTopography(timestamp):
    print("appendTopography Files " + timestamp)

    #NH
    os.system("ncks -A "+topog_NH+" "+var1+"_"+timestamp+"NH_tmp.nc")
    #rename file
    os.system("mv "+var1+"_"+timestamp+"NH_tmp.nc owdata_"+timestamp+"_NH.nc")

    #SH
    os.system("ncks -A "+topog_SH+" "+var1+"_"+timestamp+"SH_tmp.nc")
    #rename file
    os.system("mv "+var1+"_"+timestamp+"SH_tmp.nc owdata_"+timestamp+"_SH.nc")



#
#   Concat individual theshold data of individual days of year into single file (OWZ2)
#   This file will be used as an input for the OWZ Tracker
#
def concatThresholdFiles(year):

    #Remove Old OWZ and SubJ Files
    os.system("rm OWZ2tracker_"+str(year)+"*.txt")
    os.system("rm SubJ_"+str(year)+"*.txt")


    #Concat TH files for this year to singlge OWZ2 file which will be used for tracking
    os.system("cat TH_001_"+str(year)+"*_NH*.txt >OWZ2tracker_"+str(year)+"_NH.txt")
    os.system("cat STJ_"+str(year)+"*_NH*.txt >SubJ_"+str(year)+"_NH.txt")

    os.system("cat TH_001_"+str(year)+"*_SH*.txt >OWZ2tracker_"+str(year)+"_SH.txt")
    os.system("cat STJ_"+str(year)+"*_SH*.txt >SubJ_"+str(year)+"_SH.txt")

    os.system("rm STJ* TH*")

#
#   Run threshold detector with preprocessed files for NH and SH
#
def runThresholdDetector(date,time, hem,file):

    timestamp = date+time

    if hem != "both":
        Threshold.process(file,date,time,hem)
    else:
        #NH
        Threshold.process("owdata_"+timestamp+"_NH.nc",date,time,"NH")

    	#SH
        Threshold.process("owdata_"+timestamp+"_SH.nc",date,time,"SH")
#   END OF runThresholdDetector




#
#   Run tracker on the threshold from OWZ2 file and save output in Clump, S_CT, S_CT2 and S_CT3
#
def runTracker(year,hem):
    Tracker.start("OWZ2tracker_"+str(year)+"_"+hem+".txt")
    os.system("mv Clump_out_py.txt "+outdir+"Clump_out_"+str(year)+"_"+hem+".txt")
    os.system("mv S_CT_out_py.txt "+outdir+"S_CT_out_"+str(year)+"_"+hem+".txt")
    os.system("mv S_CT_out2_py.txt "+outdir+"S_CT_out2_"+str(year)+"_"+hem+".txt")
    os.system("mv S_CT_out3_py.txt "+outdir+"S_CT_out3_"+str(year)+"_"+hem+".txt")
#   END OF runTracker






#
#
#
def preprocessData(fromYear, toYear):

    time00 = 1
    time12 = 2
    for year in range(fromYear,toYear+1) :
        for month in range(1,13):
            mm = '{:02d}'.format(month)
            for day in range(1,monthrange(year, month)[1]+1):
                dd = '{:02d}'.format(day)

                timestamp00= str(year)+mm+dd+"00"
                timestamp12= str(year)+mm+dd+"12"
                #Extract RH
                extractLevelData(var1,timestamp00,time00,"RH_"+str(year)+".nc")
                extractLevelData(var1,timestamp12,time12,"RH_"+str(year)+".nc")

                #Extract SPFH
                extractLevelData(var2,timestamp00,time00,"SPFH_"+str(year)+".nc")
                extractLevelData(var2,timestamp12,time12,"SPFH_"+str(year)+".nc")

                #Extract TEMP
                extractLevelData(var3,timestamp00,time00,"TEMP_"+str(year)+".nc")
                extractLevelData(var3,timestamp12,time12,"TEMP_"+str(year)+".nc")

                #Extract UWIND
                extractLevelData(var4,timestamp00,time00,"UWND_"+str(year)+".nc")
                extractLevelData(var4,timestamp12,time12,"UWND_"+str(year)+".nc")

                #Extract VWIND
                extractLevelData(var5,timestamp00,time00,"VWND_"+str(year)+".nc")
                extractLevelData(var5,timestamp12,time12,"VWND_"+str(year)+".nc")

                #Append all variables in single file for timestamp
                appendFiles(timestamp00)
                appendFiles(timestamp12)


                #Extract NH and SH data from source file
                split_NH_SH_data(timestamp00)
                split_NH_SH_data(timestamp12)

                #Append topography data to NH and SH
                appendTopography(timestamp00)
                appendTopography(timestamp12)

                #Detect threshold for NH and SH
                runThresholdDetector(str(year)+mm+dd,"00","both","")
                runThresholdDetector(str(year)+mm+dd,"12","both","")


                time00=time00 + 2
                time12=time12 + 2

        concatThresholdFiles(year)
        runTracker(year,"NH")
        runTracker(year,"SH")



##
def processTestData(dataPath):
    listOfFiles = os.listdir(dataPath)
    years = []
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, "*.nc"):
             timestamp = entry.split("_")[1]
             year = timestamp[0:4]
             if(year not in years):
                years.append(year)
             date = timestamp[0:8]
             time = str(timestamp[-2:])
             hem = entry.split("_")[2].split(".")[0]
	     path = ""
             if path != ".":
		path = dataPath
	        if path[-1] != "/":
		   path += "/"

             path = path+entry 
             print path,timestamp,date,time,hem
             runThresholdDetector(date,time,hem,path)

    runTrackerOnTestData(years)

def runTrackerOnTestData(years):
    for year in years:
        concatThresholdFiles(year)
        runTracker(year,"NH")
        runTracker(year,"SH")






def main(argv):
   global indir, outdir
   indir = ''
   outdir = ''
   yearFrom = 0
   yearTo = 0
   test = False
   testFiles = []
   if len(argv) == 0:
        print 'usage: owz_tracker.py -i <input directory> -o <output directory> -f <year from> - t <year to>'
        print 'For Testing preprocessed owd_data file, type --test followed by name of files separated by , '
        print 'E.g --test file1,file2,file3'
        sys.exit()

   if argv[0] == "--test":

         #Run on preprocessed owz_data files for testing
         testFiles = argv[1] #.split(",")
         print testFiles
         processTestData(testFiles)
#         testOWZ(testFiles)
   else:

       try:
          opts, args = getopt.getopt(argv,"hi:o:f:t:",["ifile=","ofile=","yearFrom=","yearTo="])
       except getopt.GetoptError:
          print 'owz_tracker.py -i <input directory> -o <output directory> -f <year from> - t <year to>'
          print 'For Testing preprocessed owd_data file, type --test followed by name of files separated by , '
          print 'E.g --test file1,file2,file3'

          sys.exit(2)


           # Run OWZ on dataset
       for opt, arg in opts:
              if opt == '-h':
                 print 'usage: owz_tracker.py -i <input directory> -o <output directory> -f <year from> - t <year to>'
                 sys.exit()
              elif opt in ("-i", "--ifile"):
                 indir = arg
              elif opt in ("-o", "--ofile"):
                 outdir = arg
              elif opt in ("-f", "--yearFrom"):
                  yearFrom = arg
              elif opt in ("-t", "--yearTo"):
                  yearTo = arg



       if yearFrom == 0 or yearTo == 0 or indir == '':
                print 'Invalid arguments'
                print 'usage: owz_tracker.py -i <input directory> -o <output directory> -f <year from> - t <year to>'
                sys.exit()

       if indir != "" and os.path.isdir(indir) == False:
                print "Input directory doesn't exist"
                sys.exit()

       if outdir != "" and os.path.isdir(outdir) == False:
                print "Output directory doesn't exist"
                sys.exit()

       preprocessData(int(yearFrom),int(yearTo))

if __name__ == "__main__":
    main(sys.argv[1:])
