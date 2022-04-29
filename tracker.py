# -*- coding: utf-8 -*-
#
#  Please refer to the Documentation for citation and acknowledgement
#


import sys
import numpy as np
from netCDF4 import Dataset


##
## input_variables
##
nlns   = 0                    # number of lines in input file
date   = np.empty([1],int)    # Date
time   = np.empty([1],int)    # Time
lat    = np.empty(1,float)    # Latitude
lon    = np.empty(1,float)    # Longitude
usteer = np.empty(1,float)  # Zonal steering velocity
vsteer = np.empty(1,float)  # Meridional steering velocity
speed  = np.empty(1,float)  # Wind speed at 850 hPa
topo   = np.empty(1,float)  # Topography
owz850 = np.empty(1,float)  # OWZ at 850 hPa
owz500 = np.empty(1,float)  # OWZ at 500 hPa
rh950  = np.empty(1,float)  # Relative humidity at 950 hPa
rh700  = np.empty(1,float)  # Relative humidity at 700 hPa
wsh    = np.empty(1,float)  # Smoothed 850 - 200 hPa wind shear
sh950  = np.empty(1,float)  # Specific humidity at 950 hPa



##
## clump_variables
##

cl_num              = 0                            # Number of clumps
cl_arsz             = 1000000                            # Clump array size increment
C_size              = np.empty([1],int)            # Number of gridpoints in the clump
C_size2             = np.empty([1],int)            # Number of gridpoints in the clump that satisfy all thresholds
C_land              = np.empty([1],int)            # Number of land points in the clump
C_date              = np.empty([1],int)            # Clump date
C_time              = np.empty([1],int)            # Clump time
C_flag              = np.empty([1],bool)           # Clump flag. False = free = not yet included in a CT string
C_ign_flag          = np.empty([1],bool)           # Flag for temporarily ignoring clumps in find_link
C_thrsh             = np.empty([1],bool)           # Overide threshold flag. True = overiding thrsh satisfied
C_lat               = np.empty([1],float)          # Clump latitude
C_lon               = np.empty([1],float)          # Clump longitude
C_latp1             = np.empty([1],float)          # Estimated clump latitude at next half time interval
C_lonp1             = np.empty([1],float)          # Estimated clump longitude at next halftime interval
C_latm1             = np.empty([1],float)          # Estimated clump latitude at previous half time interval
C_lonm1             = np.empty([1],float)          # Estimated clump longitude at previous half time interval
C_latp2             = np.empty([1],float)          # Estimated clump latitude at next time interval
C_lonp2             = np.empty([1],float)          # Estimated clump longitude at next time interval
C_latm2             = np.empty([1],float)          # Estimated clump latitude at previous time interval
C_lonm2             = np.empty([1],float)          # Estimated clump longitude at previous time interval
C_ust               = np.empty([1],float)          # Clump zonal steering velocity
C_vst               = np.empty([1],float)          # Clump meridional steering velocity
C_SPmax             = np.empty([1],float)          # Maximum value of 850 wind speed in the clump
C_OWZmax            = np.empty([1],float)          # Maximum vlaue of owz850 in the clump
C_owz850            = np.empty([1],float)          # Average owz850 in the clump
C_owz500            = np.empty([1],float)          # Average owz500 in the clump
C_rh950             = np.empty([1],float)          # Average rh950 in the clump
C_rh700             = np.empty([1],float)          # Average rh700 in the clump
C_wsh               = np.empty([1],float)          # Average wsh in the clump
C_sh950             = np.empty([1],float)          # Average sh950 in the clump
C_flagfail1_owz850  = np.empty([1],int)            # Count number of owz850 failures for invalid points in a clump
C_flagfail1_owz500  = np.empty([1],int)            # Count number of owz500 failures for invalid points in a clump
C_flagfail1_rh950   = np.empty([1],int)            # Count number of rh950 failures for invalid points in a clump
C_flagfail1_rh700   = np.empty([1],int)            # Count number of rh700 failures for invalid points in a clump
C_flagfail1_wsh     = np.empty([1],int)            # Count number of wsh failures for invalid points in a clump
C_flagfail1_sh950   = np.empty([1],int)            # Count number of sh950 failures for invalid points in a clump
C_flagfail2_OT_count= np.empty([1],int)            # Flag indicating inadequate number of valid points in a clump
C_flagfail2_owz850  = np.empty([1],int)            # Flag indicating clump mean owz850 value below threshold
C_flagfail2_owz500  = np.empty([1],int)            # Flag indicating clump mean owz500 value below threshold
C_flagfail2_rh950   = np.empty([1],int)            # Flag indicating clump mean rh950 value below threshold
C_flagfail2_sh950   = np.empty([1],int)            # Flag indicating clump mean sh950 value below threshold


##
## data_info
##

dlat      = 0.0        # Latitude increment (degrees)
dlon      = 0.0        # Longitude increment (degrees)
dtim      = 0.0        # Time increment (hours)
e_min     = 0          # Minimum number of neighbouring events for a clump to be considered
TC_min    = 0          # Minimum number of consecutive True links before TC declared
sea_min   = 0          # Minimum number of sea points to make a land influenced clump True
land_lim  = 0.0        # Topography value above which the grid point is considered to be land [m]
srch_rad  = 0.0        # Search radius to determine links in CT strings (km)
clmp_rad  = 0.0        # Clump radius, distance in which two CTs are to be combined
TH_owz850 = 0.0        # Overiding 850 hPa OWZ threshold (in case stricter thresholds are desired)
TH_owz500 = 0.0        # Overiding 500 hPa OWZ threshold
TH_rh950  = 0.0        # Overiding 950 hPa RH threshold
TH_rh700  = 0.0        # Overiding 700 hPa RH threshold
TH_wsh    = 0.0        # Overiding 850 - 200 hPa windshear threshold
TH_sh950  = 0.0        # Overiding 950 hPa SH threshold


##
## work_variables
##


T_event    =  np.empty([1],dtype="S")                         # Threshold event: N=new, U=unchecked, C=checked
clbkt_num  =  1000                       # Array dimension of clbkt
clbkt      =  np.zeros([clbkt_num],int)  # Clump bucket



##
## string_variables
##

smax            = 100000                     # Max number of event strings
lmax            = 1000                     # Max number of links in an event string
S_CT            = np.empty([1],int)     # Array of CT IDs for each link in each string
                                        #S_CT(string number,link number)
iTCdeclareflag  = np.empty([1],int)     # flag = 1 for all clumps at and after TC declaration
TCflag  = np.empty([1],bool)            # Flag = true for strings satisfying TC definition


##
## constants
##

Rearth  = 6.37e6            # Earth radius [m]
pi      = 2.0*np.arcsin(1.0)  # Circle constant
deg2rad = pi/180.0         # Degrees to radians conversion factor
rad2deg = 180.0/pi         # Radians to degrees conversion factor




def read_data(infile):
    global nlns, date,time,lat,lon,usteer,vsteer,speed,topo,owz850,owz500,rh950,rh700,wsh,sh950
    global T_event



    # Open input file
    file = open(infile, "r")
    # Count the number of lines in the input file
    nlns = sum(1 for line in file)

    file.seek(0, 0)


    # Allocate space for input threshold variables

    date = np.empty([nlns],int)
    time = np.empty([nlns],int)
    lat = np.empty([nlns],float)
    lon = np.empty([nlns],float)
    usteer = np.empty([nlns],float)
    vsteer = np.empty([nlns],float)
    speed  = np.empty([nlns],float)
    topo   = np.empty([nlns],float)
    owz850 = np.empty([nlns],float)
    owz500 = np.empty([nlns],float)
    rh950  = np.empty([nlns],float)
    rh700  = np.empty([nlns],float)
    wsh    = np.empty([nlns],float)
    sh950  = np.empty([nlns],float)


    # Allocate space for "work" threshold variables
    T_event =  np.empty([nlns],dtype="S")



    # Read input file
    count = 0
    for line in file:
        #print("LINE => " + str(line))
        date[count],time[count],lat[count],lon[count],usteer[count],vsteer[count],speed[count],topo[count],owz850[count],owz500[count],rh950[count],rh700[count],wsh[count],sh950[count]  = line.split(" ")

        # Assign character value of "N" (= new) to each threshold event to indicate
        # that the event is not yet included in a "clump".
        T_event[count] = "N"
        #print(str(usteer[count]) + " "+ str(vsteer[count]))
        #      + " " + str(speed[count]) +" "+str(topo[count]) + " " + str(round(owz850[count])) + " " + str(round(owz500[count])) + " " + str(rh950[count]) + " "+ str(rh700[count])+ " " + str(wsh[count]) + " "+str(sh950[count]))
        count = count + 1
    file.close()


#-----------END Read Data-------------------------------

def read_info(data_info_file):
    global dlon, dlat,dtim,e_min,TC_min,sea_min,land_lim,srch_rad,clmp_rad
    global TH_owz850,TH_owz500,TH_rh950,TH_rh700,TH_wsh,TH_sh950

    file      = open(data_info_file, "r")
    dlon      = float(file.readline().split(" ")[0])
    dlat      = float(file.readline().split(" ")[0])
    dtim      = float(file.readline().split(" ")[0])
    e_min     = int(file.readline().split(" ")[0])
    TC_min    = int(file.readline().split(" ")[0])
    sea_min   = int(file.readline().split(" ")[0])
    land_lim  = float(file.readline().split(" ")[0])
    srch_rad  = float(file.readline().split(" ")[0])
    clmp_rad  = float(file.readline().split(" ")[0])
    TH_owz850 = float(file.readline().split(" ")[0])
    TH_owz500 = float(file.readline().split(" ")[0])
    TH_rh950  = float(file.readline().split(" ")[0])
    TH_rh700  = float(file.readline().split(" ")[0])
    TH_wsh    = float(file.readline().split(" ")[0])
    TH_sh950  = float(file.readline().split(" ")[0])

    file.close()


#-----------END Read Info-------------------------------



def clump_allocate():
    global C_size,C_size2,C_land,C_date,C_time,C_lat,C_lon,C_latp1,C_lonp1,C_latm1,C_lonm1
    global C_latp2,C_lonp2,C_latm2,C_lonm2
    global C_ust,C_vst,C_SPmax,C_OWZmax,C_owz850,C_owz500
    global C_rh950,C_rh700,C_wsh,C_sh950,C_flag,C_ign_flag,C_thrsh,C_flagfail1_owz850, C_flagfail1_owz500
    global C_flagfail1_rh950, C_flagfail1_rh700,C_flagfail1_wsh, C_flagfail1_sh950,C_flagfail2_OT_count
    global C_flagfail2_owz850, C_flagfail2_owz500,C_flagfail2_rh950, C_flagfail2_sh950


    C_size              = np.zeros([cl_arsz],int)            # Number of gridpoints in the clump
    C_size2             = np.zeros([cl_arsz],int)            # Number of gridpoints in the clump that satisfy all thresholds
    C_land              = np.zeros([cl_arsz],int)            # Number of land points in the clump
    C_date              = np.zeros([cl_arsz],int)            # Clump date
    C_time              = np.zeros([cl_arsz],int)            # Clump time
    C_flag              = np.zeros([cl_arsz],bool)           # Clump flag. False = free = not yet included in a CT string
    C_ign_flag          = np.zeros([cl_arsz],bool)           # Flag for temporarily ignoring clumps in find_link
    C_thrsh             = np.zeros([cl_arsz],bool)           # Overide threshold flag. True = overiding thrsh satisfied
    C_lat               = np.zeros([cl_arsz],float)          # Clump latitude
    C_lon               = np.zeros([cl_arsz],float)          # Clump longitude
    C_latp1             = np.zeros([cl_arsz],float)          # Estimated clump latitude at next half time interval
    C_lonp1             = np.zeros([cl_arsz],float)          # Estimated clump longitude at next halftime interval
    C_latm1             = np.zeros([cl_arsz],float)          # Estimated clump latitude at previous half time interval
    C_lonm1             = np.zeros([cl_arsz],float)          # Estimated clump longitude at previous half time interval
    C_latp2             = np.zeros([cl_arsz],float)          # Estimated clump latitude at next time interval
    C_lonp2             = np.zeros([cl_arsz],float)          # Estimated clump longitude at next time interval
    C_latm2             = np.zeros([cl_arsz],float)          # Estimated clump latitude at previous time interval
    C_lonm2             = np.zeros([cl_arsz],float)          # Estimated clump longitude at previous time interval
    C_ust               = np.zeros([cl_arsz],float)          # Clump zonal steering velocity
    C_vst               = np.zeros([cl_arsz],float)          # Clump meridional steering velocity
    C_SPmax             = np.zeros([cl_arsz],float)          # Maximum value of 850 wind speed in the clump
    C_OWZmax            = np.zeros([cl_arsz],float)          # Maximum vlaue of owz850 in the clump
    C_owz850            = np.zeros([cl_arsz],float)          # Average owz850 in the clump
    C_owz500            = np.zeros([cl_arsz],float)          # Average owz500 in the clump
    C_rh950             = np.zeros([cl_arsz],float)          # Average rh950 in the clump
    C_rh700             = np.zeros([cl_arsz],float)          # Average rh700 in the clump
    C_wsh               = np.zeros([cl_arsz],float)          # Average wsh in the clump
    C_sh950             = np.zeros([cl_arsz],float)          # Average sh950 in the clump
    C_flagfail1_owz850  = np.zeros([cl_arsz],int)            # Count number of owz850 failures for invalid points in a clump
    C_flagfail1_owz500  = np.zeros([cl_arsz],int)            # Count number of owz500 failures for invalid points in a clump
    C_flagfail1_rh950   = np.zeros([cl_arsz],int)            # Count number of rh950 failures for invalid points in a clump
    C_flagfail1_rh700   = np.zeros([cl_arsz],int)            # Count number of rh700 failures for invalid points in a clump
    C_flagfail1_wsh     = np.zeros([cl_arsz],int)            # Count number of wsh failures for invalid points in a clump
    C_flagfail1_sh950   = np.zeros([cl_arsz],int)            # Count number of sh950 failures for invalid points in a clump
    C_flagfail2_OT_count= np.zeros([cl_arsz],int)            # Flag indicating inadequate number of valid points in a clump
    C_flagfail2_owz850  = np.zeros([cl_arsz],int)            # Flag indicating clump mean owz850 value below threshold
    C_flagfail2_owz500  = np.zeros([cl_arsz],int)            # Flag indicating clump mean owz500 value below threshold
    C_flagfail2_rh950   = np.zeros([cl_arsz],int)            # Flag indicating clump mean rh950 value below threshold
    C_flagfail2_sh950   = np.zeros([cl_arsz],int)            # Flag indicating clump mean sh950 value below threshold

#-----------END Clump Allocate-------------------------------


def S_CT_allocate():
    global S_CT,TCflag,iTCdeclareflag

    S_CT            = np.empty([smax, lmax],int)
    TCflag          = np.empty([smax],int)
    iTCdeclareflag  = np.empty([smax,lmax],int)

#-----------END S_CT Allocate-------------------------------



def centroid_threshold():

#*********************************************************************************
#
#  This subroutine finds clumps of threshold events, and replaces them with a
#  single set of values at a centroid location weighted by OWZ values.  The
#  threshold components will be specified either as an average or maximum/minimum
#  whichever is most appropriate.
#
#  The code takes the first threshold event, sets it to false and adds it to a
#  "clump" bucket, before it looks for immediate neighbours. Any neighbours found
#  are set to false and added to the clump bucket. New neighbours (with a .true.
#  setting) are sought for all events in the bucket, until no more are found.  At
#  this point the bucket contains all events in the clump, and single values
#  representing the clump are calculated.
#
#  The process is repeated for the next .true. event in the list, until all events
#  have been investigated.
#
#*********************************************************************************
    global T_event,clbkt, cl_num
    loop = False

    print("Inside centroid_threshold")
    print(nlns)
    for n in range(nlns):
        print("Searching for neighbours to threshold number:"+ str(n))
        if (T_event[n] == "C"):
            print("  ...included in a previous clump.")
        #! Empty the clump bucket and set loop exit flag
        clbkt = np.zeros([clbkt_num],int)  # Clump bucket
        print("****"+str(n)+"****")

        #print(clbkt)
        cb_count = 0
        loop = True

        print(T_event[n])
        if (T_event[n] == "N"):
            #! Set to unchecked and add to bucket
            T_event[n] = "U"
            #cb_count = cb_count + 1
            #print("cb_count "+str(cb_count))
            if(cb_count > clbkt_num) :
               print("Insufficient space allocated for clbkt array.")
               print("clbkt_num must be increased. Current value = " + str(clbkt_num))
               return

            clbkt[cb_count] = n

            #print("clbkt[cb_count]" + str(clbkt[cb_count]))

            while loop:
                 icnt = cb_count + 1
                 #print("Inside do while loop: icnt = "+str(icnt))
                 for i in range(icnt):
                     #print("inside inner loop"+str(i))
                     loop = False #! Exit do while loop when all events have been checked ("C")
                     #print(str(i) + " " + str(clbkt[i]) + " " + str(T_event[clbkt[i]]))
                     if (T_event[clbkt[i]] == "U") :
                         # Search for new neighbours, add them to bucket and update cb_count.
                         cb_count = neighbour_search(clbkt[i],cb_count,date[n],time[n])
                         T_event[clbkt[i]] = "C"
                         loop = True # Stay in do while loop until all events have been checked
                #end for-loop
            #end while-loop

            # If sufficient events in the bucket, calculate individual clump values
            # Remove this condition so that established circulations under shear that
            # have only one threshold point are included.  The e_min limit is now
            # applied in clump_values
            #If (cb_count >= e_min)
            #print("cb_count "+ str(cb_count))
            cl_num = clump_values(cb_count+1,cl_num)
            print("cl_num "+ str(cl_num))
            #endif
        #endif
    #end for-loop


    # Write clump information to output file
    file = open("Clump_out_py.txt", "w")
    file.write("**************************************\n")
    file.write("Number of clumps identified = "+str(cl_num)+"\n")

    print("**************************************")
    print("Number of clumps identified = "+str(cl_num))

    for k in range(1,(cl_num+1)) :
        print(str(C_date[k])+" "+ str(C_time[k])+" "+ str(C_lat[k])+" "+ str(C_lon[k])+" "+ str(C_size[k])+" "+ str(C_OWZmax[k])+" "+
        str(C_owz850[k])+" "+ str(C_owz500[k])+" "+ str(C_rh950[k])+" "+ str(C_rh700[k])+" "+ str(C_wsh[k])+" "+ str(C_sh950[k])+" "+ str(k)+" "+ str(C_thrsh[k])+"\n")

        file.write(str(C_date[k])+" "+ str(C_time[k])+" "+ str(C_lat[k])+" "+ str(C_lon[k])+" "+ str(C_size[k])+" "+ str(C_OWZmax[k])+" "+
        str(C_owz850[k])+" "+ str(C_owz500[k])+" "+ str(C_rh950[k])+" "+ str(C_rh700[k])+" "+ str(C_wsh[k])+" "+ str(C_sh950[k])+" "+ str(k)+" "+ str(C_thrsh[k])+"\n")
    print("**************************************")

    file.close()

#---------------End centroid_threshold---------------------------

def neighbour_search(cen,cb_count,cdate,ctime):
    global latp,latm,lonp,lonm, T_event,clbkt

    #print("Inside neighbour_search")
    #print(str(cen)+", " +str(date[cen])+ ", "+ str(time[cen])+ ", "+ str(lat[cen])+ ", "+ str(lon[cen])+ ", "+ str(cb_count))

    # Set neighbouring lat and lon values to 10% greater magnitude than actual values
    #  to avoid possible round-off issues with inequalities
    latp = lat[cen] + dlat*1.1
    latm = lat[cen] - dlat*1.1
    lonp = lon[cen] + dlon*1.1
    lonm = lon[cen] - dlon*1.1



    #print("nlns " + str(nlns) + " / cb_count "+ str(cb_count))
    # Loop through all events and look for current "New" neighbours.  If found set
    #  the neighbour to "Unchecked", add 1 to cb_count, and add to the clump bucket.
    for n in range(nlns):
        if ( (T_event[n] == "N") and (date[n] == cdate) and (time[n] == ctime) ) :
            if ( (lat[n] < latp) and (lat[n] > latm) and (lon[n] < lonp) and (lon[n] > lonm) ) :
                T_event[n] = "U"
                cb_count = cb_count + 1
                if(cb_count > clbkt_num) :
                   print("Insufficient space allocated for clbkt array.")
                   print("clbkt_num must be increased. Current value = "+ str(clbkt_num ))
                   return

                clbkt[cb_count] = n
                #print( "increment cb_count "+ str(cb_count))


    return cb_count

#------------------neighbour_search------------------------


def clump_values(cb_count,cl_num_2):
   global C_size,C_date,  C_time, C_flagfail1_owz850,  C_flagfail1_owz500,  C_flagfail1_rh950,  C_flagfail1_rh700,  C_flagfail1_wsh
   global C_flagfail1_sh950,  C_flagfail2_OT_count,  C_flagfail2_owz850,  C_flagfail2_owz500,  C_flagfail2_rh950,  C_flagfail2_sh950
   global C_lat, C_lon, C_ust, C_vst, C_SPmax, C_owz850, C_owz500, C_rh950, C_rh700, C_wsh, C_sh950, C_OWZmax, C_land
   global clbkt
   global C_size2,C_latp1,C_lonp1,C_latm1,C_lonm1
   global C_latp2,C_lonp2,C_latm2,C_lonm2


    # Local variabels
   m = 0         # ID number for T_event in the clump bucket
   Xdist = 0.0
   Ydist = 0.0 # Estimated distance [m] clump will move in one time interval
   del_lat = 0.0
   del_lon = 0.0 # Estimated distance (degrees) clump will move in one time interval
   Test_OWZmax = 0.0  # Experimental maximum OWZ parameter

   O_TH = np.zeros([cb_count], dtype=bool)     # Overide threshold flag
   OT_count = 0 # Counter for points that satisfy the overide threshold

   # Update clump_count
   cl_num_2 = cl_num_2 + 1

   # Check for sufficient room in clump arrays and reallocate if not
   if (cl_num_2 > cl_arsz) :
      print("Increasing the size of clump arrays")
      print("Current number of clumps = "+ str(cl_num_2))
      print("Current clump array size = "+str(cl_arsz))
     #     call reallocate_clump
      return            # Temporary instruction until subroutine reallocate_clump is completed

    # Set simple clump values
   C_size[cl_num_2] = cb_count
   C_date[cl_num_2] = date[clbkt[0]]
   C_time[cl_num_2] = time[clbkt[0]]

   # Set fail flags to zero
   C_flagfail1_owz850[cl_num_2]=0
   C_flagfail1_owz500[cl_num_2]=0
   C_flagfail1_rh950[cl_num_2]=0
   C_flagfail1_rh700[cl_num]=0
   C_flagfail1_wsh[cl_num_2]=0
   C_flagfail1_sh950[cl_num_2]=0
   C_flagfail2_OT_count[cl_num_2]=0
   C_flagfail2_owz850[cl_num_2]=0
   C_flagfail2_owz500[cl_num_2]=0
   C_flagfail2_rh950[cl_num_2]=0
   C_flagfail2_sh950[cl_num_2]=0

# Set averaged and maximized clump values
   C_lat[cl_num_2] = 0.0
   C_lon[cl_num_2] = 0.0
   C_ust[cl_num_2] = 0.0
   C_vst[cl_num_2] = 0.0
   C_SPmax[cl_num_2] = 0.0
   C_owz850[cl_num_2] = 0.0
   C_owz500[cl_num_2] = 0.0
   C_rh950[cl_num_2] = 0.0
   C_rh700[cl_num_2] = 0.0
   C_wsh[cl_num_2] = 0.0
   C_sh950[cl_num_2] = 0.0
   C_OWZmax[cl_num_2] = 0.0
   C_land[cl_num_2] = 0


   for n in range(cb_count) :
         m = clbkt[n]
         C_lat[cl_num_2] = C_lat[cl_num_2] + lat[m]
         C_lon[cl_num_2] = C_lon[cl_num_2] + lon[m]
         #print("n = "+str(n) + " clbkt[n] = " +str(clbkt[n])+" m = "+str(m)+" C_ust= "+str(C_ust[cl_num_2]) + " + "+ str(usteer[m]) +" = "+ str(C_ust[cl_num_2] + usteer[m]) +"\n")

         C_ust[cl_num_2] = C_ust[cl_num_2] + usteer[m]

         C_vst[cl_num_2] = C_vst[cl_num_2] + vsteer[m]
         if(C_SPmax[cl_num_2] < speed[m]):
              C_SPmax[cl_num_2] = speed[m]
         C_owz850[cl_num_2] = C_owz850[cl_num_2] + owz850[m]
         C_owz500[cl_num_2] = C_owz500[cl_num_2] + owz500[m]
         C_rh950[cl_num_2] = C_rh950[cl_num_2] + rh950[m]
         #print(str(C_rh950[cl_num_2]) + " + "+ str(rh950[m]) +" = "+ str(C_rh950[cl_num_2] + rh950[m]) +"\n")
         C_rh700[cl_num_2] = C_rh700[cl_num_2] + rh700[m]
         C_wsh[cl_num_2] = C_wsh[cl_num_2] + wsh[m]
         C_sh950[cl_num_2] = C_sh950[cl_num_2] + sh950[m]
    #     Test_OWZmax = owz850[m] + owz500[m]
    #     if(C_OWZmax[cl_num] < Test_OWZmax) C_OWZmax[cl_num] = Test_OWZmax
         if(C_OWZmax[cl_num_2] < owz850[m]):
              C_OWZmax[cl_num_2] = owz850[m]
    #  C_land now calculated below, so that only points satisfying the OT are considered
    #     if(topo[m] > 0.1) C_land[cl_num] = C_land[cl_num] + 1
         #print(str(C_lat[cl_num_2]) +" " + str(C_lon[cl_num_2])  +" " + str(C_ust[cl_num_2]) +" " + str(C_vst[cl_num_2]) +" " + str(C_SPmax[cl_num_2]) +" " + str(C_owz850[cl_num_2]) +" " + str(C_owz500[cl_num_2]) +" " + str(C_rh950[cl_num_2])+" " + str(C_rh700[cl_num_2])+" " + str(C_wsh[cl_num_2])+" " + str(C_sh950[cl_num_2]) + "\n")

   C_lat[cl_num_2] = C_lat[cl_num_2] / cb_count
   C_lon[cl_num_2] = C_lon[cl_num_2] / cb_count
   C_ust[cl_num_2] = C_ust[cl_num_2] / cb_count
   C_vst[cl_num_2] = C_vst[cl_num_2] / cb_count

   #print( str(cl_num_2) + " => " + str(C_ust[cl_num_2]) + " / " + str(cb_count) + " = " + str(C_ust[cl_num_2] / cb_count) )
   #print( str(cl_num_2) + " => " + str(C_vst[cl_num_2]) + " / " + str(cb_count) + " = " + str(C_vst[cl_num_2] / cb_count) )

   C_owz850[cl_num_2] = C_owz850[cl_num_2] / cb_count
   C_owz500[cl_num_2] = C_owz500[cl_num_2] / cb_count
   C_rh950[cl_num_2] = C_rh950[cl_num_2] / cb_count
   C_rh700[cl_num_2] = C_rh700[cl_num_2] / cb_count
   C_wsh[cl_num_2] = C_wsh[cl_num_2] / cb_count
   C_sh950[cl_num_2] = C_sh950[cl_num_2] / cb_count

# Estimate future and past clump positions
   Xdist   = C_ust[cl_num_2]*3600.*dtim*0.5
   Ydist   = C_vst[cl_num_2]*3600.*dtim*0.5
   del_lat = Ydist*rad2deg/( Rearth )
   del_lon = Xdist*rad2deg/( Rearth*np.cos(C_lat[cl_num_2]*deg2rad) )


   #print(str(cl_num_2) + " Xdist => " + str(C_ust[cl_num_2]) + "*3600*"+str(dtim) + "*0.5 = " + str(Xdist))
   #print(str(cl_num_2) + " Ydist => " + str(C_vst[cl_num_2]) + "*3600*"+str(dtim) + "*0.5 = " + str(Ydist))
   #print(str(cl_num_2) + " del_lat => " + str(Ydist) + "*"+str(rad2deg) + " / " + str(Rearth) + " = " + str(Ydist*rad2deg/( Rearth )))
   #print(str(cl_num_2) + " del_lon => " + str(Xdist) + "*"+str(rad2deg) + " / (" + str(Rearth) + " * cos(" + str(C_lat[cl_num_2]) + " * "+ str(deg2rad) + ") ) = " + str(Xdist*rad2deg/( Rearth*np.cos(C_lat[cl_num_2]*deg2rad) )))
   #print("cl_num_2 = " + str(cl_num_2))
   #print(str(C_lat[cl_num_2]) + " + " + str(del_lat) + " = " + str((C_lat[cl_num_2] + del_lat)))
   #print(str(C_lat[cl_num_2]) + " - " + str(del_lat) + " = " + str((C_lat[cl_num_2] - del_lat)))

   C_latp1[cl_num_2] = C_lat[cl_num_2] + del_lat
   C_latm1[cl_num_2] = C_lat[cl_num_2] - del_lat
   C_latp2[cl_num_2] = C_lat[cl_num_2] + del_lat*2.0
   C_latm2[cl_num_2] = C_lat[cl_num_2] - del_lat*2.0
   C_lonp1[cl_num_2] = C_lon[cl_num_2] + del_lon
   C_lonm1[cl_num_2] = C_lon[cl_num_2] - del_lon
   C_lonp2[cl_num_2] = C_lon[cl_num_2] + del_lon*2.0
   C_lonm2[cl_num_2] = C_lon[cl_num_2] - del_lon*2.0


   print(str(cl_num_2)+" => C_latp1[cl_num_2] = "+str(C_latp1[cl_num_2]) )
   print(str(cl_num_2)+" => C_latm1[cl_num_2] = "+str(C_latm1[cl_num_2]) )
   print(str(cl_num_2)+" => C_latp2[cl_num_2] = "+str(C_latp2[cl_num_2]) )
   print(str(cl_num_2)+" => C_latm2[cl_num_2] = "+str(C_latm2[cl_num_2]) )
   print(str(cl_num_2)+" => C_lonp1[cl_num_2] = "+str(C_lonp1[cl_num_2]) )
   print(str(cl_num_2)+" => C_lonm1[cl_num_2] = "+str(C_lonm1[cl_num_2]) )
   print(str(cl_num_2)+" => C_lonp2[cl_num_2] = "+str(C_lonp2[cl_num_2]) )
   print(str(cl_num_2)+" => C_lonm2[cl_num_2] = "+str(C_lonm2[cl_num_2]) )



# Determine if sufficient, neighbouring points satisfy the overiding thresholds
 # First, disregard the neighbouring points requirement and set to true if 'e_min'
 # points are satisfied.
   # Find points where overiding thresholds are satisfied
   O_TH =  np.zeros([cb_count], dtype=bool)
   OT_count = 0
   for n in range(cb_count):
      m = clbkt[n]
      if( (owz850[m] > TH_owz850) and (owz500[m] > TH_owz500) and (rh950[m] > TH_rh950)
             and (rh700[m] > TH_rh700) and (wsh[m] < TH_wsh) and (sh950[m] > TH_sh950) ) :
           O_TH[n] = True
           OT_count = OT_count + 1
           if(topo[m] > land_lim):
                C_land[cl_num_2] = C_land[cl_num_2] + 1
      else:
           if(owz850[m] <= TH_owz850):
                C_flagfail1_owz850[cl_num_2]=C_flagfail1_owz850[cl_num_2]+1
           if(owz500[m] <= TH_owz500):
                C_flagfail1_owz500[cl_num_2]=C_flagfail1_owz500[cl_num_2]+1
           if(rh950[m] <= TH_rh950):
                 C_flagfail1_rh950[cl_num_2]=C_flagfail1_rh950[cl_num_2]+1
           if(rh700[m] <= TH_rh700):
                 C_flagfail1_rh700[cl_num_2]=C_flagfail1_rh700[cl_num_2]+1
           if(wsh[m] >= TH_wsh):
                 C_flagfail1_wsh[cl_num_2]=C_flagfail1_wsh[cl_num_2]+1
           if(sh950[m] <= TH_sh950):
                 C_flagfail1_sh950[cl_num_2]=C_flagfail1_sh950[cl_num_2]+1

   C_size2[cl_num_2] = OT_count
# Initial method for determining true clumps was for at least two grid points to satisfy the
#  overiding threshold (5 lines immediately below).  This was giving poor results for many
#  clumps.  Average clump values gives a better result except for intense systems where storm
#  system shear is large, and 700 RH is occasionally suspect. To get around this problem
#  the threshold can be reassessed with average clump values excluding the shear and 700 RH.
   if ( (OT_count >= e_min) and (C_owz850[cl_num_2] > TH_owz850)
         and (C_owz500[cl_num_2] > TH_owz500) and (C_rh950[cl_num_2] > TH_rh950)
         and (C_sh950[cl_num_2] > TH_sh950) ) :
     C_thrsh[cl_num_2] = True
   else :
         C_thrsh[cl_num_2] = False
         if(OT_count < e_min):
             C_flagfail2_OT_count[cl_num_2]=1
         if(C_owz850[cl_num_2] <= TH_owz850):
              C_flagfail2_owz850[cl_num_2]=1
         if(C_owz500[cl_num_2] <= TH_owz500):
              C_flagfail2_owz500[cl_num_2]=1
         if(C_rh950[cl_num_2] <= TH_rh950):
              C_flagfail2_rh950[cl_num_2]=1
         if(C_sh950[cl_num_2] <= TH_sh950):
              C_flagfail2_sh950[cl_num_2]=1

# Write tracking comments
   #print("Clump"+str(cl_num_2)+" contains:")
   for n in range(cb_count):
     m = clbkt[n]
     print(str(date[m])+" " + str(time[m])+" "+str(lat[m])+" "+ str(lon[m])+ " "+ str(m)+" "+str(O_TH[n]))

   return cl_num_2




#------------------------------------------

def CT_strings():
   global C_lat,C_lon,C_lat,C_lon,clmp_rad,C_date,C_time,C_flag,C_thrsh
   global C_size2,C_OWZmax
   global S_CT,iTCdeclareflag,TCflag
   mj = 0              # m-jump, to jump forward in the m-loop
   mj2 = 0             # Second m-jump used in the recursive call of find_link
   s_num = 1           # String number
   l_num = 0           # Link number
   nxt_dt = 0
   nxt_tm = 0          # Date and time of next time level
   mret = 0            # The value of m to be returned from find_link
   Tcount = 0          # Counter for number of consecutive Ts in a string
   last_time = False   # Flag to ensure exit of recursive subroutine
   First1  = False     # Flag to ensure first TC declaration time only goes to output

   print("Inside CT_strings\n")
   #Initialize CT flags for each CT to False (= free) to indicate not included in a string
   C_flag = np.zeros([cl_arsz],bool)           # Clump flag. False = free = not yet included in a CT string
   S_flag = False
   dist = 0.0
# Look for nearby CTs and set the smaller/weaker ones to True (= taken) to exclude them from
# the search.  Or if only one satisfies the overiding threshold use it. If the CT's differ in
# size by one or less point, exclude the CT with a smaller
# OWXmax.  Otherwise, exclude the smaller CT (CT with less grid points).
#  Use S_flag to identify when a nearby CT has been found
   for n in range(cl_num+1):
     if (C_flag[n] == False ):
       for m in range(n+1,cl_num+1):
         if (C_date[n] < C_date[m] ):
              break  # Assumes CTs are ordered in time
         #print(str(C_date[n]) + " " + str(C_date[m]) + " " + str(C_time[n]) + " " + str(C_time[m]) +" "+ str(C_flag[m]))
         if ( (C_date[n]==C_date[m]) and (C_time[n]==C_time[m]) and (C_flag[m] == False) ) :
           print("How far apart are clumps"+str(n)+" and "+str(m)+"?\n")

           S_flag,dist = proximity_search(C_lat[n],C_lon[n],C_lat[m],C_lon[m],clmp_rad,S_flag,dist)

           if(S_flag):
                if( (C_thrsh[n]) and (C_thrsh[m] == False) ) :
                  C_flag[m] = True
                elif( (C_thrsh[m]) and (C_thrsh[n] == False) ) :
                  C_flag[n] = True
                else:
                    if(C_size2[n] > C_size2[m] + 1) :
                       C_flag[m] = True
                    elif(C_size2[n] < C_size2[m] - 1) :
                       C_flag[n] = True
                    else:
                       if(C_OWZmax[n] > C_OWZmax[m]) :
                         C_flag[m] = True
                       else:
                         C_flag[n] = True

                print("After proximity_search",C_OWZmax[n],C_OWZmax[m],C_flag[n],C_flag[m],n,m)


   print("******************************************\n")
   for n in range(cl_num):
     if(C_flag[n]):
         print("Ignoring clump:"+str(n)+"\n")

   print("******************************************\n")



# Initialize S_CT array to -999 to easily identify empty array locations
   S_CT.fill(-999)


 # Loop through CTs from each time period and search for a match with CT from next period
   print(C_date)
   print("Link search begins.\n")
   #print("cl_num - 2 = " + str(cl_num-2))
   #print("cl_num = " + str(cl_num))

   for n in range(1 , cl_num+1):
        S_flag = True
        l_num = 1

        # If still free assign CT number to the first position on the string
        if (C_flag[n] == False) :
          print("n = " +str(n)+" => Looking for a string beginning with clump:"+str(n) + "\n")

          S_CT[s_num][l_num] = n
          mj = n
          while S_flag:
                # write(6,*) "   (Inside do while loop)"
                # Calculate the next date/time
                nxt_dt,nxt_tm = next_time(C_date[mj],C_time[mj],nxt_dt,nxt_tm)
                print("After next_time.  Before "+ str(C_date[mj])+" "+str(C_time[mj])+ " After "+ str(nxt_dt) + str(nxt_tm))
                # Loop through remaining CT's looking for the next CT with nxt_dt, nxt_tm
                # Exit the loop when one is found and record the CT id
                last_time = False
                mj2 = 1
                mj,mj2,nxt_dt,nxt_tm,S_flag,l_num,s_num,mret,last_time = find_link(mj,mj2,nxt_dt,nxt_tm,S_flag,l_num,s_num,mret,last_time)
                #print(str(mj)+" "+str(mj2)+" "+str(nxt_dt)+" "+str(nxt_tm)+" "+str(S_flag)+" "+str(l_num)+" "+str(s_num)+" "+str(mret)+" "+str(last_time))
                #print("cl_num = " + str(cl_num))
                if(mret>cl_num):
                    break

        # Only update the string number if a string was identified in the do while loop
        if (l_num > 1):
            s_num = s_num + 1
        if (s_num > smax):
          #call reallocate_S_CT and TCflag iTCdeclareflag
          print("s_num exceeds array dimensions.\n")
          print("s_num = " + str(s_num) +" Array = "+str(smax))
          return



   # Remove the last string if it only has one CT
   if (l_num <= 1):
       s_num = s_num - 1



   # Search for strings that satisfy the TC genesis criteria
   TCflag.fill(False)
   iTCdeclareflag.fill(0)

   print("s_num = " + str(s_num))
   for n in range(1,s_num+1):
        Tcount = 0
        #print("TC test: String number "+str(n))
        for m in range(1,lmax+1):
          #print("Inside m loop"+str(m)+" "+str(C_thrsh[S_CT[n][m]])+" " + str(S_CT[n][m]))
          print("n = "+str(n)+" m ="+str(m)+" "+str(C_thrsh[S_CT[n][m]])+" " + str(S_CT[n][m])+"\n")

          #print("=> " + str(S_CT[n][m]) + "\n")

          if(S_CT[n][m] == -999):
              break
          if(C_thrsh[S_CT[n][m]]) :   # Current link is true.
            # Check for land influence
            if( (C_land[S_CT[n][m]] > 0) and (C_size2[S_CT[n][m]] - C_land[S_CT[n][m]] < sea_min) ) :
               Tcount = 0
            else:
               Tcount = Tcount + 1
               # Exit loop if Tcount has reached TC_min
        #       print("Tcount = "+str(Tcount) + " TC_min = " + str(TC_min))
               if(Tcount == TC_min) :
                    TCflag[n] = True
                    print(str(C_date[S_CT[n][m]]) + " " + str(nxt_dt))
                    iTCdeclareflag[n][m:lmax+1].fill(1)
                    break
               #endif
            #endif
            # Check for a time gap
            #write(6,*) "S_CT(n,m+1)",S_CT(n,m+1)
            if (S_CT[n][m] != -999) :
                  nxt_dt,nxt_tm = next_time(C_date[S_CT[n][m]],C_time[S_CT[n][m]],nxt_dt,nxt_tm)
                  #write(6,*) C_date(S_CT(n,m+1)),nxt_dt
                  if ( (nxt_dt != C_date[S_CT[n][m+1]] ) or (nxt_tm != C_time[S_CT[n][m+1]]) ) :
                     Tcount = 0
                  #endif
            #endif
          else:
              Tcount = 0
          #endif
        #enddo
   #enddo

   # Print to screen S_CT
   # Write S_CT to output file


   file_out  = open("S_CT_out_py.txt","w")
   file_out2 = open("S_CT_out2_py.txt","w")
   file_out3 = open("S_CT_out3_py.txt","w")

   print("#############################################")
   print(" Begin string information (S_CT)")
   print("#############################################")

   print(" lat , lon ,OW850,OW500, wshr,RH950,RH700,SH950,OWZmx,sz,lnd,clump, date ,tm, CoreThrsh")
   print("Core Thresh vals: " + str(TH_owz850)+" " + str(TH_owz500)+" " + str(TH_wsh)+" " + str(TH_rh950)+" " + str(TH_rh700)+" " + str(TH_sh950))
   file_out.write("#############################################\n")
   file_out.write(" Begin string information (S_CT)\n")
   file_out.write("#############################################\n")
   file_out.write(" lat , lon ,OW850,OW500, wshr,RH950,RH700,SH950,OWZmx,sz,lnd,clump, date ,tm, CoreThrsh\n")
   file_out.write("Core Thresh vals: " + str(TH_owz850)+" " + str(TH_owz500)+" " + str(TH_wsh)+" " + str(TH_rh950)+" " + str(TH_rh700)+" " + str(TH_sh950) + "\n")


   #if(any(TCflag))
   file_out3.write("Date0, Tim0,Lat0,Lon0,OWZmx0,Lat-12,Lon-12,Lat-24,Lon-24\n")
   print("s_num = "+ str(s_num))

   for n in range(1,s_num + 1):
       if(TCflag[n] == False) :
         #print("String number"+ str(n)+"\n")
         file_out.write("String number "+str(n) + "\n")
       else:
         #print("String number "+ str(n)+" *** TC declared Or test conditions met ***\n")
         file_out.write("String number "+str(n)+" *** TC declared Or test conditions met ***\n")

  #     First1=.true.
       print("lmax = "+ str(lmax))
       for m in range(1,lmax+1):
         if(S_CT[n][m] == -999):
             break

         #if(int(C_date[S_CT[n][m]]/100) == 200402):
        #     print(str(C_lat[S_CT[n][m]])+" " +str(C_lon[S_CT[n][m]]) +" "+ str(C_owz850[S_CT[n][m]])+" "+str(C_owz500[S_CT[n][m]]) + " " +str(C_wsh[S_CT[n][m]]) + " "+ str(C_rh950[S_CT[n][m]])+" "+str(C_rh700[S_CT[n][m]]) + " "+ str(C_sh950[S_CT[n][m]]) + " "+str(C_OWZmax[S_CT[n][m]]) + " "+ str(C_size2[S_CT[n][m]]) + " "+ str(C_land[S_CT[n][m]]) + " "+ str(S_CT[n][m]) + " "+str(C_date[S_CT[n][m]]) + " "+ str(C_time[S_CT[n][m]]/100) + " "+ str(C_thrsh[S_CT[n][m]])+"\n")

         file_out.write(str(C_lat[S_CT[n][m]])+" " +str(C_lon[S_CT[n][m]]) +" "+ str(C_owz850[S_CT[n][m]])+" "+str(C_owz500[S_CT[n][m]]) + " " +str(C_wsh[S_CT[n][m]]) + " "+ str(C_rh950[S_CT[n][m]])+" "+str(C_rh700[S_CT[n][m]]) + " "+ str(C_sh950[S_CT[n][m]]) + " "+str(C_OWZmax[S_CT[n][m]]) + " "+ str(C_size2[S_CT[n][m]]) + " "+ str(C_land[S_CT[n][m]]) + " "+ str(S_CT[n][m]) + " "+str(C_date[S_CT[n][m]]) + " "+ str(C_time[S_CT[n][m]]/100) + " "+ str(C_thrsh[S_CT[n][m]])+"\n")

         #write(50,131) C_lat[S_CT[n][m]],C_lon[S_CT[n][m]],C_owz850[S_CT[n][m]],C_owz500[S_CT[n][m]], &
        #               C_wsh[S_CT[n][m]],C_rh950[S_CT[n][m]],C_rh700[S_CT[n][m]],C_sh950[S_CT[n][m]], &
        #               C_OWZmax[S_CT[n][m]],C_size2[S_CT[n][m]],C_land[S_CT[n][m]],S_CT(n,m), &
        #               C_date[S_CT[n][m]],C_time[S_CT[n][m]]/100,C_thrsh[S_CT[n][m]]
         #chtm = ""
         chtm = C_time[S_CT[n][m]]/100#int2ch(chtm,2,C_time[S_CT[n][m]]/100,2)

         file_out2.write(str(m)+"\t"+str(n)+"\t"+str(iTCdeclareflag[n][m])+"\t"+str(round(C_lat[S_CT[n][m]],2))+"\t"+str(round(C_lon[S_CT[n][m]],2))+"\t"+
                        str(round(C_rh950[S_CT[n][m]],2))+"\t"+str(round(C_sh950[S_CT[n][m]],2))+"\t"+
                        str(round(C_OWZmax[S_CT[n][m]],2))+"\t"+str(round(C_SPmax[S_CT[n][m]],2))+"\t"+str(C_size[S_CT[n][m]])+"\t"+str(C_size2[S_CT[n][m]])+"\t"+
                        str(C_land[S_CT[n][m]])+"\t"+str(S_CT[n][m])+"\t"+str(C_date[S_CT[n][m]])+"\t"+str(chtm)+"\t"+str(C_thrsh[S_CT[n][m]])+"\t"+
                        str(C_flagfail1_owz850[S_CT[n][m]])+"\t"+str(C_flagfail1_owz500[S_CT[n][m]])+"\t"+
                        str(C_flagfail1_rh950[S_CT[n][m]])+"\t"+str(C_flagfail1_rh700[S_CT[n][m]])+"\t"+
                        str(C_flagfail1_wsh[S_CT[n][m]])+"\t"+str(C_flagfail1_sh950[S_CT[n][m]])+"\t"+
                        str(C_flagfail2_OT_count[S_CT[n][m]])+"\t"+
                        str(C_flagfail2_owz850[S_CT[n][m]])+"\t"+str(C_flagfail2_owz500[S_CT[n][m]])+"\t"+
                        str(C_flagfail2_rh950[S_CT[n][m]])+"\t"+str( C_flagfail2_sh950[S_CT[n][m]])+"\t"+
                        str(TH_owz850)+"\t"+str(TH_owz500)+"\t"+str(TH_rh950)+"\t"+str(TH_rh700)+"\t"+str(TH_wsh)+"\t"+str(TH_sh950)+"\n")
  #       if( (iTCdeclareflag(n,m) == 1).and.(First1) ) then
         if( (iTCdeclareflag[n][m] == 1) and (S_CT[n][m+1] == -999) ) :
               if( (C_thrsh[S_CT[n][m]]) and (C_thrsh[S_CT[n][m-1]]) and (C_thrsh[S_CT[n][m-2]]) ) :
                 file_out3.write(str(C_date[S_CT[n][m]])+" "+str(C_time[S_CT[n][m]]) + " " +str(C_lat[S_CT[n][m]])+" "+ str(C_lon[S_CT[n][m]])+" "+str(C_OWZmax[S_CT[n][m]])+" "+
                             str(C_lat[S_CT[n][m-1]]) +" "+str(C_lon[S_CT[n][m-1]])+" "+ str(C_lat[S_CT[n][m-2]])+" "+str(C_lon[S_CT[n][m-2]])+"\n")
      #         First1=.False.
          #endif
        #endif
       #enddo
     #enddo
   file_out.close()
   file_out2.close()
   file_out3.close()


#--------------End CT_strings----------------------------


def int2ch(character,l,integr,num_dig):
    #----------------------------------------------------------------------
    #  This subroutine converts an integer to a character.  It is
    #  based on "intXch" created by G.S. Dietachmayer September 1992.
    #
    #  charactr is the character string representation of the integer.
    #  l        is the desired length of the character string and must
    #           satisfy the following condition (l >= num_dig).
    #  integr   is the integer number to be converted.
    #  num_dig  is the number of digits in integr.  If (l > num_dig)
    #           charactr will be padded to the left with zeroes.
    #
    #  Created: 14-1-00 K. Tory.
    #----------------------------------------------------------------------

    epsilon=1e-6
    real_bit=0.0
    # Check that l is >= num_dig

    if(l < num_dig) :
         print('Error in int2ch. l must be >=  num_dig.')
         print('l, num_dig ='+str(l)+" "+str(num_dig))
         return
    #endif

    # Build up charactr digit by digit
    store = integr
    for i in range(l,1,-1):
         real_bit=float(store)*0.1+epsilon
         new_store=int(real_bit)
         digit = store - new_store*10
         character = chr(digit+48)
         store = new_store
       #enddo

    print('integer =',integr,'character = ',character)
    return character




def proximity_search(C_latc,C_lonc,C_latn,C_lonn,radius,S_flag,dist):

    coslat = np.cos( deg2rad*(C_latc+C_latn)*0.5 )
    Ydist  = (C_latc - C_latn)*deg2rad*Rearth
    Xdist  = (C_lonc - C_lonn)*deg2rad*Rearth*coslat
    dist   = np.sqrt( Xdist**2 + Ydist**2 +1.0e-5)

    #print(str(C_latc) + " - " + str(C_latn) + " * "+ str(deg2rad) +" * "+str(Rearth) + " = " + str((C_latc - C_latn)*deg2rad*Rearth))
    #print("coslat = " + str(coslat)+" Ydist = "+str(Ydist)+" Xdist = "+str(Xdist)+" "+str(dist)+"\n")

    if(dist < radius*1000.0):
        S_flag = True
    else:
        S_flag = False

    #print("dist = " + str(dist))
    # Convert dist from m to km
    dist = dist/1000.0

    #print("   Distance = "+str(dist)+"km. SR = "+ str(radius))
    #   write(6,*) "S_flag = ",S_flag
    #   write(6,*) "C_latc,C_lonc, C_latn, C_lonn",C_latc,C_lonc, C_latn, C_lonn
    return S_flag,dist

#--------------End proximity_search----------------------------



def next_time(in_dt,in_tm,out_dt,out_tm):


       hr1  = 0.0        # Input time real hours
       ihr1 = 0        # Input time integer hours
       imin1 =0       # Input time integer minutes
       hr2   = 0.0        # Output time real hours
       ihr2  =0       # Output time integer hours
       imin2  =0      # Output time integer minutes

       dy1 = 0
       dy2 = 0     # In and out days
       mn1 = 0
       mn2 = 0     # In and out months
       yr1 = 0
       yr2 = 0     # In and out years

       add_day = 0     # Number of days to be added to in_dt
       ldm  = 0        # Day number of the last day of the month



    # Calculate in_tm in hours (real)
       ihr1 = int( float(in_tm)/100.0)
       imin1= in_tm - ihr1*100
       hr1  = float(ihr1) + float(imin1)/60.0

    #   write(6,*) "ihr1, imin1, hr1 ",ihr1, imin1, hr1

    # Add dtim to hr to get out_tm
       hr2 = dtim + hr1
       add_day = 0
       while (hr2 >= 24.0):
         hr2 = hr2 - 24.0
         add_day = add_day + 1
         if(add_day > 10) :
           print("Problem with hr2 calculation")
           return

       ihr2 = int(hr2)
       imin2= int( (hr2 - float(int(hr2)) )*60.0)
       out_tm = ihr2*100 + imin2

    # Add days to in_dt
       if (add_day > 0) :
         # Get input yr, mn, dy
         yr1 = int( float(in_dt)/10000.0)
         mn1 = int( float(in_dt)/100.0)   # Temporary value
         dy1 = in_dt - mn1*100
         mn1 = mn1 - yr1*100

    #     write(6,*) "yr1, mn1, dy1",yr1, mn1, dy1

         # Determine last day of month
         if ((mn1==1) or (mn1==3) or (mn1==5) or (mn1==7) or (mn1==8) or (mn1==10) or (mn1==12)) :
           ldm = 31
         elif ((mn1==4) or (mn1==6) or (mn1==9) or (mn1==11)):
           ldm = 30
         elif (mn1==2) :
           if(np.mod(yr1,4)==0) :
             ldm = 29
           else:
             ldm = 28

           if((yr1==1900) or (yr1==2100) or (yr1==2200) or (yr1==2300)):
                ldm = 28


         # Find dy2,mn2,yr2
         dy2 = dy1 + add_day
         if (dy2 > ldm) :
           dy2 = dy2 - ldm
           mn2 = mn1 + 1
         else:
           mn2 = mn1

         if (mn2 > 12) :
           yr2 = yr1 + 1
           mn2 = mn2 - 12
         else:
           yr2 = yr1


         out_dt = 10000*yr2 + 100*mn2 + dy2
       else:
         out_dt = in_dt

       return out_dt,out_tm



def find_link(mj,mj2,nxt_dt,nxt_tm,S_flag,l_num,s_num,mret,last_time):

   global C_size,C_size2,C_land,C_date,C_time,C_lat,C_lon,C_latp1,C_lonp1,C_latm1,C_lonm1
   global C_latp2,C_lonp2,C_latm2,C_lonm2
   global C_ust,C_vst,C_SPmax,C_OWZmax,C_owz850,C_owz500
   global C_rh950,C_rh700,C_wsh,C_sh950,C_flag,C_ign_flag,C_thrsh,C_flagfail1_owz850, C_flagfail1_owz500
   global C_flagfail1_rh950, C_flagfail1_rh700,C_flagfail1_wsh, C_flagfail1_sh950,C_flagfail2_OT_count
   global C_flagfail2_owz850, C_flagfail2_owz500,C_flagfail2_rh950, C_flagfail2_sh950
   global dlon, dlat,dtim,e_min,TC_min,sea_min,land_lim,srch_rad,clmp_rad
   global TH_owz850,TH_owz500,TH_rh950,TH_rh700,TH_wsh,TH_sh950

   m = mj+mj2
   m_false  = False               # Stored 'm' when link to 'false' C_thrsh is found
   crnt_dt  = 0                   # Value of nxt_dt before recursive call
   crnt_tm  = 0                   # Value of nxt_tm before recursive call
   dist = 0.0
   dist_store = 0.0               # Distance between points (returned from proximity_search)
   sr1 = 200.0                    # Minimum search radius (km)
   sr_inc = 100.0                 # Search radius increment
   SR   = 0.0                     # Variable search radius
   s_rad  = 0.0                   # maximum search radius (as a function of latitude) 2013-07-30
   modlat = 0.0                   # Modulus of the latitude (degrees)
   repeat_srch = True             # Flag to determine if search is to be repeated over larger area
   double_exit = False            # Flag set to true when both do loops are to be exited

   # Search for links within an increasing radius, from sr1 to s_rad, at sr_inc increments
   SR = sr1
   repeat_srch = True
   C_ign_flag.fill(False)

   while(repeat_srch):
       repeat_srch = False

       m_false = 0
       for m in range(mj+mj2,cl_num+1):
             # Set variable search radius based on latitude
             modlat = np.sqrt(C_lat[mj]**2)
             if(modlat > 30.0) :
               s_rad = 400.0  # km
             elif(modlat < 15.0) :
               s_rad = 600.0  # km
             else:
               s_rad = 800.0 - 13.3333*modlat

             s_rad = s_rad*srch_rad/600.0
             if (last_time == False) :
               print("  Try clump "+str(m)+" Search radius =" + str(SR)+" "+str(s_rad))
             else:
               print("  Try clump "+str(m)+" one time-step later.  Search radius ="+ str(SR) +" "+str(s_rad))

             #print( "m = " + str(m) + " cl_num = " + str(cl_num))
             # If C_date > nxt_dt search for a link in the next time period then exit
             #  (Assumes CTs are ordered in time)
             # write(6,*) "C_date[m],C_time[m],nxt_dt,nxt_tm,m:",C_date[m],C_time[m],nxt_dt,nxt_tm,m
             if ( (C_date[m] > nxt_dt) or ( (C_date[m] == nxt_dt) and (C_time[m] > nxt_tm) ) or (m == cl_num) ) :
                   S_flag = False
                   SR = SR + sr_inc
                   if( (SR > s_rad) and (SR < s_rad+sr_inc-0.1) ):
                        SR = s_rad
                   if (SR <= s_rad) :
                        repeat_srch = True
                #   print("m = "+ str(m) + "cl_num = " + str(cl_num) +"   SR ="+str(SR)+" repeat_srch = "+str(repeat_srch))
                   if (m == cl_num) :
                     print("   Search stopped because there are no more clumps to search.\n")
                   else:
                     print("   Search stopped because clump",m,"is from a later time period.\n")

                   # If a link was found to a C_thrsh=false clump, add it now.
                   if ( (m_false != 0) and (repeat_srch == False) ) :
                         print("   No C_thrsh = true clumps found.  Add stored clump:"+str(m_false))
                         C_flag[m_false] = True  # T = taken, i.e., it will be excluded from further searches
                         l_num = l_num + 1
                         S_CT[s_num][l_num] = m_false
                        # print("    Clumps "+str(mj)+" "+str(m_false)+" are linked, on string "+str(s_num))
                         mj = m_false
                         if (l_num > lmax) :
                           #call reallocate_S_CT
                           print("l_num exceeds array dimensions.")
                           print("l_num = "+str(l_num)+" Array = "+ str(lmax))
                           return

                         nxt_dt,nxt_tm = next_time(nxt_dt,nxt_tm,nxt_dt,nxt_tm)
                         SR = sr1
                         double_exit = True
                         C_ign_flag.fill(False)
                         S_flag = True
                   #endif
                   break
            #endif

            #write(6,*) "C_date,C_time,nxt_dt,nxt_tm",C_date[m],C_time[m],nxt_dt,nxt_tm
             if ( (C_date[m] != nxt_dt) and (C_time[m] != nxt_tm) ):
                C_ign_flag[m] = True
             if ( (C_date[m] == nxt_dt) and (C_time[m] == nxt_tm) and (C_flag[m] == False) and (C_ign_flag[m] == False) ) :
                print("  How far apart are clumps"+str(mj)+" and "+ str(m)+"?")
                #write(6,*) "  C_ign_flag[m] = ",C_ign_flag[m]
                #****************************************
                # Test tracker without half-way storm position estimates.
                #    Comment out the lines between the stars below
                #
                #   call proximity_search(C_lat(mj),C_lon(mj),C_lat[m],C_lon[m],SR,S_flag,dist)
                #   dist_store = dist
                #**************************************

                #*************************************
                #   Use half-way storm position estimates
                #      Comment out the lines between the stars above

                #print("LatMJ "+str(C_latp1[mj-1])+" "+str(C_lonp1[mj-1])+" " +str(C_latp1[mj])+" "+str(C_lonp1[mj]))
                #print("LatM "+str(C_latm1[m-1])+" "+str(C_lonm1[m-1])+" " +str(C_latm1[m])+" "+str(C_lonm1[m]))

                dist_store = 0.0
                if ( last_time == False) :
                    S_flag,dist = proximity_search(C_latp1[mj],C_lonp1[mj],C_latm1[m],C_lonm1[m],SR,S_flag,dist)
                    if (dist < dist_store):
                        dist_store = dist
                else:
                    S_flag,dist = proximity_search(C_latp2[mj],C_lonp2[mj],C_latm2[m],C_lonm2[m],SR,S_flag,dist)
                    if (dist < dist_store):
                         dist_store = dist

                #*************************************
                if (dist_store > s_rad) :
                    C_ign_flag[m] = True
                    #write(6,*) "minimum dist > s_rad for clump",m
                    #write(6,*) "dist_store =",dist_store
                #endif
                # If link exists, add to string and exit both loops. If link does not exist S_flag will
                # be returned as false. It is reset to true before trying next CT
                double_exit = False
                if (S_flag) :
                    if (C_thrsh[m] == True) :
                        C_flag[m] = True  # T = taken, i.e., it will be excluded from further searches
                        l_num = l_num + 1
                        S_CT[s_num][l_num] = m
                        print("    Clumps"+str(mj)+" "+str(m)+" are linked, on string "+ str(s_num))
                        mj = m
                        if (l_num > lmax) :
                             #call reallocate_S_CT
                             print("l_num exceeds array dimensions.")
                             print("l_num = "+str(l_num)+" Array = "+ str(lmax))
                             return
                         #endif
                        nxt_dt,nxt_tm = next_time(nxt_dt,nxt_tm,nxt_dt,nxt_tm)
                        SR = sr1
                        double_exit = True
                        C_ign_flag.fill(False)
                        break
                    else:
                         # A link to a C_thrsh=false clump has been found. Add it later if no 'true' links found
                         # Add only the first one found for now.  May need to try something more clever later.
                         if (m_false == 0):
                               m_false = m
                         print("    Link found, but C_thrsh=false.  Link added later if no true C_thrsh found")
                    #endif
                else:
                    S_flag = True
                #endif
            #endif
       #end-for


       #print("*** M = " + str(m))
       if (double_exit) :
           break

       if ( (C_date[m] > nxt_dt) and (repeat_srch == False) and ( last_time == False) ) :
             last_time = True
             crnt_dt = nxt_dt
             crnt_tm = nxt_tm
             mj2 = m - mj
             nxt_dt,nxt_tm = next_time(crnt_dt,crnt_tm,nxt_dt,nxt_tm)
             mj,mj2,nxt_dt,nxt_tm,S_flag,l_num,s_num,mret,last_time = find_link(mj,mj2,nxt_dt,nxt_tm,S_flag,l_num,s_num,mret,last_time)
             nxt_dt = crnt_dt
             nxt_tm = crnt_tm
             if(S_flag == False) :
               print("   Second search stopped because clump"+str(m)+" is from a later time period.")
               break
             #endif
       #endif
   #enddo
   C_ign_flag.fill(False)

   mret = m
   print("*******************************  mret ="+str(mret))
   return mj,mj2,nxt_dt,nxt_tm,S_flag,l_num,s_num,mret,last_time



def start(inputFile):

    read_data(inputFile)

    print("reading data_info_file")
    read_info("data_info_file")


    # Allocate space for clump variables
    print("Before clump_allocate")
    clump_allocate()

# Find threshold clumps and calculate centroid threshold (CT)
#  values for each clump
    print("Before centroid_threshold")
    centroid_threshold()



    # Allocate space for S_CT array
    S_CT_allocate()

    # String together CTs in time
    CT_strings()


#if __name__ == "__main__":
#    main()
