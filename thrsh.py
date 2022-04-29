# -*- coding: utf-8 -*-
#
#  Please refer to the Documentation for citation and acknowledgement
#

import sys
import numpy as np
from netCDF4 import Dataset



owz    = [] # Normalized Okubo-Weiss * abs vort (s^-1)
rh     = [] # Relative humidity (%)
mrsh   = [] # Mixing ratio or spec hum (In kg/kg, Out g/kg)
wsh    = [] # Smoothed wind shear (m/s)
avor   = [] # Absolute vorticity (s^-1)
thrsh  = []   # Combined threshold (value 1.0 when threshold met)
fld1   =  []  #Temporary threshold variable arrays
fld2   =  []  #Temporary threshold variable arrays
fld3   =  []  #Temporary threshold variable arrays
fld4   =  []  #Temporary threshold variable arrays
fld5   =  []  #Temporary threshold variable arrays
fld6   =  []  # Temporary threshold variable arrays
uwnd   =  [] # Zonal wind (m/s)
vwnd   =  [] # Meridional wind (m/s)
temp   =  [] # Air temperature (K)
topo3d = []                   # Topography (when determined from a 3-d field)
topo   =  []     # Topography

lat = []       # latitude (degrees)
lon = []      # longitude (degrees)
lvl = []     # pressure levels
date = []      # date (YYYYMMDD)
time = []      # time (HHMM)

xd = 0                        # length of longitude dimension
yd = 0                        # length of latitude dimension
zd = 0                        # length of height dimension
td = 0                        # length of time dimension

k950 = -1                     # pressure level integers
k850 = -1                     # pressure level integers
k700 = -1                     # pressure level integers
k500 = -1                     # pressure level integers
k200 = -1                     # pressure level integers
i1 = 0
il = 0
j1 = 0
jl = 0                  #  First and last i and j indices that define the output  window


Hem = "SH"


Lowz = True
Lwsh = True
Lrh  = True #            ! Flags=T if fields (OWZ, wind shear, relative humidity) in input file
Lmr = True
Lsh = True			  #! Flags=T if fields (mixing ratio, specific humidity) in input file
Ltop = True
Lsfg = True          #       ! Flags=T if fields (topography, surface geopotential) in input file
Lav  = True       # Flags=T if fields (absolute vorticity) in input file


num_th_comb = 0 #    ! Number of threshold combinations
Towz850 = 0.0 #          ! Threshold owz value for 850 hPa
Towz500 = 0.0 #         ! Threshold owz value for 500 hPa
Trh950  = 0.0 #         ! Threshold rh or sd value for 950 hPa
Trh700  = 0.0 #         ! Threshold rh or sd value for 700 hPa
Twsh    = 0.0 #        ! Threshold wsh value for 850 - 200 hPa wind shear
Tsh950  = 0.0 #         ! Threshold sh value for 950 hPa

Wltmin = 0.0  #            ! Minimum latitude for output Window
Wltmax = 0.0  #         ! Maximum latitude for output Window
Wlnmin = 0.0  #        ! Minimum longitude for output Window
Wlnmax = 0.0  #          ! Maximum longitude for output Window

num_smth = 0  #      ! Number of smoothing operations for wind shear

Lwrite  = False #      ! If true, write threshold or other arrays to input NetCDF file
Loutput = False      #! If true, write lat/lon of every thrsh satisfied gridpoint to ASCII file
Lsd   = False   # ! If true, RH arrays will carry saturation deficit instead of RH




#
#  Load configuration from input_thrsh_NH/SH file
#  Read the input file to determine: the "window of interest" in which the combined
#  thresholds are to be calculated

def input_read(input_thrsh_file):
    global num_th_comb,num_smth, Towz850, Towz850, Towz500, Trh950,Trh700, Twsh, Tsh950, Wltmin,Wltmax,Wlnmin ,Wlnmax,num_smth
    global Lwrite,Loutput,Lsd


    file = open(input_thrsh_file, "r")
    Wltmin = float(file.readline().split(" ")[0])   #Minimum latitude for output Window
    Wltmax = float(file.readline().split(" ")[0])   #Maximum latitude for output Window
    Wlnmin = float(file.readline().split(" ")[0])   #Minimum longitude for output Window
    Wlnmax = float(file.readline().split(" ")[0])   #Maximum longitude for output Window
    num_smth = int(file.readline().split(" ")[0])   #Number of smoothing operations for wind shear
    Lwrite = bool(file.readline().split(" ")[0])    #If true, write threshold or other arrays to input NetCDF file
    Loutput = bool(file.readline().split(" ")[0])   #If true, write lat/lon of every thrsh satisfied gridpoint to ASCII file
    num_th_comb = int(file.readline().split(" ")[0])   #  Number of threshold combinations


    Towz850 = np.zeros(num_th_comb)  #Threshold owz value for 850 hPa
    Towz500 = np.zeros(num_th_comb)  #Threshold owz value for 500 hPa
    Trh950  = np.zeros(num_th_comb)  #Threshold rh or sd value for 950 hPa
    Trh700  = np.zeros(num_th_comb)  #Threshold rh or sd value for 700 hPa
    Twsh    = np.zeros(num_th_comb)  #Threshold wsh value for 850 - 200 hPa wind shear
    Tsh950  = np.zeros(num_th_comb)  #Threshold sh value for 950 hPa


    for i in range(num_th_comb):
        line = file.readline()
        Towz850[i] = float(line.split(",")[0])
        Towz500[i] = float(line.split(",")[1])
        Trh950[i] = float(line.split(",")[2])
        Trh700[i] = float(line.split(",")[3])
        Twsh[i] = float(line.split(",")[4])
        Tsh950[i] = float(line.split(",")[5])

    print("Check threshold values")
    print("Towz850 "+str(Towz850)+", Towz500 "+str(Towz500)+", Trh950 "+str(Trh950)+", Trh700 "+str(Trh700)+", Twsh "+str(Twsh)+", Tsh950 "+str(Tsh950))

    print("If RH thresholds < 0, SD will be used instead.")
    if(Trh950[0] < 0 and Trh700[0] < 0):
        Lsd = True
        Lwrite = False
        print("Sturation deficit thresholds in use.")
    else:
        Lsd = False
        print("Relative humidity thresholds in use.")

    file.close()


#
#  Write to an ascii text file threshold information everywhere the combined
#  threshold is satisfied
#
def output(idcmb):
   print "** writing TH output ** "+date+" - "+time

   st_u  = 0.0
   st_v  = 0.0

   # Calculate yr and md
   store = int(date)/10000.
   yr = int(store)
   md = int(date) - yr*10000

   print("Date " + str(date))

   if(time[0] < 24):
       time[0] = time[0]*100

   avlat =  np.sum(lat)/yd

   # Create output file
   output_file = "TH_00"+str(idcmb+1)+"_"+str(date)+"_"+str(time)+"_"+Hem+"_00.txt"
   file = open(output_file, "w")

   # Loop through thrsh array and write out location where threshold is met

   for n in range(td):
      for i in range(i1, il+1):
          for j in range(j1, jl+1):
               if(thrsh[i][j][n] == 1.0):
                  print("i = " + str(i) + " j = " + str(j) + " n = " + str(n))
                  st_u,st_v = steer_vel(st_u,st_v,i,j,n)
                  # Check for unrealistic steering velocities caused by high topography in climate data
                  if( (st_u < 1.e5) and (st_u > -1.e5) and (st_v < 1.e5) and (st_v > -1.e5) ):
                      speed = np.sqrt( (uwnd[i][j][k850][n])**2 + (vwnd[i][j][k850][n])**2 )

                      file.write(str(date)+" "
                           +str(time)+" "
                           +str(lat[j])+" "+str(lon[i])+" "
                           +str(st_u)+" "
                           +str(st_v)+" "
                           +str(speed)+" "
                           +str(topo[i][j])+" "
                           +str(fld1[i][j][n])+" "
                           +str(fld2[i][j][n])+" "
                           +str(fld3[i][j][n])+" "
                           +str(fld4[i][j][n]) + " "
                           +str(fld5[i][j][n]) +" "
                           +str(fld6[i][j][n])+"\n")
   file.close()


# The steering velocity is approximated by the average wind at 700 hPa
# over box of size box_size
#
#  st_u and st_v = Zonal and meridional components of steering velocity
#  ii,jj = Grid positions where steering velocity is to be calculated
#  nn =  Time level when steering velocity is to be calculated
def steer_vel(st_u,st_v,ii,jj,nn):
    box_size = 4.1  #Size of averaging box (degrees)

    #Distance between neighbouring grid points (degrees)
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]


    #Number of grid points surrounding ii,jj in averaging box
    ibox = int(box_size * 0.5 / dlon) + 1
    jbox = int(box_size * 0.5 / dlat) + 1


    #The number of grid points plus/minus to be included in averaging
    ip = ii + ibox
    im = ii - ibox
    jp = jj + jbox
    jm = jj - jbox
    if(ip > xd):
        ip = xd

    if(im < 0):
        im = 0

    if(jp > yd):
        jp = yd


    if(jm < 0):
        jm = 0

    count = 0   #Counter to determine how many grid points used in wind averaging
    st_u  = 0.0
    st_v  = 0.0


    # Calculate average winds in box
    for i in range(im, ip+1):
        for j in range(jm, jp+1):
            if(i != xd):
                if(j != yd):
                    count = count + 1
                    st_u = st_u + uwnd[i][j][k700][nn]
                    st_v = st_v + vwnd[i][j][k700][nn]




    st_u = st_u/float(count)
    st_v = st_v/float(count)

    print("st_u = " + str(st_u) + " --- st_v = " + str(st_v) )


    return st_u,st_v


# Put threshold fields into 3-d arrays for efficient array manipulation
def fill_fields():
    global fild1
    global fild2
    global fild3
    global fild4
    global fild5
    global fild6

    print("Inside fill_fields")

    #x!  Fill fields

    smth = np.zeros([xd, yd],float)
    for i in range(xd):
        for j in range(yd):
            for n in range(td):
                fld1[i][j][n] = owz[i][j][k850][n]
                fld2[i][j][n] = owz[i][j][k500][n]
                fld3[i][j][n] = rh[i][j][k950][n]
                fld4[i][j][n] = rh[i][j][k700][n]
                fld5[i][j][n] = wsh[i][j][k200][n]
                fld6[i][j][n] = mrsh[i][j][k950][n]



    #If large values of OWZ are encountered the format specification of f5.1 means "*****" is written
    # to the output file.  To avoid this set a maximum value of 999.9 for OWZ.   12-05-2016

    np.where(fld1 > 999.9, 999.9,fld1)
    np.where(fld2 > 999.9, 999.9, fld2)

    print("End of fill_fields")



# By first setting the combined threshold array (thrsh) to everywhere
# represent a satisfied threshold (thrsh = 1.0), the unsatisfied locations
# for each field can be progressively set to 0.0, to effectively mask all
# unsatisfied grid points.  This leaves only those points with all threshold
# variables satisfied, with a value of 1.0.
# idcbm = Threshold combination ID
def find_threshold(idcmb):
   print "** finding threshold ** "+date+" - "+time

   global thrsh

   RHmax=120.0     # Maximum RH limit to ensure no highly unrealistic values appear
   RHmin=-10.0     # Minimum RH limit to ensure no highly unrealistic values appear
   SDmax=1.0       # Maximum SD limit to ensure no highly unrealistic values appear
   SDmin=-10.0     # Minimum SD limit to ensure no highly unrealistic values appear

   fldmax = 0      # Maximum 950 hPa RH or SD limit
   fldmin = 0      # Minimum 950 hPa RH or SD limit

   print("Inside find_threshold")

   # Reset thrsh array to 1.0
   thrsh = np.ones([xd,yd,td],float)

   print(thrsh)


# Search for threshold failure locations and set thrsh values to zero
   thrsh[fld1 < Towz850[idcmb]] = 0.0
   thrsh[fld2 < Towz500[idcmb]] = 0.0
   thrsh[fld3 < Trh950[idcmb]] = 0.0
   thrsh[fld4 <  Trh700[idcmb]] = 0.0
   thrsh[fld5 >  Twsh[idcmb]] = 0.0
   thrsh[fld6 <  Tsh950[idcmb]] = 0.0


# Search for unreasonable fld values and set thrsh values to zero
# Climate data have missing values where fields intersect topography.
# The windowing program creates dodgy numbers when encountering these fields.

# Set fld3max etc.
   if(Lsd):
     fldmax = SDmax
     fldmin = SDmin
   else:
     fldmax = RHmax
     fldmin = RHmin


   thrsh[fld1 > 1e5] = 0.0
   thrsh[fld2 > 1e5] = 0.0
   thrsh[fld3 > fldmax] = 0.0
   thrsh[fld4 > fldmax] = 0.0
   thrsh[fld5 > 1e5] = 0.0
   thrsh[fld6 > 120] = 0.0



   thrsh[fld1 < -1e5] = 0.0
   thrsh[fld2 < -1e5] = 0.0
   thrsh[fld3 < fldmin] = 0.0
   thrsh[fld4 < fldmin] = 0.0
   thrsh[fld5 < -1e5] = 0.0
   thrsh[fld6 < -10.0] = 0.0

   print("End of find_threshold")



#
# Calculate OWZ
#
def owz_calc():
   print "** Inside owz_calc ** "+date+" - "+time

   global owz
   global avor

   DVDX = 0.0        #
   DUDY = 0.0        # Wind gradients
   DVDY = 0.0        #
   DUDX = 0.0        #

   XX  = 0.0                    # Normalized OW

   scale = 1.0e6                # Scale parameter for OWZ

   #Reciprocal of two times the grid separation distance (m)
   idx = 0.0
   idy = 0.0

   #Calculate corioils limit equivalent to a latitude of 20 degrees
   Cormax = 2.0 * 7.292116e-5 * np.sin(20.0 * 0.01745329) #Maximum contribution of the coriolis paramter to OWZ

   #dlat/dlon = Longitude/Latitude increments
   dlon = lon[1] - lon[0]
   dlat = lat[1] - lat[0]
   idy = 1.0 / (2.0 * 111.2 * 1000.0 * dlat)


   for j in range(yd):
       Cor = 2.0 * 7.292116e-5 * np.sin(lat[j] * 0.01745329)
       if(lat[j] < 0.0 ):
           SgnCor =  -1.0
           if (Cor < -Cormax):
               Cor = -Cormax
       else:
           if (Cor > Cormax):
               Cor = Cormax
           SgnCor = 1.0

       idx = 1.0 / (2.0 * 111.2 * 1000.0 * dlon * np.cos(lat[j] * 0.01745329))


       for n in range(td):
          for k in range(zd):
              for i in range(xd):

                 if(j == 0):
                     DUDY = 2.0 * idy * (uwnd[i][j+1][k][n] - uwnd[i][j][k][n])
                     DVDY = 2.0 * idy *(vwnd[i][j+1][k][n] - vwnd[i][j][k][n] )
                 elif(j+1 == yd):
                     DUDY = 2.0 * idy * (uwnd[i][j][k][n] - uwnd[i][j-1][k][n])
                     DVDY = 2.0 * idy * (vwnd[i][j][k][n] - vwnd[i][j-1][k][n])
                 else:
                     DVDY = idy * (vwnd[i][j+1][k][n] - vwnd[i][j-1][k][n])
                     DUDY = idy * (uwnd[i][j+1][k][n] - uwnd[i][j-1][k][n])


                 if(i == 0 ):
                     DVDX = 2.0 * idx * ( vwnd[i+1][j][k][n] - vwnd[i][j][k][n])
                     DUDX = 2.0 * idx * (uwnd[i+1][j][k][n] - uwnd[i][j][k][n])

                 elif( i + 1 == xd ):
                     DVDX = 2.0 * idx * (vwnd[i][j][k][n] - vwnd[i-1][j][k][n])
                     DUDX = 2.0 * idx * (uwnd[i][j][k][n] - uwnd[i-1][j][k][n])

                 else:
                     DVDX = idx * (vwnd[i+1][j][k][n] - vwnd[i-1][j][k][n])
                     DUDX = idx * (uwnd[i+1][j][k][n] - uwnd[i-1][j][k][n])




                 ZZ  = (DVDX - DUDY)**2  # Relativevorticity
                 EE  = (DUDX - DVDY)**2  # EE & EF = Shear and strain defrmation
                 FF  = (DVDX + DUDY)**2  #
                 ZZA = DVDX - DUDY + Cor #absolute vorticity

                 XX  = (ZZ - EE - FF)/(ZZ+1.0e-20)

#                  XX  = np.amax(XX,0)
                 if(XX < 0):
                    XX = 0.0

                 avor[i][j][k][n] = ZZA*scale
                 owz[i][j][k][n] = XX * SgnCor * ZZA * scale

   print "end of owz_calc"



# Calculate the magnitude of the wind shear at every level relative to the 850 hPa level.
# Sub-tropical jet (STJ) calculation
def wsh_calc():
    print "** Inside wsh_calc ** "+date+" - "+time

    global uwnd
    global vwnd
    global wsh

    smth = np.zeros([xd, yd],float) # 2-d array to be passed to smoothing subroutine

    u850 = 0
    v850 = 0                # wind components at 850 hPa

    yr = 0
    md = 0                  # Year number, monthday number
    store = 0               # Temporary variable for calculating yr and md
    avlat = 0.0             # Average of all latitudes in lat array (used to find hemisphere)
    latjet = 0.0            # Latitude of the STJ
    speed = np.zeros(yd)    # Wind speed at 200 hPa
    shear = 0.0             # meridional shear of the 200 hPa wind


    print("num_smth = " + str(num_smth))

    #Smooth uwnd
    if (num_smth > 0):
        for n in range(td):
            for k in range(zd):
                for i in range(xd):
                     for j in range(yd):
                         smth[i][j] = uwnd[i][j][k][n]

                smth = smoother(smth,xd,yd,num_smth)

                for i in range(xd):
                    for j in range(yd):
                        uwnd[i][j][k][n] = smth[i][j]

    #! Smooth vwnd
    print("**** Smooth vwnd *****")
    if (num_smth > 0):
        print("Smooth vwnd")
        for n in range(td):
            for k in range(zd):
                for i in range(xd):
                     for j in range(yd):
                         smth[i][j] = vwnd[i][j][k][n]
                smth = smoother(smth,xd,yd,num_smth)

                for i in range(xd):
                    for j in range(yd):
                        vwnd[i][j][k][n] = smth[i][j]

    #! Calculate wind shear
    for n in range(td):
        for i in range(xd):
            for j in range(yd):
                u850 = uwnd[i][j][k850][n]
                v850 = vwnd[i][j][k850][n]
                for k in range(zd):
                    wsh[i][j][k][n] = np.sqrt( (uwnd[i][j][k][n]-u850)**2 + (vwnd[i][j][k][n]-v850)**2 )



    # Determine the position of the sub-tropical jet (STJ) if it exists or the position
    # of the closest jet to the equator (more than 7 degrees), that is stronger than
    # 25 m/s, and the zonal component is greater than 15 m/s.

    # Calculate yr and md
    store = int(date)/10000
    yr = int(store)
    md = int(date) - yr*10000
    if(time[0] < 24):
         time[0] = time[0]*100


    # Determine the hemisphere (assume the sign of the average latitude determines the hemisphere)
    avlat =  np.sum(lat)/yd
    print("sum of lat = " + str(avlat))


    #! Open output file to write STJ data
    file = open("STJ_"+str(yr)+str(md)+str(time[0])+"_"+Hem+"_00.txt", "w")

    for n in range(td):
       for i in range(xd):
           for j in range(yd):
               speed[j] = np.sqrt( (uwnd[i][j][k200][n])**2 + (vwnd[i][j][k200][n])**2 )

           if (Hem == "NH"):
               for j in range(yd-1):
                   latjet = lat[yd-1]
                   if ( (lat[j] > 7.0) and (uwnd[i][j][k200][n] > 15.0) ) :
                       shear = speed[j+1]-speed[j-1]
                       if ( (shear < 0.0) and (speed[j] > 25.0) ):
                           latjet = lat[j]
           elif (Hem == "SH"):
               for j in range(yd-1,0, -1):
                   latjet = lat[0]
                   if ( (lat[j] < -7.0) and (uwnd[i][j][k200][n] > 15.0) ) :
                       shear = speed[j-1]-speed[j+1]
                       if ( (shear < 0.0) and (speed[j] > 25.0) ):
                           latjet = lat[j]
           file.write(str(lon[i])+" , "+str(latjet)+" "+str(yr)+str(md)+str(time[0])+"\n")
    file.close()


# This subroutine smooths by averaging values using a 2-dimensional 2-4-2 weighting, plus 1
# for the diagonals.  That is, the grid points to the north, south, east and west have a weighting
# of 2, the central grid point has a weigting of 4 and the NW, NE, SE, SW grid points have a
# weighting of 1.
#
# smth     = Field to be smoothed
# xd yd    = X and Y dimensions
# num_smth = Number of smoothing operations
def smoother(smth,xd,yd,num_smth):
    print "** inside smoother ** "+date+" - "+time

    tmp = np.zeros([xd, yd],float) # 2-d array to be passed to smoothing subroutine

    #Fractions used in averaging
    third   = 1.0/3.0
    sixth   = 1.0/6.0
    ninth   = 1.0/9.0
    twelfth = 1.0/12.0

    for n in range(num_smth):
        tmp = np.zeros([xd, yd],float) # 2-d array to be passed to smoothing subroutine

        # Interior
        for i in range(1,xd-1):
           for j in range(1,yd-1):
                 tmp[i][j] =   (0.25*smth[i][j]
                 + 0.125 * ( smth[i-1][j] + smth[i+1][j] + smth[i][j-1] + smth[i][j+1] )
                 + 0.0625 * ( smth[i-1][j-1]+smth[i-1][j+1]+smth[i+1][j+1]+smth[i+1][j-1] ))

        #! Edges
        for i in range(1,xd-1):
           tmp[i,0]   = (third  * smth[i][0]
           + sixth *  ( smth[i-1][0]+smth[i+1][0]+smth[i][1] )
           + twelfth * ( smth[i-1][1]+smth[i+1][1] ))

           tmp[i][yd-1] = (third * smth[i][yd-1]
           + sixth *  ( smth[i-1][yd-1]+smth[i+1][yd-1]+smth[i][yd-2] )
           + twelfth * ( smth[i-1][yd-2]+smth[i+1][yd-2] ))



        for j in range(1,yd-1):
           tmp[0][j] =   (third * smth[0][j]
           + sixth *  ( smth[0][j-1] + smth[0][j+1]+smth[1][j] )
           + twelfth * ( smth[1][j-1]+smth[1][j+1] ))


           tmp[xd-1][j] = (third * smth[xd-1][j]
           + sixth *  ( smth[xd-1][j-1]+smth[xd-1][j+1]+smth[xd-2][j] )
           + twelfth * ( smth[xd-2][j-1]+smth[xd-2][j+1] ))


        #! Corners
        tmp[0][0] = (4*ninth*smth[0][0] + 2*ninth*( smth[0][1]+smth[1][0] ) + ninth*smth[1][1])
        tmp[0][yd-1]= (4*ninth*smth[0][yd-1] + 2*ninth*( smth[0][yd-2]+smth[1][yd-1] ) +ninth*smth[1][yd-2])
        tmp[xd-1][0]= (4*ninth*smth[xd-1][0]+ 2*ninth*( smth[xd-2][0]+smth[xd-1][1] ) +ninth*smth[xd-2][1])
        tmp[xd-1][yd-1]= (4*ninth*smth[xd-1][yd-1]+2*ninth*( smth[xd-2][yd-1]+smth[xd-1][yd-2] ) +ninth*smth[xd-2][yd-2])


        smth = tmp

    return smth




#    This subroutine ensures the relative humidity (RH) and specific humidity (SH) arrays are filled.
#    The input moisture variable can be RH, SH or mixing ratio (MR).  The logical flag corresponding to
#    the input moisture variable will be "true" and the other two logical flags "false".
#
#    Case:
#    Lrh = true --> SH is calculated and returned in "mrsh" array.
#    Lmr = true --> SH is calculated and returned in "mrsh", and RH is calculated and returned in "rh" array.
#    Lsh = true --> RH is calculated and returned in "rh" array.
def rh_calc(Lrh,Lmr,Lsh):
   print "** inside rh_calc ** "+date+" - "+time

   global rh
   global mrsh

   vps1=6.112    #! Saturation vapour pressure constants
   vps2=17.67
   vps3=243.5
   vps4=273.15
   vp1=0.622     #! Vapour pressure constant = Rd/Rv

   fac = 1.0     #! Factor that distinguishes vapour pressure as a function of
                 #! mixing ratio or specific humidity


   vpsat = np.zeros([xd,yd,zd,td],float) #! Saturation vapour pressure



    #! Set mr vs. sh factor
   if (Lmr):
     fac = 1.0
   else:
     fac = 1.0 - vp1

   # ! Calculate saturation vapour pressure
   print("Inside rh_calc - Lrh, Lmr, Lsh",Lrh,Lmr,Lsh)
   for n in range(td):
       for k in range(zd):
           for i in range(xd):
               for j in range(yd):
                   TdegC = temp[i][j][k][n] - vps4
                   vpsat[i][j][k][n] = vps1 * np.exp( (vps2*TdegC)/(TdegC+vps3) )
                   #if(TdegC > -20.0) write(6,*) "TdegC = ",TdegC,i,j,k

   #! If rel_hum exists calculate spec_hum
   if (Lrh):
      print("In LRH")
      for n in range(td):
          for k in range(zd):
              for i in range(xd):
                  for j in range(yd):
                     vp = rh[i][j][k][n] * vpsat[i][j][k][n]/100.0
                     mrsh[i][j][k][n]= vp1*vp/(lvl[k] - vp*fac)


   else:
      for n in range(td):
          for k in range(zd):
              for i in range(xd):
                  for j in range(yd):
                     vp = mrsh[i][j][k][n] * lvl[k] / ( vp1 + mrsh[i][j][k][n] * fac )
                     rh[i][j][k][n] = 100.0 * vp / vpsat[i][j][k][n]
                     #!             if(rh(i,j,k,n) > 5.0) write(6,*) "RH ",rh(i,j,k,n),i,j,k


   #! Calculate saturation deficit if required and put in rh array (values should be negative)
   if (Lsd):
       print("in LSD")
       for n in range(td):
          for k in range(zd):
              for i in range(xd):
                   for j in range(yd):
                       qsat = vp1*vpsat[i][j][k][n]/(lvl[k] - vpsat[i][j][k][n]*fac)
                       rh[i][j][k][n] = mrsh[i][j][k][n] - qsat

       rh = np.multiply(rh, 1000.0)   #! Change to g/kg



   #! Replace mixing ratio with specific humidity if required
   if (Lmr):
       mrsh = mrsh/(1 + mrsh)

   #! Change units of mrsh to g/kg
   mrsh = np.multiply(mrsh, 1000.0)


#
#   Set Date/time of this data
#
def setDateTime(d,t):
    global date, time
    date = d
    time = t


# Open and read the NetCDF file
def read_data(nc_file):
    print "** reading data ** "+date+" - "+time
    global yd, xd,zd,td
    global Lowz,Lav,Lrh,Lwsh,Lmr,Lsh,Ltop,Lsfg
    global owz,rh,wsh,thrsh,uwnd,vwnd,temp,fld1,fld2,fld3,fld4,fld5,fld6,topo,topo3d,avor,mrsh
    global lat, lon,lvl#, date, time
    global i1, il, j1, jl
    global k950,k850,k700,k500,k200

    i = 0
    j = 0
    k = 0
    max_lvl = 0

    # Set logical variables to the default value of TRUE.  They are set to false when found to not exist
    Lowz = False
    Lrh  = True
    Lwsh = False
    Lmr  = False
    Lsh  = True
    Ltop = True
    Lsfg = True

    nc_fid = Dataset(nc_file, 'r')  # Dataset is the class behavior to open the file
                             # and create an instance of the ncCDF4 class

    print("*************************")
    # print(nc_dims)
    #ncdump(nc_fid)

    print("Getting dimension ID's")

    if("lat" not in nc_fid.dimensions):
        print("Error: lat dimension not available")
        return -1

    latid = nc_fid.dimensions["lat"]

    if("lon" not in nc_fid.dimensions):
        print("Error: lon dimension not available")
        return -1

    lonid = nc_fid.dimensions["lon"]


    if("time" not in nc_fid.dimensions):
        print("Error: time dimension not available")
        return -1

    timid = nc_fid.dimensions["time"]


    if("lvl" not in nc_fid.dimensions):
        print("Error: lvl dimension not available")
        return -1

    lvlid = nc_fid.dimensions["lvl"]


    print("Getting dimension lengths")

    yd = latid.size
    xd = lonid.size
    zd = lvlid.size
    td = timid.size



    print("Getting variable ID's")

    if("rh" not in nc_fid.variables):
        print("rh variable not available")
        Lrh = False
        if("mrsh" not in nc_fid.variables):
            print("spec_hum variable not available")
            Lmr = False
            print("No recognised moisture parameter available.")
            print("Program can use relative humidity(%),")
            print("specific humidity (kg/kg) and mixing ratio (kg/kg).")
            print("Check variable names match those in the data file.")
            return -1
        else:
            mrshid = nc_fid.variables["mrsh"]
            Lsh = False
    else:
        Lmr = False
        Lsh = False


    rhid = nc_fid.variables["rh"]

    if("mrsh" not in nc_fid.variables):
        print("Error: mrsh variable not available")
        return -1
    mrshid = nc_fid.variables["mrsh"]


    if("uwnd" not in nc_fid.variables):
        print("Error: uwnd variable not available")
        return -1
    uid = nc_fid.variables["uwnd"]


    if("vwnd" not in nc_fid.variables):
        print("Error: vwnd variable not available")
        return -1
    vid = nc_fid.variables["vwnd"]


    if("temp" not in nc_fid.variables):
        print("Error: temp variable not available")
        return -1
    tdid = nc_fid.variables["temp"]

    if("topog" not in nc_fid.variables):
        Ltop = False
        if("topo" not in nc_fid.variables):
            print("Error: topog variable not available")
            Ltop = False
            print("No recognised topography field available.")
            print("Program can use 'topog' or 'sfc_geop'.")
            return -1
        else:
            topid = nc_fid.variables["topo"]
    else:
        Lsfg = False
        topid = nc_fid.variables["topog"]




    if("lat" not in nc_fid.variables):
        print("Error: lat variable not available")
        return -1

    if("lon" not in nc_fid.variables):
        print("Error: lon variable not available")
        return -1


    if("lvl" not in nc_fid.variables):
        print("Error: lvl variable not available")
        return -1


    print("Allocating space for field variables")
    owz = np.zeros([xd,yd,zd,td],float)
    rh = np.zeros([xd,yd,zd,td],float)
    wsh = np.zeros([xd,yd,zd,td],float)
    thrsh = np.zeros([xd,yd,td],float)
    uwnd = np.zeros([xd,yd,zd,td],float)
    vwnd = np.zeros([xd,yd,zd,td],float)
    temp = np.zeros([xd,yd,zd,td],float)
    fld1 = np.zeros([xd,yd,td],float)
    fld2 = np.zeros([xd,yd,td],float)
    fld3 = np.zeros([xd,yd,td],float)
    fld4 = np.zeros([xd,yd,td],float)
    fld5 = np.zeros([xd,yd,td],float)
    fld6 = np.zeros([xd,yd,td],float)
    topo = np.zeros([xd,yd],float)
    avor = np.zeros([xd,yd,zd,td],float)
    mrsh = np.zeros([xd,yd,zd,td],float)



    print("Allocating space for dimension variables")
    lat  = np.zeros([yd],float)
    lon  = np.zeros([xd],float)
    lvl  = np.zeros([zd],float)


    print("Filling dimension arrays")
    #! Fill the dimension arrays (note if a subset of the array is required use nf_get_vara_real)

    #for i in nc_fid.variables:
    #    print(i, nc_fid.variables[i][:], nc_fid.variables[i].shape)

    lat = nc_fid.variables["lat"][:]
    lon = nc_fid.variables["lon"][:]
    lvl = nc_fid.variables["lvl"][:]
    #date = nc_fid.variables["valid_date"][:]
    #time = nc_fid.variables["valid_time"][:]

    print(date)
    print(time)

    #! Convert lvl units from Pa to hPa if necessary
    max_lvl = 0.0
    for j in range(zd):
         if (lvl[j] > max_lvl):
             max_lvl = lvl[j]

    if (max_lvl > 1100.0):
        lvl = lvl * 0.01

    print("Filling 4-d variable arrays")
   # if (Lowz):
   #     owz = nc_fid.variables["n_ow_zta"][:,:,:,:]
    _rh = np.zeros([xd,yd,zd,td],float)
    if (Lrh):
        _rh =  rhid[:]#:,:,:,:]


    _mrsh = np.zeros([xd,yd,zd,td],float)
    if (Lmr or Lsh):
       _mrsh = mrshid[:]
    _uwnd = uid[:]


    _vwnd = vid[:]#,:,:,:]

    _temp = tdid[:]#,:,:,:]

    for i in range(xd):
        for j in range(yd):
            for k in range(zd):
                for n in range(td):
                    rh[i][j][k][n] = _rh[n][k][j][i]
                    uwnd[i][j][k][n] = _uwnd[n][k][j][i]
                    vwnd[i][j][k][n] = _vwnd[n][k][j][i]
                    temp[i][j][k][n] = _temp[n][k][j][i]
                    if (Lmr or Lsh):
                        mrsh[i][j][k][n] = _mrsh[n][k][j][i]


    print("xd " + str(xd))
    print("yd " + str(yd))
    print("zd " + str(zd))
    print("td " + str(td))

    #! Find standard level integers
    for k in range(zd):
         if(np.abs(lvl[k]-950.) < 0.001):
             k950=k
         if(np.abs(lvl[k]-850.) < 0.001):
             k850=k
         if(np.abs(lvl[k]-700.) < 0.001):
             k700=k
         if(np.abs(lvl[k]-500.) < 0.001):
             k500=k
         if(np.abs(lvl[k]-200.) < 0.001):
             k200=k

    print('k950 = ',k950, lvl[k950])
    print('k850 = ',k850, lvl[k850])
    print('k700 = ',k700, lvl[k700])
    print('k500 = ',k500, lvl[k500])
    print('k200 = ',k200, lvl[k200])


    #! If the 950 hPa pressure level does not exist, look for 925 hPa (climate data)
    if (k950 == -1):
         print("950 hPa pressure level does not exist in this data file.")
         print("Replace with 925 hPa level.")
         for k in range(zd):
            if(np.abs(lvl[k]-925.) < 0.001):
               k950=k
         if (k950 == 0):
           print("Could not find 925 hPa pressure level either.")
         else:
           print("Found 925 hPa pressure level. k950 = ",k950, lvl[k950])

    #! Fill missing arrays
    if(not Lowz):
        owz_calc()
    rh_calc(Lrh,Lmr,Lsh)

    if(not Lwsh):
        wsh_calc()

    print("Filling 3-d variable arrays")

    print("Lsfg" + str(Lsfg))
    if(Lsfg):
        topo3d = topid[:]
        print(topo3d.shape)

        for i in range(xd):
            for j in range(yd):
                topo[i][j] = topo3d[j][i]



    print("Filling 2-d variable arrays")
    #! Fill the 2-d variable arrays
    if(Ltop):
        print("Ltop")
        topoTemp = topid[:,:]
        for i in range(xd):
            for j in range(yd):
                topo[i][j] = topoTemp[j][i]
      # if(status /= nf_noerr) call error_handle(status)

    print("Topo check: topo(20,20) ="+ str(topo[19][19]))
    #print("Extended topo check")
    #for i in range(xd):
    #    for j in range(yd):
    #       print(str(topo[i][j])+ " " + str(lon[i]) + " " + str(lat[j]) +" "+ str(i) + " " +str(j))


#! Close the NetCDF file
    nc_fid.close()

#   if(status /= nf_noerr) call error_handle(status)

#! Get first and last output window indices
    i1 = 0
    il = xd -1
    for i in range(xd) :
      if (np.abs(lon[i]-Wlnmin) < 0.001): i1 = i
      if (np.abs(lon[i]-Wlnmax) < 0.001): il = i

    j1 = 0
    jl = yd -1
    print("Wltmin " + str(Wltmin))
    print("Wltmax " + str(Wltmax))
    for j in range(yd):
      #print(str(lat[j]) + " - " + str(Wltmin) + "=" + str(np.abs(lat[j]-Wltmin)))
      if (np.abs(lat[j]-Wltmin) < 0.001): j1 = j
      if (np.abs(lat[j]-Wltmax) < 0.001): jl = j
      #print(str("j1 = " + str(j1) + " - jl = " + str(jl)))



#
# Module to print ncdump
#
def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print ("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print ('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print ("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print ("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print ('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print ("NetCDF dimension information:")
        for dim in nc_dims:
            print ("\tName:", dim)
            print ("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print ("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print ('\tName:', var)
                print ("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print ("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars



#
#   Main Entry point of threshold module
#
def process(inputFile,date,time,hemisphere):
    global Hem
    Hem = hemisphere

    if(Hem == "SH"):
        input_read("input_thresh_SH")
    elif(Hem == "NH"):
        input_read("input_thresh_SH")
    else:
        print("Invaid Hemisphere")
        return

    setDateTime(date,time)

#    print(sys.argv[1])
    read_data(inputFile)

    print("xd,yd,zd,td",xd,yd,zd,td)


   # ! Check arrays
    print("owz - 850 " + str(owz[1][1][k850][0]))
    print('owz - 500' + str(owz[1][1][k500][0]))
    print( 'rh  - 950 '+ str(rh[1][1][k950][0]))
    print( 'rh  - 700 ' + str(rh[1][1][k700][0]))
    print( 'wsh - 850-200 '+ str(wsh[1][1][k200][0]))
    print( 'sh  - 950 ' + str(mrsh[1][1][k950][0]))

    fill_fields()


    for n in range(num_th_comb):
        find_threshold(n)
        if (Loutput):
            output(n)
#        if (Lwrite and (n == 1)):
#            write_data()


#if __name__ == "__main__":
#    main()
