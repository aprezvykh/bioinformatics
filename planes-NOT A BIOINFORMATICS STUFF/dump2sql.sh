#!/bin/bash
#Shell script to listen parse and aggregate ADSB data
#Run in foreground to see progress bar
#run  in background with: dump2sql.sh >/dev/null &

Timestamp=$(date +"%Y-%m-%d")   #separate logfile every day
ADSBhost="127.0.0.1"
ADSBport="30003"
ADSBlog="/home/alexander/adsb-log-$Timestamp.txt"
#now the mySQL credentials...use certificates in next version
SQLsrv='XXXXXX.YYYYYY.eu-central-1.rds.amazonaws.com'
SQLusr="ZZZZZ"
SQLpwd="*****"
SQLdbs="DDDDDD"
SQLqueries="/home/alexander/dump.sql"
counter=5      #for loop control
countmax=60    #script can stop after so many loops

echo "DUMP1090 Aggregator2mySQL by Matthias Gemelli 2016"
echo "listening to server: $ADSBhost writing to log $ADSBlog"
echo "mySQL as $SQLusr to $SQLsrv"
echo "progress bar: . for every message, X for every location, Y for missing cal                                                                                                                lsign"
echo "exit with Ctrl-C or set message limit in countmax"
echo "-------"
echo "MsgCount,HexIdent,Date,Time,Lat,Long,Callsign,Altitude,Speed,Track,Vertica                                                                                                                l" >> "$ADSBlog"

#now declare the arrays
declare -a arr_call
declare -a arr_alti
declare -a arr_sped
declare -a arr_trck
declare -a arr_vert

#loop through the Netcat data
nc -d $ADSBhost $ADSBport | while IFS="," read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10                                                                                                                 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22
do     #loop until a break is thrown

#first update the filename
Timestamp=$(date +"%Y-%m-%d")   #separate logfile every day
ADSBlog="/home/alexander/adsb-log-$Timestamp.txt"

#now read the relevant data fields in every ADSB record
#echo "Field 05 HexIdent         :$f2"
#echo "Field 07 Date message gen :$f7"
#echo "Field 08 Time message gen :$f8"
#echo "Field 11 Callsign         :$f11"
#echo "Field 12 Altitude         :$f12"
#echo "Field 13 GroundSpeed      :$f13"
#echo "Field 14 Track            :$f14"
#echo "Field 15 Latitude         :$f15"
#echo "Field 16 Longitude        :$f16"
#echo "Field 17 Vertical Rate    :$f17"

#now save the data into array, using HexIdent as index
ident=$((0x${f5}))
if [ "$f11" != "" ];  then arr_call[ident]="$f11"; fi
if [ "$f12" != "" ];  then arr_alti[ident]="$f12"; fi
if [ "$f13" != "" ];  then arr_velo[ident]="$f13"; fi
if [ "$f14" != "" ];  then arr_trck[ident]="$f14"; fi
if [ "$f17" != "" ];  then arr_vert[ident]="$f17"; fi

#if position and if callsign is broadcast
if [ "$f15" != "" ]; then  #if f15 not empty
if [ "${arr_call[ident]}" != "" ]; then #if callsign is already known

punkt="$counter,$f5,$f7,$f8,$f15,$f16,${arr_call[ident]},${arr_alti[ident]}"
punkt="$punkt,${arr_velo[ident]},${arr_trck[ident]},${arr_vert[ident]}"
echo "$punkt" >> "$ADSBlog"      #write to log file

#now compose SQL statement
QUERY="INSERT INTO esp8266data.flights "
QUERY="$QUERY(msgcount,hexident,date,time,lat,lon,sign,alti,speed,trck,vert) VAL                                                                                                                UES"
QUERY="$QUERY ($counter,\"$f5\",\"$f7\",\"$f8\",\"$f15\",\"$f16\","
QUERY="$QUERY \"${arr_call[ident]}\",${arr_alti[ident]},${arr_velo[ident]},"
QUERY="$QUERY ${arr_trck[ident]},${arr_vert[ident]});"
echo "$QUERY" >$SQLqueries  #write SQL to a file before executing
#echo "$QUERY"
mysql -h $SQLsrv -u $SQLusr -p$SQLpwd $SQLdbs <$SQLqueries

#progress bar on shell - X for position, Y for pos without callsign
printf "X"
else  #if no callsign is known at position
printf "Y $f5"
fi #if callsign is aleady known
else  #what to do if no position is given
printf "."  #progress bar
fi #if f15 not empty


#reset the array if it is midnight (fewer planes)
#if reached max counter then break from loop
#if [ "$counter" -gt "$countmax" ]; then break; fi
((counter++))   #increase counter

done            #netcat listener loop
#attention: variables set within the loop stay in the loop
echo "done"