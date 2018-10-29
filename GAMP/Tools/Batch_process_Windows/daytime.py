#!/usr/bin/env python
# coding:utf-8

################################################################################
# Import Python modules
import datetime, time


################################################################################
# FUNCTION: 
################################################################################
def date_diff_seconds(date):
	date_1970jan1 = datetime.datetime(1970, 1, 1, 0, 0, 0)
	time_delta    = date - date_1970jan1
	return time_delta.days * 24 * 3600 + time_delta.seconds


################################################################################
# FUNCTION: to convert 2-digit year to 4-digit year
################################################################################
def yy2yyyy(yy):
	if yy >= 1900:
		yyyy = yy
	elif yy > 50 and yy < 1900:
		yyyy = yy + 1900
	elif yy <= 50:
		yyyy = yy + 2000

	return yyyy


################################################################################
# FUNCTION: to convert year, day of year to GPS week, day of week
################################################################################
def yrdoy2gpst(year, doy):
	year0 = yy2yyyy(year)
	date_1980jan6 = datetime.datetime(1980, 1, 6, 0, 0, 0)
	date          = datetime.datetime(year0, 1, 1, 0, 0, 0)
	time_delta    = date - date_1980jan6
	days_delta    = time_delta.days + doy - 1
	gps_week      = days_delta / 7
	gps_dow       = days_delta - gps_week * 7

	return gps_week, gps_dow


################################################################################
# FUNCTION: to convert year, day of year to GPS week, day of week
################################################################################
def gpst2yrdoy(gps_week, gps_sow):
	mjd = gpst2mjd(gps_week, gps_sow)
	year, doy = mjd2yrdoy(mjd)

	return year, doy


################################################################################
# FUNCTION: to convert year, day of year to month, day
################################################################################
def yrdoy2ymd(year,doy):
	year0 = yy2yyyy(year)
	days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
	if ((year0 % 4 == 0) and ((year0 % 100 != 0) or (year0 % 400 == 0))):
		days_in_month[1] = 29
	id0 = doy
	for imonth in range(1,12):
		id0 = id0 - days_in_month[imonth-1]
		if (id0 > 0):
			continue
		iday = id0 + days_in_month[imonth-1]
		break

	return imonth, iday


################################################################################
# FUNCTION: to convert GPS week, seconds of week to modified Julian date (MJD)
################################################################################
def gpst2mjd(gps_week, gps_sow):
	mjd = gps_week * 7 + 44244 + gps_sow / 86400

	return mjd


################################################################################
# FUNCTION: to convert year, month, day or year, day of year to modified Julian date (MJD)
################################################################################
def ymd2mjd(year, month, day):
	doy_of_month = [0,31,59,90,120,151,181,212,243,273,304,334]

	year0 = yy2yyyy(year)
	# check the input data
	if (year0 < 0 or month < 0 or month > 12 or day > 366) or (month != 0 and day > 31):
		print ' *** ERROR (ymd2mjd): Incorrect date (year,month,day): ', year0, month, day
		return

	# doy to month, day
	if month == 0:
		im, id0 = yrdoy2ymd(year0,day)
	else:
		im = month
		id0 = day

	year1 = year0
	if im <= 2:
		year1 = year1 - 1

	mjd = 365 * year0 - 678941 + int(year1 / 4) - int(year1 / 100) + int(year1 / 400) + id0
	im = im - 1
	if im != -1:
		mjd = mjd + doy_of_month[im]

	return mjd


################################################################################
# FUNCTION: to convert modified Julian date (MJD) to year, day of year
################################################################################
def mjd2yrdoy(mjd):
	year = int((mjd + 678940) / 365)
	doy = mjd - ymd2mjd(year,1,1)
	while doy <= 0:
		year = year - 1
		doy = mjd - ymd2mjd(year,1,1) + 1

	return year, doy


################################################################################
# FUNCTION: to get the week day
################################################################################
def weekd_str(year, month, day):
	weekds = datetime.datetime(year,month,day).strftime("%w")
	weekd = int(weekds)
	if weekd == 0:
		weekd_str = 'Sunday'
	elif weekd == 1:
		weekd_str = 'Monday'
	elif weekd == 2:
		weekd_str = 'Tuesday'
	elif weekd == 3:
		weekd_str = 'Wednesday'
	elif weekd == 4:
		weekd_str = 'Thursday'
	elif weekd == 5:
		weekd_str = 'Friday'
	elif weekd == 6:
		weekd_str = 'Saturday'

	return weekd_str