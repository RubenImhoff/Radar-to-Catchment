#
# Radar to Catchment v3.3 for Windows
#
# Copyright Ruben Imhoff
#
# Contact:	Ruben Imhoff (ruben.imhoff@wur.nl)
#				Aart Overeem (aart.overeem@wur.nl)
#				
#######################################################################################

print("This script requires some Python packages in order to run. If one of these packages is not yet installed, it is recommended to install them prior to running the script. For most systems, this is a matter of typing the following text in a command window: pip install [name of package] ")

import h5py
import numpy as np
import shapefile
import urllib
import ftplib
from osgeo import gdal
from osgeo import gdal_array
from osgeo import ogr, osr
import os
import sys
import osgeo.gdal
import subprocess
from osgeo import gdal, gdalnumeric, ogr, osr
from PIL import Image, ImageDraw
from subprocess import call
import Tkinter
import tkFileDialog
from Tkinter import *
import inspect
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from numpy import linspace
from numpy import meshgrid
from osgeo import gdalconst
from PIL import Image
from netCDF4 import Dataset
import netCDF4
import csv
import glob
from datetime import datetime

startTime = datetime.now()

# Ask what the user wants to do
while True:
	root = Tk()
	FirstQues = IntVar()
	FirstQues.set(1)

	options = [
		("Get precipitation values for a pixel", 1),
		("Get precipitation values for an area / a catchment", 2)
	]

	def ShowChoice():
		print(FirstQues.get())
		root.quit()	
			
	Label(root, text="""What would you like to do? Click on your choice and close the window.""", justify = LEFT, padx = 20).pack()

	for txt, val in options:
		Radiobutton(root, text=txt, padx = 20, variable = FirstQues, command = ShowChoice, value=val,tristatevalue=0).pack(anchor=W)
	
	mainloop()
	FirstQues = FirstQues.get()
	
	break

if FirstQues == 1:
	
	def safe():
		print(e1.get())	
		DataQuestion = e1.get()		
		root2.quit()
				
	root2 = Tk()
	root2.title("Indicate which type of dataset will be supplied. Fill out the number of your choice and click on enter.")
	Label(root2, text = "1 NLRadarhdf5_NL25", justify = LEFT).grid(row=0)
	Label(root2, text = "2 NLRadarhdf5_NL21", justify = LEFT).grid(row=1)
	Label(root2, text = "3 DERadar_Radolan", justify = LEFT).grid(row=2)		
	Label(root2, text = "4 NL KDC data in NETCDF", justify = LEFT).grid(row=3)
	Label(root2, text = "5 OPERA - Radar Europe", justify = LEFT).grid(row=4)
	Label(root2, text = "6 NASA GPM IMERG_HQprecipitation level 3", justify = LEFT).grid(row=5)	
		
	e1 = Entry(root2)
	e1.grid(row=2, column=3)
		
	Button(root2, text = "Enter", command = safe).grid(row=3,column=3,sticky=W,pady=4)
						
	mainloop()							
				
	DataQuestion = e1.get()
	print(DataQuestion)		
	
	if DataQuestion == "1":
		########################################################
		# Open the hdf5-files for the Netherlands, in this case the NL_25 dataset.
		########################################################
		root = Tkinter.Tk()
		root.withdraw()
		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.h5'):
					list.append(os.path.join(root,filename))	
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################

		while True:
			def safe():
				print(e1.get())
				print(e2.get())	
				Col = e1.get()
				Row = e2.get()				
				master.quit()
				
			master = Tk()
			Label(master, text = "Column number of the target cell", justify = LEFT).grid(row=0)
			Label(master, text = "Row number of the target cell", justify = LEFT).grid(row=1)
		
			e1 = Entry(master)
			e2 = Entry(master)
			e1.grid(row=0, column=1)
			e2.grid(row=1, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			Col = e1.get()
			Row = e2.get()			
			break
			
		# Ask if the user wants to change reflectivities into mm.
		
		def safe():
			print(e1.get())
			CalcQues = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "Do you want to convert reflectivities [dBZ] into [mm] of rain? Note: only use this").grid(row=0)
		Label(master, text = "tool when the provided data has reflectivities in dBZ. Please give the number of your answer: [1] Yes, [2] No and press enter").grid(row=1)
		
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		CalcQues = e1.get()				
		
		if CalcQues == "1":
			def safe():
				print(e1.get())
				NumMin = e1.get()		
				master.quit()
				
			master = Tk()
			Label(master, text = "What is the temporal resolution of your data (e.g. 5 min. data)? Please give the amount in minutes.").grid(row=0)
		
			e1 = Entry(master)
			e1.grid(row=0, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			NumMin = e1.get()
			NumMin = float(NumMin)
			
		else:
			print("No conversion will take place.")	

		########################################################
		# It is time to start the 'For Loop'
		########################################################
		Outfile = os.path.join(dirpath,'Results','PixelPrec.csv') 				
		outcsv = open(Outfile,"w")
		writer = csv.writer(outcsv, delimiter = ",")		
		writer.writerow(['Name','Precipitation [mm]'])	        
                
		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/image1/image_data'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
		
			if CalcQues == "1":
				filedata = 0.5 * filedata - 32
				filedata = ((10.**((filedata-23)/16))/(60/NumMin))
			if CalcQues == "2":
				filedata = filedata * 0.01	 
	
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]		
			print(basenametxt)		
	
			# Multiply the data with 0.01 to get mm.
			filedata = filedata * 0.01	
			
			Value = filedata[Row,Col]
			writer.writerow([basenametxt, Value])
			
		print("End of this script.")

	if DataQuestion == "2":
		########################################################
		# Open the hdf5-files for the Netherlands, in this case the NL_21 dataset.
		########################################################
		root = Tkinter.Tk()
		root.withdraw()
		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.h5'):
					list.append(os.path.join(root,filename))
	
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		while True:
			def safe():
				print(e1.get())
				print(e2.get())	
				Col = e1.get()
				Row = e2.get()				
				master.quit()
				
			master = Tk()
			Label(master, text = "Column number of the target cell", justify = LEFT).grid(row=0)
			Label(master, text = "Row number of the target cell", justify = LEFT).grid(row=1)
		
			e1 = Entry(master)
			e2 = Entry(master)
			e1.grid(row=0, column=1)
			e2.grid(row=1, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			Col = e1.get()
			Row = e2.get()			
			break
			
		# Ask if the user wants to change reflectivities into mm.
		
		def safe():
			print(e1.get())
			CalcQues = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "Do you want to convert reflectivities [dBZ] into [mm] of rain? Note: only use this").grid(row=0)
		Label(master, text = "tool when the provided data has reflectivities in dBZ. Please give the number of your answer: [1] Yes, [2] No and press enter").grid(row=1)
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		CalcQues = e1.get()				
		
		if CalcQues == "1":
			def safe():
				print(e1.get())
				NumMin = e1.get()		
				master.quit()
				
			master = Tk()
			Label(master, text = "What is the temporal resolution of your data (e.g. 5 min. data)? Please give the amount in minutes.").grid(row=0)
		
			e1 = Entry(master)
			e1.grid(row=0, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			NumMin = e1.get()
			NumMin = float(NumMin)
			
		else:
			print("No conversion will take place.")	
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		Outfile = os.path.join(dirpath,'Results','PixelPrec.csv') 				
		outcsv = open(Outfile,"wb")
		writer = csv.writer(outcsv, delimiter = ",")		
		writer.writerow(['Name','Precipitation [mm]'])
		
		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/image1/image_data'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
	
			if CalcQues == "1":
				filedata = 0.5 * filedata - 32
				filedata = ((10.**((filedata-23)/16))/(60/NumMin))
			if CalcQues == "2":
				filedata = filedata * 0.01	
	
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]		
			print(basenametxt)		
	
			# Multiply the data with 0.01 to get mm.
			filedata = filedata * 0.01	
			
			Value = filedata[Row,Col]
			writer.writerow([basenametxt, Value])
			
		print("End of this script.")

	if DataQuestion == "3":
		########################################################
		# Open the radolan radar dataset for Germany.
		########################################################
		root = Tkinter.Tk()
		root.withdraw()
		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with ASCII-files and press OK')
		radolanlist = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.asc'):
					radolanlist.append(os.path.join(root,filename))
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		Outfile = os.path.join(dirpath,'Results','PixelPrec.csv') 				
		outcsv = open(Outfile,"wb")
		writer = csv.writer(outcsv, delimiter = ",")		
		writer.writerow(['Name','Precipitation [mm]'])
		
		while True:
			def safe():
				print(e1.get())
				print(e2.get())	
				Col = e1.get()
				Row = e2.get()				
				master.quit()
				
			master = Tk()
			Label(master, text = "Column number of the target cell", justify = LEFT).grid(row=0)
			Label(master, text = "Row number of the target cell", justify = LEFT).grid(row=1)
		
			e1 = Entry(master)
			e2 = Entry(master)
			e1.grid(row=0, column=1)
			e2.grid(row=1, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			Col = e1.get()
			Row = e2.get()			
			break
			
		for x in range(0, len(radolanlist)):
			radolanfilename = radolanlist[x]
			filedata = np.loadtxt(radolanlist[x], skiprows=6)
				
			basename = os.path.basename(radolanfilename)
			basenametxt = os.path.splitext(basename)[0]		
			print(basenametxt)		
			
			Value = filedata[Row,Col]
			writer.writerow([basenametxt, Value])
			
		print("End of this script.")

	if DataQuestion == "4":
		########################################################
		# Open the NETCDF-files for the Netherlands. In this case that is the monthly averaged precipitation, but it could also be hourly or five minute data.
		########################################################
		root = Tkinter.Tk()
		root.withdraw()
		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with NETCDF-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.nc'):
					list.append(os.path.join(root,filename))
	
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		while True:
			def safe():
				print(e1.get())
				print(e2.get())	
				Col = e1.get()
				Row = e2.get()				
				master.quit()
				
			master = Tk()
			Label(master, text = "Column number of the target cell", justify = LEFT).grid(row=0)
			Label(master, text = "Row number of the target cell", justify = LEFT).grid(row=1)
		
			e1 = Entry(master)
			e2 = Entry(master)
			e1.grid(row=0, column=1)
			e2.grid(row=1, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			Col = e1.get()
			Row = e2.get()			
			break

		Outfile = os.path.join(dirpath,'Results','PixelPrec.csv') 				
		outcsv = open(Outfile,"wb")
		writer = csv.writer(outcsv, delimiter = ",")		
		writer.writerow(['Name','Precipitation [mm]'])
		
		for x in range(0, len(list)):
			netcdffilename = list[x]
			file = Dataset(list[x],'r')
			values = file.variables['prediction'][:,:]		
		
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(netcdffilename)
			basename = os.path.basename(netcdffilename)
			basenametxt = os.path.splitext(basename)[0]	
	
			print(basename)
			# The data has to be converted into numpy arrays to make them ready for use.
			np_data = np.array(values)
			np_data = np_data.flatten()
			np_data = np_data[::-1]	
			New = np_data.reshape(315,266)
			New = np.fliplr(New)		
			
			Value = New[Row,Col]
			writer.writerow([basenametxt, Value])
			
		print("End of this script.")

	if DataQuestion == "5":
		########################################################
		# Open the hdf5-files for the OPERA-dataset of Europe.
		########################################################
		root = Tkinter.Tk()
		root.withdraw()
		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.h5'):
					list.append(os.path.join(root,filename))
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		while True:
			def safe():
				print(e1.get())
				print(e2.get())	
				Col = e1.get()
				Row = e2.get()				
				master.quit()
				
			master = Tk()
			Label(master, text = "Column number of the target cell", justify = LEFT).grid(row=0)
			Label(master, text = "Row number of the target cell", justify = LEFT).grid(row=1)
		
			e1 = Entry(master)
			e2 = Entry(master)
			e1.grid(row=0, column=1)
			e2.grid(row=1, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			Col = e1.get()
			Row = e2.get()			
			break
			
		# Ask if the user wants to change reflectivities into mm.
		
		def safe():
			print(e1.get())
			CalcQues = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "Do you want to convert reflectivities [dBZ] into [mm] of rain? Note: only use this").grid(row=0)
		Label(master, text = "tool when the provided data has reflectivities in dBZ. Please give the number of your answer: [1] Yes, [2] No and press enter").grid(row=1)
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		CalcQues = e1.get()				
		
		if CalcQues == "1":
			def safe():
				print(e1.get())
				NumMin = e1.get()		
				master.quit()
				
			master = Tk()
			Label(master, text = "What is the temporal resolution of your data (e.g. 5 min. data)? Please give the amount in minutes.").grid(row=0)
		
			e1 = Entry(master)
			e1.grid(row=0, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			NumMin = e1.get()
			NumMin = float(NumMin)
			
		else:
			print("No conversion will take place.")				
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		Outfile = os.path.join(dirpath,'Results','PixelPrec.csv') 				
		outcsv = open(Outfile,"wb")
		writer = csv.writer(outcsv, delimiter = ",")		
		writer.writerow(['Name','Precipitation [mm]'])
		
		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/image1/image_data'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
		
			if CalcQues == "1":
				filedata = 0.5 * filedata - 32
				filedata = ((10.**((filedata-23)/16))/(60/NumMin))
			if CalcQues == "2":
				filedata = filedata * 0.01	
	
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]		
			print(basenametxt)		
	
			# Multiply the data with 0.01 to get mm.
			filedata = filedata * 0.01	
			
			Value = filedata[Row,Col]
			writer.writerow([basenametxt, Value])
			
		print("End of this script.")
	
	if DataQuestion == "6":
		########################################################
		# Open the hdf5-files of the NASA GPM IMERG High Quality Precipitation level 3 dataset.
		########################################################

		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.HDF5'):
					list.append(os.path.join(root,filename))

		root = Tkinter.Tk()
		root.withdraw()	
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		while True:
			def safe():
				print(e1.get())
				print(e2.get())	
				Col = e1.get()
				Row = e2.get()				
				master.quit()
				
			master = Tk()
			Label(master, text = "Column number of the target cell", justify = LEFT).grid(row=0)
			Label(master, text = "Row number of the target cell", justify = LEFT).grid(row=1)
		
			e1 = Entry(master)
			e2 = Entry(master)
			e1.grid(row=0, column=1)
			e2.grid(row=1, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			Col = e1.get()
			Row = e2.get()			
			break

		Outfile = os.path.join(dirpath,'Results','PixelPrec.csv') 				
		outcsv = open(Outfile,"wb")
		writer = csv.writer(outcsv, delimiter = ",")		
		writer.writerow(['Name','Precipitation [mm]'])
		
		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/Grid/HQprecipitation'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
		
			filedata = filedata * 0.5 # The units are mm/h, while the temporal resolution is 30 min.
		
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(hdf5filename)
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]	
		
			print(basename)
			
			Value = filedata[Row,Col]
			writer.writerow([basenametxt, Value])
			
		print("End of this script.")


if FirstQues == 2:
	#########################################################
	# Time to read the shapefile
	#########################################################

	root = Tkinter.Tk()
	root.withdraw()
	filename = tkFileDialog.askopenfilename(parent=root,filetypes=[("ESRI Shapefile","*.shp"),("all files","*.*")],title='Choose a Shapefile of your study area')
	sf = shapefile.Reader(filename)

	#########################################################
	# Get the projection of the shapefiles
	##########################################################

	driver = ogr.GetDriverByName('ESRI Shapefile')
	dataset = driver.Open(filename)

	# from Layer
	layer = dataset.GetLayer()
	spatialRef = layer.GetSpatialRef()
	# from Geometry
	feature = layer.GetNextFeature()
	geom = feature.GetGeometryRef()
	spatialRef = geom.GetSpatialReference()
	
	print(spatialRef)

	#########################################################
	# Indicate which type of data set will be supplied
	#########################################################	

	def safe():
		print(e1.get())	
		DataQuestion = e1.get()		
		root2.quit()
				
	root2 = Tk()
	root2.title("Indicate which type of dataset will be supplied. Fill out the number of your choice and click on enter.")
	Label(root2, text = "1 NLRadarhdf5_NL25", justify = LEFT).grid(row=0)
	Label(root2, text = "2 NLRadarhdf5_NL21", justify = LEFT).grid(row=1)
	Label(root2, text = "3 DERadar_Radolan", justify = LEFT).grid(row=2)		
	Label(root2, text = "4 NL KDC data in NETCDF", justify = LEFT).grid(row=3)
	Label(root2, text = "5 OPERA - Radar Europe", justify = LEFT).grid(row=4)
	Label(root2, text = "6 NASA GPM IMERG_HQprecipitation level 3", justify = LEFT).grid(row=5)	
		
	e1 = Entry(root2)
	e1.grid(row=2, column=3)
		
	Button(root2, text = "Enter", command = safe).grid(row=3,column=3,sticky=W,pady=4)
						
	mainloop()							
				
	DataQuestion = e1.get()
	print(DataQuestion)		
	
	#########################################################
	# Reproject the geometry of the shapefile
	#########################################################

	# set file names
	infile = filename
	dirname = os.path.dirname(filename)
	basename = os.path.basename(filename)
	basenametxt = os.path.splitext(basename)[0]
	outfile = os.path.join(dirname, basenametxt+'_Reprojected.shp')			
	
	# The projection will be based on the projection of the given data set	
	
	if DataQuestion == "1":
		ShapefileProj = '+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0'
	if DataQuestion == "2":
		ShapefileProj = '+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0'
	if DataQuestion == "3":
		ShapefileProj = '+proj=stere +lat_0=90.0 +lat_ts=60.0 +lon_0=10.0 +a=6370040 +b=6370040 +units=m'
	if DataQuestion == "4":
		ShapefileProj = '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs'
	if DataQuestion == "5":
		ShapefileProj = '+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0'
	if DataQuestion == "6":
		ShapefileProj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

	
	subprocess.call(['ogr2ogr', '-t_srs', ShapefileProj, outfile, infile])
	
	#########################################################
	# Define the new, reprojected shapefile as sfnew
	# We are going to rasterize this shapefile, since the resulting radar raster (after
	# rasterizing the hdf5-dataset) should be exactly the same as this raster of this shapefile
	#########################################################
	
	sfnew = shapefile.Reader(outfile)

	# 1. Define pixel_size and NoData value of new raster
	NoData_value = -9999
	# Let's ask for the cell size and use that for the xres, yres and pixel_size

	def safe():
		print(e1.get())
		Res = e1.get()		
		master.quit()
				
	master = Tk()
	Label(master, text = "Please give the desired pixel resolution in km or", justify = LEFT).grid(row=0)
	Label(master, text = "in degrees, depending on the provided data set.", justify = LEFT).grid(row=1)	
	
	e1 = Entry(master)
	e1.grid(row=0, column=1)
			
	Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
	mainloop()							
						
	Res = e1.get()
	Res = float(Res)		
	
	if DataQuestion == "3":
		Res = Res * 1000
	if DataQuestion == "4":
		Res = Res * 1000		
	
	#Res = input("Please give the desired pixel resolution in km2")
	x_res = Res
	y_res = Res
	Pixel_size = Res

	# 2. Filenames for in- and output
	_out = os.path.join(dirname, basenametxt+'_Reprojected.tiff')

	# 3. Open Shapefile
	driver = ogr.GetDriverByName('ESRI Shapefile')
	shapefile = driver.Open(outfile)
	source_layer = shapefile.GetLayer()
	# Obtain xmin, xmax, ymin and ymax as floats
	x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
	# Round the values for xmin, xmax, ymin and ymax (make a rough rounding)
	# in order to not get a too small extent. Finally, make integer values of it.
	xmin = int(round(x_min, 0) - 1)
	xmax = int(round(x_max, 0) + 1)
	ymin = int(round(y_min, 0) - 1)
	ymax = int(round(y_max, 0) + 1)

	# 4. Create Target - TIFF
	cols = int( (xmax - xmin) / x_res )
	rows = int( (ymax - ymin) / y_res )

	# Save raster
	# TO DO: Make this an option, it's not necessary to save it in all cases.
	_raster = gdal.GetDriverByName('GTiff').Create(_out, cols, rows, 1, gdal.GDT_Byte)
	_raster.SetGeoTransform((xmin, x_res, 0, ymax, 0, -y_res))
	projection = osr.SpatialReference()
	projection.ImportFromProj4(ShapefileProj)
	_raster.SetProjection(projection.ExportToWkt())
	_band = _raster.GetRasterBand(1)
	_band.SetNoDataValue(NoData_value)

	# 5. Rasterize
	gdal.RasterizeLayer(_raster, [1], source_layer, None, None, [1], ['ALL_TOUCHED=TRUE'])

	inraster = _out

	########################################################
	# Open the dataset
	########################################################

	if DataQuestion == "1":
		########################################################
		# Open the hdf5-files for the Netherlands, in this case the NL_25 dataset.
		########################################################

		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.h5'):
					list.append(os.path.join(root,filename))

		print(dirpath)

		# Ask if the user wants to change reflectivities into mm.
		
		def safe():
			print(e1.get())
			CalcQues = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "Do you want to convert reflectivities [dBZ] into [mm] of rain? Note: only use this").grid(row=0)
		Label(master, text = "tool when the provided data has reflectivities in dBZ. Please give the number of your answer: [1] Yes, [2] No and press enter").grid(row=1)
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		CalcQues = e1.get()				
		
		if CalcQues == "1":
			def safe():
				print(e1.get())
				NumMin = e1.get()		
				master.quit()
				
			master = Tk()
			Label(master, text = "What is the temporal resolution of your data (e.g. 5 min. data)? Please give the amount in minutes.").grid(row=0)
		
			e1 = Entry(master)
			e1.grid(row=0, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			NumMin = e1.get()
			NumMin = float(NumMin)
			
		else:
			print("No conversion will take place.")	

		########################################################
		# It is time to start the 'For Loop'
		########################################################

		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/image1/image_data'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
		
			if CalcQues == "1":
				filedata = 0.5 * filedata - 32
				filedata = ((10.**((filedata-23)/16))/(60/NumMin))
			if CalcQues == "2":
				filedata = filedata * 0.01	
		
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(hdf5filename)
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]	
		
			print(basename)
			#########################################################
			# Create a raster from the hdf5 file
			#########################################################
	
			# Make a raster out of the hdf5 dataset
			# First convert the data into numpy arrays
			np_data = np.array(filedata)
	
			# Get the spatial extent of the data
			num_cols = float(np_data.shape[1])
			num_rows = float(np_data.shape[0])
			xmin = 0
			xmax = 699
			ymin = 0
			ymax = -3650 # offset is 3650
			xres = (xmax - xmin)/num_cols
			yres = (ymax - ymin)/num_rows

			# Time to set up the transformation which is necessary to create the raster we want
			geotransform = (xmin, xres, 0, ymin, 0, -yres)

			# We are going to create the raster with the right coordinate encoding and projection
			def array_to_raster(array):
				dst_filename = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')
				x_pixels = 700
				y_pixels = 765
				pixel_size = 1
				x_min = 0
				y_max = -3650
	
				driver = gdal.GetDriverByName('GTiff')
	
				dataset = driver.Create(
					dst_filename,
					x_pixels,
					y_pixels,
					1,
					gdal.GDT_Float32, )
		
				dataset.SetGeoTransform((
					x_min,
					pixel_size,
					0,
					y_max,
					0,
					-pixel_size))
		
				projection = osr.SpatialReference()
				projection.ImportFromProj4('+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0')
				dataset.SetProjection(projection.ExportToWkt())
				dataset.GetRasterBand(1).WriteArray(array)
				dataset.FlushCache() # This is the way to write it to the disk
				return dataset, dataset.GetRasterBand(1)

			array_to_raster(np_data)
		
			rasterfile = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')

			#########################################################
			# Time to clip the raster with the rasterized shapefile
			#########################################################

			#	inraster = _out --> This is already defined in the rasterizing of the shapefile.
			outraster = os.path.join(dirpath,'RasterIn', basenametxt+'Reprojected.tiff') # Here, the clip is saved as an ASCII file, but a .tiff extension is also possible.

			# Source
			src_filename = rasterfile
			src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
			src_proj = src.GetProjection()
			src_geotrans = src.GetGeoTransform()

			# We want a section of source that matches this:
			match_filename = inraster
			match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
			match_proj = match_ds.GetProjection()
			match_geotrans = match_ds.GetGeoTransform()
			wide = match_ds.RasterXSize
			high = match_ds.RasterYSize

			# Output / destination
			dst_filename = outraster
			dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
			dst.SetGeoTransform( match_geotrans )
			dst.SetProjection( match_proj)

			# Do the work. Check whether the pixel resolution is the same as the hdf5 dataset
			# resolution. If yes, use a nearest neighbour conversion. If not, use a cubic spline
			# conversion for downscaling and averaging for upscaling to get the best possible 
			# reprojection on the new grid.
			if Res == 1:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
			elif Res > 1:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Average)
			else:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)

			# We do have the hdf5 dataset as GeoTiff in the correct raster grid and extent. Now, it's
			# time to only keep the pixels that overlap with the rasterized shapefile.

			Result =  os.path.join(dirpath,'Results', basenametxt+'.tiff') 
			Catchment = gdal.Open(_out)	
	
			xSize = dst.RasterXSize

			ySize = dst.RasterYSize
	
			InputBand = dst.ReadAsArray(0,0,xSize,ySize)
			
			if CalcQues == "1":
				NDvalue = 0.5 * 255 - 32
				NDvalue = ((10.**((NDvalue-23)/16))/(60/NumMin))
				InputBand = np.where(InputBand < NDvalue, InputBand, -999)
			if CalcQues == "2":
				InputBand = np.where(InputBand < 655.35, InputBand, -999)					
			
			CatchmentBand = _raster.ReadAsArray(0,0,xSize,ySize)
	
			OutputData = np.where(CatchmentBand > 0, InputBand, -999)	
	
			Driver = gdal.GetDriverByName('GTiff')	
		
			OutputDataSet = Driver.CreateCopy(Result,dst)
	
			OutputDataSet.GetRasterBand(1).WriteArray(OutputData)

			src = None			
			
			try:
				os.remove(rasterfile)
			except WindowsError:
				print("One file could not be deleted, the tool continues.")

	if DataQuestion == "2":
		########################################################
		# Open the hdf5-files for the Netherlands, in this case NL21 dataset.
		########################################################

		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.h5'):
					list.append(os.path.join(root,filename))
		print(dirpath)

		# Ask if the user wants to change reflectivities into mm.
		
		def safe():
			print(e1.get())
			CalcQues = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "Do you want to convert reflectivities [dBZ] into [mm] of rain? Note: only use this").grid(row=0)
		Label(master, text = "tool when the provided data has reflectivities in dBZ. Please give the number of your answer: [1] Yes, [2] No and press enter").grid(row=1)
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		CalcQues = e1.get()				
		
		if CalcQues == "1":
			def safe():
				print(e1.get())
				NumMin = e1.get()		
				master.quit()
				
			master = Tk()
			Label(master, text = "What is the temporal resolution of your data (e.g. 5 min. data)? Please give the amount in minutes.").grid(row=0)
		
			e1 = Entry(master)
			e1.grid(row=0, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			NumMin = e1.get()
			NumMin = float(NumMin)
			
		else:
			print("No conversion will take place.")	
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################

		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/image1/image_data'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
	
			if CalcQues == "1":
				filedata = 0.5 * filedata - 32
				filedata = ((10.**((filedata-23)/16))/(60/NumMin))
			if CalcQues == "2":
				filedata = filedata * 0.01	
	
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(hdf5filename)
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]	
	
			print(basename)
			#########################################################
			# Create a raster from the hdf5 file
			#########################################################

			# Make a raster out of the hdf5 dataset
			# First convert the data into numpy arrays
			np_data = np.array(filedata)

			# Get the spatial extent of the data
			num_cols = float(np_data.shape[1])
			num_rows = float(np_data.shape[0])
			xmin = 0
			xmax = 255
			ymin = 0
			ymax = -3727.265 # offset is 1490.906, so with a pixel size of 2.5 km, the ymax becomes
			# 1490.906 * 2.5. 
			xres = (xmax - xmin)/num_cols
			yres = (ymax - ymin)/num_rows

			# Time to set up the transformation which is necessary to create the raster we want
			geotransform = (xmin, xres, 0, ymin, 0, -yres)

			# We are going to create the raster with the right coordinate encoding and projection
			def array_to_raster(array):
				dst_filename = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')
				x_pixels = 256
				y_pixels = 256
				pixel_size = 2.5
				x_min = 0
				y_max = -3727.265
	
				driver = gdal.GetDriverByName('GTiff')
	
				dataset = driver.Create(
					dst_filename,
					x_pixels,
					y_pixels,
					1,
					gdal.GDT_Float32, )
		
				dataset.SetGeoTransform((
					x_min,
					pixel_size,
					0,
					y_max,
					0,
					-pixel_size))
		
				projection = osr.SpatialReference()
				projection.ImportFromProj4('+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0')
				dataset.SetProjection(projection.ExportToWkt())
				dataset.GetRasterBand(1).WriteArray(array)
				dataset.FlushCache() # This is the way to write it to the disk
				return dataset, dataset.GetRasterBand(1)

			array_to_raster(np_data)
	
			rasterfile = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')

			#########################################################
			# Time to clip the raster with the rasterized shapefile
			#########################################################

			#	inraster = _out --> This is already defined in the rasterizing of the shapefile.
			hdf5raster = rasterfile
			outraster = os.path.join(dirpath,'RasterIn', basenametxt+'Reprojected.tiff') # Here, the clip is saved as an ASCII file, but a .tiff extension is also possible.

			# Source
			src_filename = hdf5raster
			src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
			src_proj = src.GetProjection()
			src_geotrans = src.GetGeoTransform()

			# We want a section of source that matches this:
			match_filename = inraster
			match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
			match_proj = match_ds.GetProjection()
			match_geotrans = match_ds.GetGeoTransform()
			wide = match_ds.RasterXSize
			high = match_ds.RasterYSize

			# Output / destination
			dst_filename = outraster
			dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
			dst.SetGeoTransform( match_geotrans )
			dst.SetProjection( match_proj)

			# Do the work. Check whether the pixel resolution is the same as the hdf5 dataset
			# resolution. If yes, use a nearest neighbour conversion. If not, use a cubic spline
			# conversion for downscaling and averaging for upscaling to get the best possible 
			# reprojection on the new grid.d.
			if Res == 2.5:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
			elif Res > 2.5:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Average)
			else:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)

			# We do have the hdf5 dataset as GeoTiff in the correct raster grid and extent. Now, it's
			# time to only keep the pixels that overlap with the rasterized shapefile.

			Result =  os.path.join(dirpath,'Results', basenametxt+'.tiff') 
	
			Catchment = gdal.Open(_out)	
	
			xSize = dst.RasterXSize
			ySize = dst.RasterYSize
	
			InputBand = dst.ReadAsArray(0,0,xSize,ySize)
			
			if CalcQues == "1":
				NDvalue = 0.5 * 255 - 32
				NDvalue = ((10.**((NDvalue-23)/16))/(60/NumMin))
				InputBand = np.where(InputBand < NDvalue, InputBand, -999)
			if CalcQues == "2":
				InputBand = np.where(InputBand < 655.35, InputBand, -999)				
			
			CatchmentBand = _raster.ReadAsArray(0,0,xSize,ySize)
	
			OutputData = np.where(CatchmentBand > 0, InputBand, -999)	
	
			Driver = gdal.GetDriverByName('GTiff')	
	
			OutputDataSet = Driver.CreateCopy(Result,dst)
	
			OutputDataSet.GetRasterBand(1).WriteArray(OutputData)
			
			src = None			
			
			try:
				os.remove(rasterfile)
			except WindowsError:
				print("One file could not be deleted, the tool continues.")			


	if DataQuestion == "3":
		########################################################
		# Open the radolan radar dataset for Germany.
		########################################################

		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with ASCII-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.asc'):
					list.append(os.path.join(root,filename))
		print(dirpath)

		########################################################
		# It is time to start the 'For Loop'
		########################################################

		for x in range(0, len(list)):
			radolanfilename = list[x]
			filedata = np.loadtxt(list[x], skiprows=6)
	
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(radolanfilename)
			basename = os.path.basename(radolanfilename)
			basenametxt = os.path.splitext(basename)[0]	
	
			print(basename)
			#########################################################
			# Create a raster from the hdf5 file
			#########################################################

			# Make a raster out of the hdf5 dataset
			# First convert the data into numpy arrays
			np_data = filedata

			# Get the spatial extent of the data
			num_cols = 900
			num_rows = 900
			xmin = -523462
			xmax = 376538
			ymin = -4658645
			ymax = -3758645 
			xres = (xmax - xmin)/num_cols
			yres = (ymax - ymin)/num_rows

			# Time to set up the transformation which is necessary to create the raster we want
			geotransform = (xmin, xres, 0, ymin, 0, -yres)

			# We are going to create the raster with the right coordinate encoding and projection
			def array_to_raster(array):
				dst_filename = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')
				x_pixels = 900
				y_pixels = 900	
				pixel_size = 1000
				x_min = -523462
				y_max = -3758645
	
				driver = gdal.GetDriverByName('GTiff')
	
				dataset = driver.Create(
					dst_filename,
					x_pixels,
					y_pixels,
					1,
					gdal.GDT_Float32, )
		
				dataset.SetGeoTransform((
					x_min,
					pixel_size,
					0,
					y_max,
					0,
					-pixel_size))
		
				projection = osr.SpatialReference()
				projection.ImportFromProj4('+proj=stere +lat_0=90.0 +lat_ts=60.0 +lon_0=10.0 +a=6370040 +b=6370040 +units=m')
				dataset.SetProjection(projection.ExportToWkt())
				dataset.GetRasterBand(1).WriteArray(array)
				dataset.FlushCache() # This is the way to write it to the disk
				return dataset, dataset.GetRasterBand(1)

			array_to_raster(np_data)
	
			rasterfile = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')

			#########################################################
			# Time to clip the raster with the rasterized shapefile
			#########################################################

			#	inraster = _out --> This is already defined in the rasterizing of the shapefile.
			radolanraster = rasterfile
			outraster = os.path.join(dirpath,'RasterIn', basenametxt+'Reprojected.tiff') # Here, the clip is saved as an ASCII file, but a .tiff extension is also possible.

			# Source
			src_filename = radolanraster
			src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
			src_proj = src.GetProjection()
			src_geotrans = src.GetGeoTransform()

			# We want a section of source that matches this:
			match_filename = inraster
			match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
			match_proj = match_ds.GetProjection()
			match_geotrans = match_ds.GetGeoTransform()
			wide = match_ds.RasterXSize
			high = match_ds.RasterYSize

			# Output / destination
			dst_filename = outraster
			dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
			dst.SetGeoTransform( match_geotrans )
			dst.SetProjection( match_proj)

			# Do the work. Check whether the pixel resolution is the same as the hdf5 dataset
			# resolution. If yes, use a nearest neighbour conversion. If not, use a cubic spline
			# conversion for downscaling and averaging for upscaling to get the best possible 
			# reprojection on the new grid.
			if Res == 1:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
			elif Res > 1:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Average)
			else:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)
	
			# We do have the hdf5 dataset as GeoTiff in the correct raster grid and extent. Now, it's
			# time to only keep the pixels that overlap with the rasterized shapefile.

			Result =  os.path.join(dirpath,'Results', basenametxt+'.tiff') 
	
			Catchment = gdal.Open(_out)	
	
			xSize = dst.RasterXSize
			ySize = dst.RasterYSize
	
			InputBand = dst.ReadAsArray(0,0,xSize,ySize)
			InputBand = np.where(InputBand >= 0.0, InputBand, -999)
			CatchmentBand = _raster.ReadAsArray(0,0,xSize,ySize)
	
			OutputData = np.where(CatchmentBand > 0, InputBand, -999)	
	
			Driver = gdal.GetDriverByName('GTiff')	
	
			OutputDataSet = Driver.CreateCopy(Result,dst)
	
			OutputDataSet.GetRasterBand(1).WriteArray(OutputData)

			src = None

			try:
				os.remove(rasterfile)
			except WindowsError:
				print("One file could not be deleted, the tool continues.")

	if DataQuestion == "4":
		########################################################
		# Open the NETCDF-files for the Netherlands. In this case that is the monthly averaged precipitation, but it could also be hourly or five minute data.
		########################################################

		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with NETCDF-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.nc'):
					list.append(os.path.join(root,filename))
		print(dirpath)

		########################################################
		# It is time to start the 'For Loop'
		########################################################

		for x in range(0, len(list)):
			netcdffilename = list[x]
			file = Dataset(list[x],'r')
			lons = file.variables['lon'][:]
			lats = file.variables['lat'][:]
			values = file.variables['prediction'][:,:]				
		
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(netcdffilename)
			basename = os.path.basename(netcdffilename)
			basenametxt = os.path.splitext(basename)[0]	
	
			print(basename)
			#########################################################
			# Create a raster from the NETCDF file
			#########################################################

			### Make a raster out of the NETCDF dataset

			# First convert the data into numpy arrays

			np_data = np.array(values)
			np_data = np_data.flatten()
			np_data = np_data[::-1]	
			New = np_data.reshape(315,266)
			New = np.fliplr(New)	
		
			# Get the spatial extent of the data
			num_cols = 266
			num_rows = 315
			xmin = 12621.630033977
			xmax = 278621.630033977
			ymin = 305583.0457758
			ymax = 620583.0457758
			xres = (xmax - xmin)/num_cols
			yres = (ymax - ymin)/num_rows


			# Time to set up the transformation which is necessary to create the raster we want
			geotransform = (xmin, xres, 0, ymin, 0, -yres)

			# We are going to create the raster with the right coordinate encoding and projection
			def array_to_raster(array):
				dst_filename = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')
				x_pixels = 266
				y_pixels = 315
				pixel_size = 1000
				x_min = xmin
				y_max = ymax
	
				driver = gdal.GetDriverByName('GTiff')
	
				dataset = driver.Create(
					dst_filename,
					x_pixels,
					y_pixels,
					1,
					gdal.GDT_Float32, )
		
				dataset.SetGeoTransform((
					x_min,
					pixel_size,
					0,
					y_max,
					0,
					-pixel_size))
		
				projection = osr.SpatialReference()
				projection.ImportFromProj4('+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs')
				dataset.SetProjection(projection.ExportToWkt())
				dataset.GetRasterBand(1).WriteArray(array)
				dataset.FlushCache() # This is the way to write it to the disk
				return dataset, dataset.GetRasterBand(1)

			array_to_raster(New)
	
			rasterfile = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')

			#########################################################
			# Time to clip the raster with the rasterized shapefile
			#########################################################

			#	inraster = _out --> This is already defined in the rasterizing of the shapefile.
			hdf5raster = rasterfile
			outraster = os.path.join(dirpath,'RasterIn', basenametxt+'Reprojected.tiff') # Here, the clip is saved as an ASCII file, but a .tiff extension is also possible.
	
			# Source
			src_filename = hdf5raster
			src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
			src_proj = src.GetProjection()
			src_geotrans = src.GetGeoTransform()

			# We want a section of source that matches this:
			match_filename = inraster
			match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
			match_proj = match_ds.GetProjection()
			match_geotrans = match_ds.GetGeoTransform()
			wide = match_ds.RasterXSize
			high = match_ds.RasterYSize

			# Output / destination
			dst_filename = outraster
			dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
			dst.SetGeoTransform( match_geotrans )
			dst.SetProjection( match_proj)

			# Do the work. Check whether the pixel resolution is the same as the hdf5 dataset
			# resolution. If yes, use a nearest neighbour conversion. If not, use a cubic spline
			# conversion for downscaling and averaging for upscaling to get the best possible 
			# reprojection on the new grid.
			if Res == 1:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
			elif Res > 1:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Average)
			else:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)

			# We do have the hdf5 dataset as GeoTiff in the correct raster grid and extent. Now, it's
			# time to only keep the pixels that overlap with the rasterized shapefile.

			Result =  os.path.join(dirpath,'Results', basenametxt+'.tiff') 

			Catchment = gdal.Open(_out)	
	
			xSize = dst.RasterXSize
			ySize = dst.RasterYSize
		
			InputBand = dst.ReadAsArray(0,0,xSize,ySize)
			InputBand = np.where(InputBand >= 0.0, InputBand, -999)
			CatchmentBand = _raster.ReadAsArray(0,0,xSize,ySize)
	
			OutputData = np.where(CatchmentBand > 0, InputBand, -999)	
	
			Driver = gdal.GetDriverByName('GTiff')	
	
			OutputDataSet = Driver.CreateCopy(Result,dst)
	
			OutputDataSet.GetRasterBand(1).WriteArray(OutputData)

			src = None

			try:
				os.remove(rasterfile)
			except WindowsError:
				print("One file could not be deleted, the tool continues.")

	if DataQuestion == "5":
		########################################################
		# Open the hdf5-files of the OPERA-dataset for Europe
		########################################################

		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.h5'):
					list.append(os.path.join(root,filename))
		print(dirpath)
		
		# Ask if the user wants to change reflectivities into mm.
		
		def safe():
			print(e1.get())
			CalcQues = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "Do you want to convert reflectivities [dBZ] into [mm] of rain? Note: only use this").grid(row=0)
		Label(master, text = "tool when the provided data has reflectivities in dBZ. Please give the number of your answer: [1] Yes, [2] No and press enter").grid(row=1)
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		CalcQues = e1.get()				
		
		if CalcQues == "1":
			def safe():
				print(e1.get())
				NumMin = e1.get()		
				master.quit()
				
			master = Tk()
			Label(master, text = "What is the temporal resolution of your data (e.g. 5 min. data)? Please give the amount in minutes.").grid(row=0)
		
			e1 = Entry(master)
			e1.grid(row=0, column=1)
			
			Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
			mainloop()							
						
			NumMin = e1.get()
			NumMin = float(NumMin)
			
		else:
			print("No conversion will take place.")				
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/image1/image_data'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
		
			if CalcQues == "1":
				filedata = 0.5 * filedata - 32
				filedata = ((10.**((filedata-23)/16))/(60/NumMin))
			if CalcQues == "2":
				filedata = filedata * 0.01	
		
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(hdf5filename)
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]	
		
			print(basename)
			#########################################################
			# Create a raster from the hdf5 file
			#########################################################
	
			# Make a raster out of the hdf5 dataset
			# First convert the data into numpy arrays
			np_data = np.array(filedata)
	
			# Get the spatial extent of the data
			num_cols = float(np_data.shape[1])
			num_rows = float(np_data.shape[0])
			xmin = -1799.963710073
			xmax = -600.963710073
			ymin = 0
			ymax = -1598.99739 # offset is -399.99738 - 1199
			xres = (xmax - xmin)/num_cols
			yres = (ymax - ymin)/num_rows

			# Time to set up the transformation which is necessary to create the raster we want
			geotransform = (xmin, xres, 0, ymin, 0, -yres)

			# We are going to create the raster with the right coordinate encoding and projection
			def array_to_raster(array):
				dst_filename = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')
				x_pixels = 1200
				y_pixels = 1200
				pixel_size = 4
				x_min = -1799.963710073
				y_max = -1598.99739
	
				driver = gdal.GetDriverByName('GTiff')
	
				dataset = driver.Create(
					dst_filename,
					x_pixels,
					y_pixels,
					1,
					gdal.GDT_Float32, )
		
				dataset.SetGeoTransform((
					x_min,
					pixel_size,
					0,
					y_max,
					0,
					-pixel_size))
		
				projection = osr.SpatialReference()
				projection.ImportFromProj4('+proj=stere +lat_0=90 +lon_0=0 +lat_ts=60 +a=6378.14 +b=6356.75 +x_0=0 y_0=0')
				dataset.SetProjection(projection.ExportToWkt())
				dataset.GetRasterBand(1).WriteArray(array)
				dataset.FlushCache() # This is the way to write it to the disk
				return dataset, dataset.GetRasterBand(1)

			array_to_raster(np_data)
		
			rasterfile = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')

			#########################################################
			# Time to clip the raster with the rasterized shapefile
			#########################################################

			#	inraster = _out --> This is already defined in the rasterizing of the shapefile.
			hdf5raster = rasterfile
			outraster = os.path.join(dirpath,'RasterIn', basenametxt+'Reprojected.tiff') # Here, the clip is saved as an ASCII file, but a .tiff extension is also possible.
	
			# Source
			src_filename = hdf5raster
			src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
			src_proj = src.GetProjection()
			src_geotrans = src.GetGeoTransform()

			# We want a section of source that matches this:
			match_filename = inraster
			match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
			match_proj = match_ds.GetProjection()
			match_geotrans = match_ds.GetGeoTransform()
			wide = match_ds.RasterXSize
			high = match_ds.RasterYSize

			# Output / destination
			dst_filename = outraster
			dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
			dst.SetGeoTransform( match_geotrans )
			dst.SetProjection( match_proj)

			# Do the work. Check whether the pixel resolution is the same as the hdf5 dataset
			# resolution. If yes, use a nearest neighbour conversion. If not, use a cubic spline
			# conversion for downscaling and averaging for upscaling to get the best possible 
			# reprojection on the new grid.
			if Res == 4:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
			elif Res > 4:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Average)
			else:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)

			# We do have the hdf5 dataset as GeoTiff in the correct raster grid and extent. Now, it's
			# time to only keep the pixels that overlap with the rasterized shapefile.

			Result =  os.path.join(dirpath,'Results', basenametxt+'.tiff') 

			Catchment = gdal.Open(_out)	
	
			xSize = dst.RasterXSize
			ySize = dst.RasterYSize
		
			InputBand = dst.ReadAsArray(0,0,xSize,ySize)
			
			if CalcQues == "1":
				NDvalue = 0.5 * 255 - 32
				NDvalue = ((10.**((NDvalue-23)/16))/(60/NumMin))
				InputBand = np.where(InputBand < NDvalue, InputBand, -999)
			if CalcQues == "2":
				InputBand = np.where(InputBand < 655.35, InputBand, -999)							
			
			CatchmentBand = _raster.ReadAsArray(0,0,xSize,ySize)
	
			OutputData = np.where(CatchmentBand > 0, InputBand, -999)	
	
			Driver = gdal.GetDriverByName('GTiff')	
	
			OutputDataSet = Driver.CreateCopy(Result,dst)
	
			OutputDataSet.GetRasterBand(1).WriteArray(OutputData)			

			src = None

			try:
				os.remove(rasterfile)
			except WindowsError:
				print("One file could not be deleted, the tool continues.")
			
	if DataQuestion == "6":
		########################################################
		# Open the hdf5-files of the NASA GPM IMERG High Quality Precipitation level 3 dataset.
		########################################################

		dirpath = tkFileDialog.askdirectory(parent=root, title='Open the main folder with hdf5-files and press OK')
		list = []
		for root, dirnames, filenames in os.walk(dirpath):
			for filename in filenames:
				if filename.endswith('.HDF5'):
					list.append(os.path.join(root,filename))
		print(dirpath)			
		
		########################################################
		# It is time to start the 'For Loop'
		########################################################
		for x in range(0, len(list)):
			hdf5filename = list[x]
			file = h5py.File(list[x],'r+')
			filedata = file['/Grid/HQprecipitation'][()] # Open "image1" dataset under the root group, 
			# this dataset consists of the accumulated 5 min precipitation amounts. 
		
			filedata = filedata * 0.5 # The units are mm/h, while the temporal resolution is 30 min.
		
			# Get the name of the file path and the file itself
			dirname = os.path.dirname(hdf5filename)
			basename = os.path.basename(hdf5filename)
			basenametxt = os.path.splitext(basename)[0]	
		
			print(basename)
			#########################################################
			# Create a raster from the hdf5 file
			#########################################################
	
			# Make a raster out of the hdf5 dataset
			# First convert the data into numpy arrays
			np_data = np.array(filedata)
			
			# The data has to be rotated counterclockwise due to the way the data provided.
			np_data = np.rot90(np_data, k=1)			
			
			
			# Get the spatial extent of the data
			num_cols = float(np_data.shape[0])
			num_rows = float(np_data.shape[1])
			xmin = -180
			xmax = 180
			ymin = -90
			ymax = 90
			xres = (xmax - xmin)/num_cols
			yres = (ymax - ymin)/num_rows

			# Time to set up the transformation which is necessary to create the raster we want
			geotransform = (xmin, xres, 0, ymin, 0, -yres)

			# We are going to create the raster with the right coordinate encoding and projection
			def array_to_raster(array):
				dst_filename = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')
				x_pixels = 3600
				y_pixels = 1800
				pixel_size = 0.1
				x_min = -180
				y_max = 90
	
				driver = gdal.GetDriverByName('GTiff')
	
				dataset = driver.Create(
					dst_filename,
					x_pixels,
					y_pixels,
					1,
					gdal.GDT_Float32, )
		
				dataset.SetGeoTransform((
					x_min,
					pixel_size,
					0,
					y_max,
					0,
					-pixel_size))
		
				projection = osr.SpatialReference()
				projection.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
				dataset.SetProjection(projection.ExportToWkt())
				dataset.GetRasterBand(1).WriteArray(array)
				dataset.FlushCache() # This is the way to write it to the disk
				return dataset, dataset.GetRasterBand(1)

			array_to_raster(np_data)
		
			rasterfile = os.path.join(dirpath,'RasterIn', basenametxt+'.tiff')


			#########################################################
			# Time to clip the raster with the rasterized shapefile
			#########################################################

			#	inraster = _out --> This is already defined in the rasterizing of the shapefile.
			hdf5raster = rasterfile
			outraster = os.path.join(dirpath,'RasterIn', basenametxt+'Reprojected.tiff') # Here, the clip is saved as an ASCII file, but a .tiff extension is also possible.

			# Source
			src_filename = hdf5raster
			src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
			src_proj = src.GetProjection()
			src_geotrans = src.GetGeoTransform()
			src_x = src.RasterXSize

			# We want a section of source that matches this:
			match_filename = inraster
			match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
			match_proj = match_ds.GetProjection()
			match_geotrans = match_ds.GetGeoTransform()
			wide = match_ds.RasterXSize
			high = match_ds.RasterYSize

			# Output / destination
			dst_filename = outraster
			dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
			dst.SetGeoTransform( match_geotrans )
			dst.SetProjection( match_proj)

			# Do the work. Check whether the pixel resolution is the same as the hdf5 dataset
			# resolution. If yes, use a nearest neighbour conversion. If not, use a cubic spline
			# conversion for downscaling and averaging for upscaling to get the best possible 
			# reprojection on the new grid.
			if src_x == wide:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
			elif wide > src_x:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Average)
			else:
				gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)

			# We do have the hdf5 dataset as GeoTiff in the correct raster grid and extent. Now, it's
			# time to only keep the pixels that overlap with the rasterized shapefile.

			Result =  os.path.join(dirpath,'Results', basenametxt+'.tiff') 

			Catchment = gdal.Open(_out)	
	
			xSize = dst.RasterXSize
			ySize = dst.RasterYSize
	
			InputBand = dst.ReadAsArray(0,0,xSize,ySize)
			InputBand = np.where(InputBand >= 0.0, InputBand, -999)
			CatchmentBand = _raster.ReadAsArray(0,0,xSize,ySize)
	
			OutputData = np.where(CatchmentBand > 0, InputBand, -999)	
	
			Driver = gdal.GetDriverByName('GTiff')	
		
			OutputDataSet = Driver.CreateCopy(Result,dst)
	
			OutputDataSet.GetRasterBand(1).WriteArray(OutputData)
	
			src = None	
	
			try:
				os.remove(rasterfile)
			except WindowsError:
				print("One file could not be deleted, the tool continues.")


	print("Done with this part.")
 	

	########################################################
	# Give the possibility to get statistics, such as the mean, for the output files
	########################################################

	while True:
		
		def safe():
			print(e1.get())
			StatQues = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "The clipping was successful, would you like to calculate the statistics per time interval (such as the mean) for", justify = LEFT).grid(row=0)
		Label(master, text = "your resulting files? Please give the number of your answer: [1] Yes, [2] No and press enter.", justify = LEFT).grid(row=1)

		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		StatQues = e1.get()					
	
		if StatQues == "1":
			# Ask for the correct path with the results.
			Resultpath = tkFileDialog.askdirectory(title='Open the main folder with the resulting GeoTiff-files and press OK')
			ListResult = []
			for root, dirnames, filenames in os.walk(Resultpath):
				for filename in filenames:
					if filename.endswith('.tiff'):
						ListResult.append(os.path.join(root,filename))

			########################################################
			# It is time to start the 'For Loop'
			########################################################
		
			Outfile = os.path.join(Resultpath,'Stats.csv') 				
			outcsv = open(Outfile,"wb")
			writer = csv.writer(outcsv, delimiter = ",")		
			writer.writerow(['Name','Minimum','Maximum','Mean'])

			Loopnum = 0
		
			for x in range(0, len(ListResult)):
				Loopnum = Loopnum + 1
				Percentage = (float(Loopnum)/float(len(ListResult)))*100
				print(Percentage, "% of statistics completed.")
				Infile = ListResult[x]
				Basename = os.path.basename(Infile)
				Basenametxt = os.path.splitext(Basename)[0]
			
				gtif = gdal.Open(Infile)
				srcband = gtif.GetRasterBand(1)
				Values = srcband.ReadAsArray()
				Values_new = np.where(Values < 0, np.NaN, Values)
				statsMean = np.nanmean(Values_new)
				statsMin = np.nanmin(Values_new)
				statsMax = np.nanmax(Values_new)
				writer.writerow([Basenametxt,statsMin,statsMax,statsMean])		
		
			print("The calculation was successful. Note that the rows in the csv-files might not be ordered chronologically.")		
			outcsv.flush()
			outcsv.close()
			break
			
		if StatQues == "2":
			print("You have chosen not to calculate statistics for your resulting files.")
			break
		
		else:
			print("The answer you gave, was not one of the possible answers. Please try again.")

			
	########################################################
	# Give the possibility to change the projection of the resulting files
	########################################################

	# Ask whether or not the user wants to change the projection
	while True:
	
		def safe():
			print(e1.get())
			ProjChange = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "The clipping was successful, the resulting GeoTiff files have the same projection as the input data set.", justify = LEFT).grid(row=0)
		Label(master, text = "Do you want to reproject your resulting files [note that this takes some additional calculation time and will delete the ", justify = LEFT).grid(row=1)
		Label(master, text = "files with the current projection after a successful reprojection]?", justify = LEFT).grid(row=2)
		Label(master, text = "Please give the number of your answer: [1] Yes, [2] No and press enter.", justify = LEFT).grid(row=3)
		
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		ProjChange = e1.get()	
	
	
		if ProjChange == "1":
			# Ask for desired projection

			def safe():
				print(e1.get())	
				DesProj = e1.get()		
				root2.quit()
				
			root2 = Tk()
			root2.title("Please give the desired spatial projection. Fill out the number of your choice and click on enter.")
			Label(root2, text = "1 WGS 84", justify = LEFT).grid(row=0)
			Label(root2, text = "2 RD_New Amersfoort", justify = LEFT).grid(row=1)
			Label(root2, text = "3 ETRS89 for Europe", justify = LEFT).grid(row=2)		
			Label(root2, text = "4 ETRS89/LCC Germany (N-E)", justify = LEFT).grid(row=3)
			Label(root2, text = "5 OSR-ORG:6971 UK National Grid", justify = LEFT).grid(row=4)
			Label(root2, text = "6 EPSG:2062 Madrid 1870 Spain", justify = LEFT).grid(row=5)	
			Label(root2, text = "7 WGS 84 Arctic Polar Stereographic", justify = LEFT).grid(row=6)	
		
			e1 = Entry(root2)
			e1.grid(row=2, column=3)
		
			Button(root2, text = "Enter", command = safe).grid(row=3,column=3,sticky=W,pady=4)
						
			mainloop()							
				
			DesProj = e1.get()		
			
			if DesProj == "1":
				ShapeProj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
			if DesProj == "2":
				ShapeProj = '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs' 
			if DesProj == "3":
				ShapeProj = '+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs'
			if DesProj == "4":
				ShapeProj = '+proj=lcc +lat_1=48.66666666666666 +lat_2=53.66666666666666 +lat_0=51 +lon_0=10.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' 
			if DesProj == "5":
				ShapeProj = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +a=6377563 +b=6356256.161 +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894 +units=m +no_defs' 
			if DesProj == "6":
				ShapeProj = '+proj=lcc +lat_1=40 +lat_0=40 +lon_0=0 +k_0=0.9988085293 +x_0=600000 +y_0=600000 +a=6378298.3 +b=6356657.142669561 +pm=madrid +units=m +no_defs ' 
			if DesProj == "7":
				ShapeProj == '+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0'	
	
			Resultpath = tkFileDialog.askdirectory(title='Open the main folder with the resulting GeoTiff-files and press OK')	
			ListResult = []
			for root, dirnames, filenames in os.walk(Resultpath):
				for filename in filenames:
					if filename.endswith('.tiff'):
						ListResult.append(os.path.join(root,filename))

			########################################################
			# It is time to start the 'For Loop'
			########################################################

			for x in range(0, len(ListResult)):
				Infile = ListResult[x]
				Basename = os.path.basename(Infile)
				Basenametxt = os.path.splitext(Basename)[0]
				Outfile = os.path.join(Resultpath, Basenametxt+'_Reprojected.tiff') 		

				subprocess.call(['gdalwarp', Infile, Outfile, '-s_srs', ShapefileProj, '-t_srs', ShapeProj])
				try:
					os.remove(Infile)
				except WindowsError:
					print("One file could not be deleted, the tool continues.")
			break
	
		if ProjChange == "2":
			print("You have chosen for no reprojection of the resulting files.")
			break
		
		else:
			print("The answer you gave, was not one of the possible answers. Please try again.")


	print(datetime.now() - startTime)

	########################################################
	# Make it possible to convert the GeoTiff into an ASCII file
	########################################################
	while True:
		# Ask whether or not the user wants to change the GeoTiff files into ASCII files
		def safe():
			print(e1.get())
			TypeChange = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "Do you want to change the resulting GeoTiff files into ASCII files [note that this", justify = LEFT).grid(row=1)
		Label(master, text = "takes some additional calculation time and will delete the GeoTiff files]?", justify = LEFT).grid(row=2)
		Label(master, text = "Please give the number of your answer: [1] Yes, [2] No and press enter.", justify = LEFT).grid(row=3)
	
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		TypeChange = e1.get()	
	
		if TypeChange == "1":	
	
			Resultpath = tkFileDialog.askdirectory(title='Open the main folder with the resulting GeoTiff-files and press OK')	
			ListResult = []
			for root, dirnames, filenames in os.walk(Resultpath):
				for filename in filenames:
					if filename.endswith('.tiff'):
						ListResult.append(os.path.join(root,filename))

			# It is time to start the 'For Loop'

			for x in range(0, len(ListResult)):
				Infile = ListResult[x]
				Basename = os.path.basename(Infile)
				Basenametxt = os.path.splitext(Basename)[0]
				Outfile = os.path.join(Resultpath, Basenametxt+'_Reprojected.asc') 		
		
				subprocess.call(['gdal_translate', '-of', 'AAIGrid', Infile, Outfile])
		
				try:
					os.remove(Infile)
				except WindowsError:
					print("One file could not be deleted, the tool continues.")
			
			break
	
		else:
			print("The GeoTiff files will not be changed into ASCII files")
			break
			
	########################################################
	# Give the possibility to plot a raster
	########################################################
	
	while True:
	
		def safe():
			print(e1.get())
			Question = e1.get()		
			master.quit()
				
		master = Tk()
		Label(master, text = "This is the end of the script. Do you want a visualisation on a map? Please give the number of your answer: [1] Yes, [2] No and press enter.").grid(row=0)
		
		e1 = Entry(master)
		e1.grid(row=0, column=1)
			
		Button(master, text = "Enter", command = safe).grid(row=3,column=1,sticky=W,pady=4)
						
		mainloop()							
						
		Question = e1.get()

		# One step further, the user gave an answer, so the script has to do something with it.
		# Yes, means give the user the possibility to pick a .tiff files from the results and show
		# this result on a map. No, means that we stop the script.

		if Question == "1":

			# We use sfnew as shapefile

			plt.figure()
			for shape in sfnew.shapeRecords():
				x = [i[0] for i in shape.shape.points[:]]
				y = [i[1] for i in shape.shape.points[:]]
				plt.plot(x,y)
	
			rasterfile = tkFileDialog.askopenfilename(title='Choose a resulting file for visualisation')
			dirname = os.path.dirname(rasterfile)
			basename = os.path.basename(rasterfile)
			basenametxt = os.path.splitext(basename)[0]
			outfile = os.path.join(dirname, basenametxt+'.png')
	
			ds = gdal.Open(rasterfile)			
			data = ds.ReadAsArray()
			Highest = np.amax(data)
			gt = ds.GetGeoTransform()
			proj = ds.GetProjection()
			im = plt.imshow(data, vmin=0, vmax = Highest, clim=(0, Highest))	
			
			plt.savefig(outfile, dpi=75)
			plt.show()
	
		if Question == "2":
			print("End of this script.")
			break

