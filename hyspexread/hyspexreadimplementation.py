# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 19:41:23 2019
This class allows to Hyspex data using python.
@author: Giacomo Marchioro
"""
import struct
import ast
from array import array
import numpy as np
import os 



class HySpex():
    def __init__(self,path):
        self.path = path
        self.name = path.split(os.sep)[0]
        self.binheader = None
        self.txtheader = None
        self.RE = None
        self.QE = None
        self.background = None
        self.spectral_calib = None
        self.bad_pixels = None
        self.data = None
        self.inizialize()

  



    numberofbytes = [6,2,4,4,200,120,8,4,4,56,4,4,64,200,200,
                    200,200,200,32,32,200,200,1] + [4]*41 + [8]*9

    position = np.cumsum([0] + numberofbytes)
    
    def inizialize(self):
        def su(dtype,f):
            """Struct unpack value from file handler.
            return a value not a tuple"""
            types = {'I':4,'B':1,'d':8}
            return struct.unpack(dtype,f.read(types[dtype]))[0]

        def b2t(bytesn,f):
            """Struct unpack value from file handler.
            return a value not a tuple"""
            types = {4:'I',1:'B',8:'d'}
            return struct.unpack(types[bytesn],f.read(bytesn))[0]
        
        lens = {0:"window",
            1: "20cm",
            2 : "1m",
            3 : "2m",
            4 : "3m"}


        ScanningMode = { 0 : "none",
                    1 : "rotation",
                    2 : "translation"}
        binfile = "".join(self.path.split('.')[:-1]+[".hyspex"])
        with open(binfile,'rb') as f:
            binheader = {
            "header" : f.read(6), # 0 - 6
            "code" : "%i%i" %(su('B',f),su('B',f)), # (6) 7-8
            "header offset" : su('I',f), # (8) 7-8
            "Serialnumber" : su('I',f), # (12)
            "configgile" : f.read(200), # (16)
            "settingfile" : f.read(120), # (216)
            "scaling_factor" : su('d',f), # (336) 
            "electronics":su('I',f), # (338)
            "comsetting_electronics":su('I',f),
            "comport_electronics":f.read(56),
            "fan_speed": su('I',f),
            "backTemperature":su('I',f),
            "comport":f.read(64),
            "detectstring": f.read(200),
            "sensor": f.read(200),
            "framegrabber": f.read(200),
            "ID": f.read(200).rstrip('\x00'),
            "supplier": f.read(200),
            "left_gain": f.read(32),
            "right_gain": f.read(32),
            "Comment":f.read(200).rstrip('\x00'),
            "backgroundfile": f.read(200),
            "recordHD":f.read(1), #´ 1945
            # no idea what the following value is...
            "xcamera": f.read(4), 
            "serverindex":su('I',f),
            "comsettings":su('I',f),
            "Number of background":su('I',f),
            "bands":su('I',f), # spectral size
            "samples":su('I',f), # spatial size
            "Binning":su('I',f),
            "detected":su('I',f),
            "Integration time":su('I',f),
            "Frameperiod":su('I',f),
            "default_R":su('I',f), # default bands
            "default_G":su('I',f), # default bands
            "default_B":su('I',f), # default bands
            "bitshift":su('I',f),
            "temperature_offset":su('I',f),
            "shutter":su('I',f),
            "background_present":su('I',f),
            "power":su('I',f),
            "current":su('I',f),
            "bias":su('I',f),
            "bandwidth":su('I',f),
            "vin":su('I',f),
            "vref":su('I',f),
            "sensor_vin":su('I',f),
            "sensor_vref":su('I',f),
            "cooling_temperature":su('I',f),
            "window_start":su('I',f),
            "window_stop":su('I',f),
            "readout_time":su('I',f),
            "p":su('I',f),
            "i":su('I',f),
            "d":su('I',f),
            # it seems it is different from the lines in the hdr file
            "lines":su('I',f), # Number of frames
            "nobp":su('I',f),
            "dw":su('I',f),
            "EQ":su('I',f),
            "Lens":lens[su('I',f)], 
            "FOVexp":su('I',f),
            "Scanningmode":ScanningMode[su('I',f)], 
            "CalibAvailable":su('I',f),
            "NumberOfAvg":su('I',f), #   1945 + 41*4 = 2109
            "SF":su('d',f),  # 2117
            "Aperture size":su('d',f), # 2125
            "Pixelsize x":su('d',f),  
            "Pixelsize y":su('d',f), # start at 2133 end  2141
            "temperature":su('d',f), # 
            "max_framerate":su('d',f),  #  start at 2149 end  2157
            # it seems there are three empty slots we asign them assuming to be 
            # double to be sure anything is hiding here.
            "empty_1":su('d',f),
            "empty_2":su('d',f),
            "empty_3":su('d',f),
            }
            # Most important values
            bandsn = binheader["bands"]
            pixelsn = binheader["samples"]
            lines = binheader["lines"]
            # Start at 2181
            spectral_calib = array('d')
            spectral_calib.fromfile(f,bandsn)
            RE = array('d')
            RE.fromfile(f,bandsn*pixelsn)
            QE = array('d')
            QE.fromfile(f,bandsn)
            background = array('d')
            background.fromfile(f,bandsn*pixelsn)
            bad_pixels = array('I')
            bad_pixels.fromfile(f,binheader["nobp"])
            f.seek(binheader["header offset"])
            # from here start the data
            #band = array('H')
            # read one line of the first band
            #for i in range(lines):
            #    band.fromfile(f,pixelsn)
            #    f.seek((bandsn-1)*pixelsn*2,1) # 1 from current pos
        # Uniform binheader with txtheader notation
        binheader['Number of frames'] = binheader['lines']
        binheader['wavelength'] =  spectral_calib
        binheader['default bands']= [binheader['default_R'],
                                    binheader['default_G'],
                                    binheader['default_B']]

        # we read the ASCII header
        hdr = "".join(self.path.split('.')[:-1]+[".hdr"]) 
        with open(hdr,'r') as f:
            next(f)
            next(f)
            headertxt = dict()
            for i in f:
                if i !=" }\n":
                    values = i.split('=')
                    if len(values) == 1:
                        field,value  = values[0],''
                    else:
                        field,value  = values
                    value = value.strip().replace('{','[').replace('}',']')
                    try:
                        value = ast.literal_eval(value)
                    # except(SyntaxError, ValueError) as g:
                    except:
                        pass
                    headertxt[field.strip()] = value
                    
        self.binheader = binheader
        self.txtheader = headertxt
        self.RE = np.array(RE).reshape(bandsn,pixelsn)
        self.QE = np.array(QE)
        self.background = np.array(background).reshape(bandsn,pixelsn)
        self.spectral_calib = np.array(spectral_calib)
        self.bad_pixels = np.array(bad_pixels)
        self.data = np.memmap(binfile,dtype=np.uint16,mode='c', offset=binheader["header offset"],shape =(lines,bandsn,pixelsn))
                

    def test_heders(self):
        binmetadata = set(self.binheader.keys())
        textmetadata = set(self.txtheader.keys())            
        sharedmetadata = binmetadata.intersection(textmetadata)
        binhdrunique = binmetadata - textmetadata
        txthdrunique = textmetadata - binmetadata

        print("Differences bisides wavelength aproximation.")
        print("------------------------------------------")
        print("Field         |BinHeader     |TextHeader    ")
        print("------------------------------------------")
        for i in sharedmetadata:   
            if self.binheader[i] != self.txtheader[i] and i != "wavelength":
                print( '%-14s|%-14s|%-14s' %(i,self.binheader[i],self.txtheader[i]))
        return txthdrunique,binhdrunique
        

