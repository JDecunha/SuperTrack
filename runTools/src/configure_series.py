#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#internal
import templates
import utils
import command
import commandcondition
#external
import os

"""
Created on Fri Aug 13 13:41:33 2021

@author: joseph
"""

class SeriesPropertiesStruct:
    
    trackLibrary = None
    nThreads = None
    nOversamples = None
    sideLength = None
    sphereDiameter = None
    binType = None
    numBins = None
    lowerBin = None
    upperBin = None
    outputPathAppendix = None
    
    walltime = None
    jobfiledir = None
    
    Built = False
    
    def __init__(self):
        return
    
    def Build(self, templateString):
        
        if self.Built == False:
            self.FirstBuild(templateString)
        else:
            self.ReBuild(templateString)
        
    def ReBuild(self, templateString):
        
        print "\n Current build values:"
        print "Track Library: ", self.trackLibrary
        print "Number of threads: ", self.nThreads
        print "Number of oversamples: ", self.nOversamples
        print "Cube half length [nm]: ", self.sideLength
        print "1.) Target diameter: [nm]", self.sphereDiameter
        print "Histogram bin type ", self.binType
        print "Histogram number of bins: ", self.numBins
        print "Smallest bin: ", self.lowerBin
        print "Greatest bin: ", self.upperBin
        print "Extra info appended to output: ", self.outputPathAppendix
        print "Walltime: ", self.walltime
        print "Run file directory: ", self.jobfiledir
        
        rebuild_loop = 1
        
        while rebuild_loop == 1:
            
           rebuildCommand = command.Command("Number to change property, or rebuild: ")
           rebuildCommand.AddCondition(commandcondition.StringCondition())
           rebuild = str(rebuildCommand.GetInput())
           
           if rebuild == "1":
               
               sphereDiameterCommand = command.Command("Input scoring sphere diameter in nm: ")
               sphereDiameterCommand.AddCondition(commandcondition.FloatCondition())
               self.sphereDiameter = str(sphereDiameterCommand.GetInput())
               
           elif rebuild == "rebuild":
               
               energies = [name for name in os.listdir(self.trackLibrary)]
               particle_name = os.path.basename(self.trackLibrary);
               macronames = []
               macrofilepaths = []
               for energy in energies:  
                   if self.outputPathAppendix != " ":
                       macroFileAdditionalInformation = self.outputPathAppendix.replace("/","")
                       macroname = "%s_%s_%snm_diameter_%s_x%s_oversamples" % (particle_name,energy,self.sphereDiameter,macroFileAdditionalInformation,self.nOversamples)
                   else:
                       macroname = "%s_%s_%snm_diameter_x%s_oversamples" % (particle_name,energy,self.sphereDiameter,self.nOversamples)
                       
                   macronames.append(macroname)
                   macrofilepaths.append("../macros/%s.ini" % macroname)
                   
               if self.outputPathAppendix != " ":  
                   runFileName = "%s_%snm_diameter_%s_x%s_oversamples" % (particle_name,self.sphereDiameter,macroFileAdditionalInformation,self.nOversamples)
               else:
                   runFileName = "%s_%snm_diameter_x%s_oversamples" % (particle_name,self.sphereDiameter,self.nOversamples)
                       
               generate_macrofile_series(self.trackLibrary, energies, self.nThreads, self.nOversamples, self.sideLength, self.sphereDiameter, self.binType, self.numBins, self.lowerBin, self.upperBin, self.outputPathAppendix)
               generate_runfile_series(macrofilepaths,self.walltime,runFileName,self.jobfiledir,templateString)
               rebuild_loop = 0
               
           else:
               print "Try again mon ami."
        
        return 0
        
    def FirstBuild(self, templateString):
       print "\n ** MACRO PROPERTIES: ** \n"
       
       trackLibraryCommand = command.Command("Base folder for track library: ")
       trackLibraryCommand.AddCondition(commandcondition.StringCondition())
       self.trackLibrary = str(trackLibraryCommand.GetInput())
       
       energies = [name for name in os.listdir(self.trackLibrary)]
        
       nThreadsCommand = command.Command("Input number of CPU threads: ")
       nThreadsCommand.AddCondition(commandcondition.IntCondition())
       self.nThreads = int(nThreadsCommand.GetInput())
       
       nOversamplesCommand = command.Command("Input number of track oversamples: ")
       nOversamplesCommand.AddCondition(commandcondition.IntCondition())
       self.nOversamples = int(nOversamplesCommand.GetInput())
       
       sideLengthCommand = command.Command("Input world box half-length in nm: ")
       sideLengthCommand.AddCondition(commandcondition.FloatCondition())
       self.sideLength = str(sideLengthCommand.GetInput())
       
       sphereDiameterCommand = command.Command("Input scoring sphere diameter in nm: ")
       sphereDiameterCommand.AddCondition(commandcondition.FloatCondition())
       self.sphereDiameter = str(sphereDiameterCommand.GetInput())
       
       binTypeCommand = command.Command("Input histogram bin type (lin or log): ")
       binTypeCommand.AddCondition(commandcondition.StringCondition())
       self.binType = str(binTypeCommand.GetInput())
       
       numBinsCommand = command.Command("Input number of histogram bins: ")
       numBinsCommand.AddCondition(commandcondition.IntCondition())
       self.numBins = int(numBinsCommand.GetInput())
            
       lowerBinCommand = command.Command("Input lowest histogram bin (lower edge): ")
       lowerBinCommand.AddCondition(commandcondition.FloatCondition())
       self.lowerBin = float(lowerBinCommand.GetInput())
       
       upperBinCommand = command.Command("Input greatest histogram bin (upper edge): ")
       upperBinCommand.AddCondition(commandcondition.FloatCondition())
       self.upperBin = float(upperBinCommand.GetInput())
       
       outputPathAppendixCommand = command.Command("Include anything you would like to append to the output path (enter a space to skip): ")
       outputPathAppendixCommand.AddCondition(commandcondition.StringCondition())
       self.outputPathAppendix = str(outputPathAppendixCommand.GetInput())
       
       print "\n ** RUNFILE PROPERTIES: ** \n"
                        
       walltimeCommand = command.Command("Input walltime request: ")
       walltimeCommand.AddCondition(commandcondition.StringCondition())
       self.walltime = str(walltimeCommand.GetInput())
    
       fileDirCommand = command.Command("Input directory for jobfile relative to the SuperTrack main directory: ")
       fileDirCommand.AddCondition(commandcondition.StringCondition())
       self.jobfiledir = "../" + str(fileDirCommand.GetInput())
          
       particle_name = os.path.basename(self.trackLibrary);
       macronames = []
       macrofilepaths = []
       for energy in energies:  
           if self.outputPathAppendix != " ":
               macroFileAdditionalInformation = self.outputPathAppendix.replace("/","")
               macroname = "%s_%s_%snm_diameter_%s_x%s_oversamples" % (particle_name,energy,self.sphereDiameter,macroFileAdditionalInformation,self.nOversamples)
           else:
               macroname = "%s_%s_%snm_diameter_x%s_oversamples" % (particle_name,energy,self.sphereDiameter,self.nOversamples)
               
           macronames.append(macroname)
           macrofilepaths.append("../macros/%s.ini" % macroname)
           
       if self.outputPathAppendix != " ":  
           runFileName = "%s_%snm_diameter_%s_x%s_oversamples" % (particle_name,self.sphereDiameter,macroFileAdditionalInformation,self.nOversamples)
       else:
           runFileName = "%s_%snm_diameter_x%s_oversamples" % (particle_name,self.sphereDiameter,self.nOversamples)
               
       generate_macrofile_series(self.trackLibrary, energies, self.nThreads, self.nOversamples, self.sideLength, self.sphereDiameter, self.binType, self.numBins, self.lowerBin, self.upperBin, self.outputPathAppendix)
       generate_runfile_series(macrofilepaths,self.walltime,runFileName,self.jobfiledir,templateString)
       
       print "build complete."       
       self.Built = True
       
       return  0

def determine_series_properties(templateString):
    
   print "\n ** MACRO PROPERTIES: ** \n"
   
   trackLibraryCommand = command.Command("Base folder for track library: ")
   trackLibraryCommand.AddCondition(commandcondition.StringCondition())
   trackLibrary = str(trackLibraryCommand.GetInput())
   
   energies = [name for name in os.listdir(trackLibrary)]
    
   nThreadsCommand = command.Command("Input number of CPU threads: ")
   nThreadsCommand.AddCondition(commandcondition.IntCondition())
   nThreads = int(nThreadsCommand.GetInput())
   
   nOversamplesCommand = command.Command("Input number of track oversamples: ")
   nOversamplesCommand.AddCondition(commandcondition.IntCondition())
   nOversamples = int(nOversamplesCommand.GetInput())
   
   sideLengthCommand = command.Command("Input world box half-length in nm: ")
   sideLengthCommand.AddCondition(commandcondition.FloatCondition())
   sideLength = str(sideLengthCommand.GetInput())
   
   sphereDiameterCommand = command.Command("Input scoring sphere diameter in nm: ")
   sphereDiameterCommand.AddCondition(commandcondition.FloatCondition())
   sphereDiameter = str(sphereDiameterCommand.GetInput())
   
   binTypeCommand = command.Command("Input histogram bin type (lin or log): ")
   binTypeCommand.AddCondition(commandcondition.StringCondition())
   binType = str(binTypeCommand.GetInput())
   
   numBinsCommand = command.Command("Input number of histogram bins: ")
   numBinsCommand.AddCondition(commandcondition.IntCondition())
   numBins = int(numBinsCommand.GetInput())
        
   lowerBinCommand = command.Command("Input lowest histogram bin (lower edge): ")
   lowerBinCommand.AddCondition(commandcondition.FloatCondition())
   lowerBin = float(lowerBinCommand.GetInput())
   
   upperBinCommand = command.Command("Input greatest histogram bin (upper edge): ")
   upperBinCommand.AddCondition(commandcondition.FloatCondition())
   upperBin = float(upperBinCommand.GetInput())
   
   outputPathAppendixCommand = command.Command("Include anything you would like to append to the output path (enter a space to skip): ")
   outputPathAppendixCommand.AddCondition(commandcondition.StringCondition())
   outputPathAppendix = str(outputPathAppendixCommand.GetInput())
   
   print "\n ** RUNFILE PROPERTIES: ** \n"
                    
   walltimeCommand = command.Command("Input walltime request: ")
   walltimeCommand.AddCondition(commandcondition.StringCondition())
   walltime = str(walltimeCommand.GetInput())

   fileDirCommand = command.Command("Input directory for jobfile relative to the SuperTrack main directory: ")
   fileDirCommand.AddCondition(commandcondition.StringCondition())
   jobfiledir = "../" + str(fileDirCommand.GetInput())
      
   particle_name = os.path.basename(trackLibrary);
   macronames = []
   macrofilepaths = []
   for energy in energies:  
       if outputPathAppendix != " ":
           macroFileAdditionalInformation = outputPathAppendix.replace("/","")
           macroname = "%s_%s_%snm_diameter_%s_x%s_oversamples" % (particle_name,energy,sphereDiameter,macroFileAdditionalInformation,nOversamples)
       else:
           macroname = "%s_%s_%snm_diameter_x%s_oversamples" % (particle_name,energy,sphereDiameter,nOversamples)
           
       macronames.append(macroname)
       macrofilepaths.append("../macros/%s.ini" % macroname)
       
   if outputPathAppendix != " ":  
       runFileName = "%s_%snm_diameter_%s_x%s_oversamples" % (particle_name,sphereDiameter,macroFileAdditionalInformation,nOversamples)
   else:
       runFileName = "%s_%snm_diameter_x%s_oversamples" % (particle_name,sphereDiameter,nOversamples)
           
   generate_macrofile_series(trackLibrary, energies, nThreads, nOversamples, sideLength, sphereDiameter, binType, numBins, lowerBin, upperBin, outputPathAppendix)
   generate_runfile_series(macrofilepaths,walltime,runFileName,jobfiledir,templateString)
   
   print "build complete."        
   
   return  0

def generate_macrofile_series(trackLibrary, energies, nThreads, nOversamples, sideLength, sphereDiameter, binType, numBins, lowerBin, upperBin, outputPathAppendix):
    
   for energy in energies:

       #determine the path for the input tracks to analyze
       inputPath = trackLibrary + "/" + str(energy) + "/"
        
       #determine the macro filename and path
       particle_name = os.path.basename(trackLibrary);
       
       macroFileAdditionalInformation = ""
       if outputPathAppendix != " ":
           macroFileAdditionalInformation = outputPathAppendix.replace("/","")
           macro_name = "%s_%s_%snm_diameter_%s_x%s_oversamples" % (particle_name,energy,sphereDiameter,macroFileAdditionalInformation,nOversamples)
       else:
           macro_name = "%s_%s_%snm_diameter_x%s_oversamples" % (particle_name,energy,sphereDiameter,nOversamples)
      
       
       macro_filepath = "../macros/%s.ini" % macro_name
        
       #determine the SuperTrack output file path
       outputPath = "../output/%s" % particle_name
        
       #make the necessary folders and render the template
       utils.make_directory("../output/")
       utils.make_directory("../output/%s" % particle_name)
       
       if outputPathAppendix != " ":
           for newFolder in outputPathAppendix.split("/"):
               utils.make_directory("%s/%s" % (outputPath,newFolder))
               outputPath=outputPath+"/"+newFolder    
           outputPath+="/"
       else:
           outputPath+="/"
           
       macro_template_filled = templates.macro_template.format(inputPath=inputPath,nThreads=nThreads,nOversamples=nOversamples,scoringHalfLength=sideLength,scoringSphereDiameter=sphereDiameter,histType=binType,nBins=numBins,lowerBin=lowerBin,upperBin=upperBin,outputPath=outputPath,outputName=(macro_name+"_"),macro=macro_filepath)
        
       with file(macro_filepath, "w") as f:
           f.write(macro_template_filled)
        
   return 0

def generate_runfile_series(macros,walltime,jobname,jobdir,templateString):
        
    try:
        os.mkdir(jobdir)
    except OSError as err: 
        if err[0] == 17: #error no. 17 is file exists
            pass
        else: #any other type of error then raise it
            raise
            
    utils.make_directory(jobdir)
    with file("%s/%s%s" % (jobdir,jobname,templateString[1]) , "w") as f:
        seadragon_template_filled = templateString[0].format(walltime_request=walltime,job_name=jobname,jobdir=jobdir,run_command = "")
        f.write(seadragon_template_filled)
    
        for macro in macros:
    
            # Render the template
            run_command_template_filled = templates.run_command_template.format(macro=macro)
            f.write(run_command_template_filled)

    return 0

