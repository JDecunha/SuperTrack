#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#internal
import utils
import templates
import command
import commandcondition
#external
import random
import sys
import os


"""
Created on Fri Aug 13 13:41:33 2021

@author: joseph
"""

def determine_single_properties(templateString):
      
   print "\n ** MACRO PROPERTIES: ** \n"
   
   trackLibraryCommand = command.Command("Base folder for track library: ")
   trackLibraryCommand.AddCondition(commandcondition.StringCondition())
   trackLibrary = str(trackLibraryCommand.GetInput())
   
   trackEnergyCommand = command.Command("Which energy do you want to Superimpose in MeV: ")
   trackEnergyCommand.AddCondition(commandcondition.FloatCondition())
   trackEnergy = float(trackEnergyCommand.GetInput())
    
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
           
   particle_name = os.path.basename(trackLibrary);
   macroFileAdditionalInformation = ""
   if outputPathAppendix != " ":
       macroFileAdditionalInformation = outputPathAppendix.replace("/","")
       jobname = "%s_%sMeV_%snm_diameter_%s_x%s_oversamples" % (particle_name,trackEnergy,sphereDiameter,macroFileAdditionalInformation,nOversamples)
   else:
       jobname = "%s_%sMeV_%snm_diameter_x%s_oversamples" % (particle_name,trackEnergy,sphereDiameter,nOversamples)
   
   macro_filepath = "../macros/%s.ini" % jobname

   fileDirCommand = command.Command("Input directory for jobfile relative to the SuperTrack main directory: ")
   fileDirCommand.AddCondition(commandcondition.StringCondition())
   jobfiledir = "../" + str(fileDirCommand.GetInput())
    
   generate_macrofile(trackLibrary,trackEnergy,nThreads,nOversamples,sideLength,sphereDiameter,binType,numBins,lowerBin,upperBin, outputPathAppendix)
   generate_runfile(macro_filepath,walltime,jobname,jobfiledir,templateString)

   print "build complete."
   
   return 1

def generate_macrofile(trackLibrary, trackEnergy, nThreads, nOversamples, sideLength, sphereDiameter, binType, numBins, lowerBin, upperBin, outputPathAppendix):
    
    #determine the path for the input tracks to analyze
    inputPath = trackLibrary + "/" + str(trackEnergy) + "MeV/"
    
    #determine the macro filename and path
    particle_name = os.path.basename(trackLibrary);
    macroFileAdditionalInformation = ""
    if outputPathAppendix != " ":
        macroFileAdditionalInformation = outputPathAppendix.replace("/","")
        macro_name = "%s_%sMeV_%snm_diameter_%s_x%s_oversamples" % (particle_name,trackEnergy,sphereDiameter,macroFileAdditionalInformation,nOversamples)
    else:
        macro_name = "%s_%sMeV_%snm_diameter_x%s_oversamples" % (particle_name,trackEnergy,sphereDiameter,nOversamples)
        
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
        

    return macro_filepath

def generate_runfile(macro,walltime,jobname,jobdir,templateString):
        
    utils.make_directory(jobdir)

    # Render the template
    runfile_template_filled = templates.run_command_template.format(macro=macro)
    seadragon_template_filled = templateString[0].format(walltime_request=walltime,job_name=jobname,jobdir=jobdir,run_command = runfile_template_filled)

    with file("%s/%s%s" % (jobdir,jobname,templateString[1]) , "w") as f:
        f.write(seadragon_template_filled)

    return 1
