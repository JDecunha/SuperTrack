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
        
   lowerBinCommand = command.Command("Input order of magnitude for lowest histogram bin: ")
   lowerBinCommand.AddCondition(commandcondition.IntCondition())
   lowerBin = int(lowerBinCommand.GetInput())
   
   upperBinCommand = command.Command("Input order of magnitude for greatest histogram bin: ")
   upperBinCommand.AddCondition(commandcondition.IntCondition())
   upperBin = int(upperBinCommand.GetInput())
   
   print "\n ** RUNFILE PROPERTIES: ** \n"
                    
   walltimeCommand = command.Command("Input walltime request: ")
   walltimeCommand.AddCondition(commandcondition.StringCondition())
   walltime = str(walltimeCommand.GetInput())
           
   particle_name = os.path.basename(trackLibrary);
   jobname = "%s_%sMeV_%snm_diameter_x%s_oversamples" % (particle_name,trackEnergy,sphereDiameter,nOversamples)
   macro_filepath = "../macros/%s.ini" % jobname

   fileDirCommand = command.Command("Input directory for jobfile relative to the SuperTrack main directory: ")
   fileDirCommand.AddCondition(commandcondition.StringCondition())
   jobfiledir = "../" + str(fileDirCommand.GetInput())
    
   generate_macrofile(trackLibrary,trackEnergy,nThreads,nOversamples,sideLength,sphereDiameter,lowerBin,upperBin)
   generate_runfile(macro_filepath,walltime,jobname,jobfiledir,templateString)

   print "build complete."
   
   return 1

def generate_macrofile(trackLibrary, trackEnergy, nThreads, nOversamples, sideLength, sphereDiameter, lowerBin, upperBin):
    
    #determine the path for the input tracks to analyze
    inputPath = trackLibrary + "/" + str(trackEnergy) + "MeV/"
    
    #determine the macro filename and path
    particle_name = os.path.basename(trackLibrary);
    macro_name = "%s_%sMeV_%snm_diameter_x%s_oversamples" % (particle_name,trackEnergy,sphereDiameter,nOversamples)
    macro_filepath = "../macros/%s.ini" % macro_name
    
    #determine the SuperTrack output file path
    outputPath = "../output/%s/" % particle_name
    
    #make the necessary folders and render the template
    utils.make_directory("../output/")
    utils.make_directory("../output/%s" % particle_name)    
    macro_template_filled = templates.macro_template.format(inputPath=inputPath,nThreads=nThreads,nOversamples=nOversamples,scoringHalfLength=sideLength,scoringSphereDiameter=sphereDiameter,lowerBin=lowerBin,upperBin=upperBin,outputPath=outputPath,outputName=(macro_name+"_"),macro=macro_filepath)

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
