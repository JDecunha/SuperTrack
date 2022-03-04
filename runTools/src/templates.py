#!/usr/bin/env python2
# -*- coding: utf-8 -*-

macro_template = """; Example macro file for SuperTrack

[Run]
Type = Folder
Path = {inputPath}
FirstFile = -1 ; #-1 is the default case, will analyze whole folder
LastFile = -1 ;
Threads = {nThreads}
Oversamples = {nOversamples}

[Simulation]             
Method = VoxelConstrainedSphere

[VoxelConstrainedSphere]
;Geometry parameters
ScoringRegionHalfLength = {scoringHalfLength}    ; side length in nanometers
ScoringSphereDiameter = {scoringSphereDiameter}        ; sphere diameter in nanometers

;Meta parameters
SuggestedCudaBlocks = 256
SuggestedCudaThreads = 32


[Histogram]
Type = {histType} 
NBins = {nBins}
BinLower = {lowerBin}
BinUpper = {upperBin}

;Meta Parameters
SuggestedCudaAccumulateBlocks = 64
SuggestedCudaAccumulateThreads = 32

[Output]
Path = {outputPath}
NamePrefix = {outputName}"""

run_command_template = """../build/SuperTrack {macro} \n"""
