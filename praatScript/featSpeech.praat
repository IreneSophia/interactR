####################################################
########## Irene Sophia Plank (c)2026 ##############
####### Includes code by Afton Nelson (c)2020 ######
# Includes the uhm-o-meter by de Jong et al (2021) #
####################################################
########## Preprocess audio .wav files #############
######## Each file containing one speaker ##########
####################################################

# CAUTION: if the output will be further processed 
# using the interactR script, then the filenames
# should contain the Dyad and Identifier. 
# E.g.: Recording_[Identifier1]-[Identifier2]_[Identifier].wav
# Thus, the information can be extracted in the
# interactR scripts. 


####################################################
############### Initialise settings ################
####################################################

#### [!PARAMETERS TO ADJUST!] ####
# specify directory containing the data and the output
folder$ = "audio"
# please enter the file separator from your OS
filesep$ = "/" 
# set a prefix for all the output files
prefix$ = "audio"

# set which processing steps should be run
procpint = 1
procuhm = 1

#### set variables

# set the pitch floor and ceiling - used to determine personal limits
pitch_floor = 50
pitch_ceil  = 700
pitch_step = 0.015

# set the timestep and pitch minimum for intensity
int_step = 0.01
int_pmin = 100

# add file separator to folder
folder$ = folder$ + filesep$

#### set uhm-o-meter parameters

# set file specification of the files to be included
fileSpec$ = folder$ + filesep$ + "*.wav"

# set the preprocessing options ["Reduce noise", "None", "Band pass (300..3300 Hz)"]
pre_processing$ = "Reduce noise"

# silence threshold: the more negative, the fewer sounding instances are detected
silence_threshold = -25
# recommended to set between 2 and 4
minimum_dip_near_peak = 2
minimum_pause_duration = 0.3
pitch_floor       = 30
# voicing threshold affects the syllable detection, not the sounding instances
voicing_threshold =  0.25
parser$           = "dip-peak-dip"

# which language? ["Dutch", "English"]
language$ = "English"

# set the options for filled pause detection
# if set to 1, download the script from 
# https://sites.google.com/view/uhm-o-meter/scripts/filledpauses-praat
# and save in the same folder
detect_Filled_Pauses = 0
filled_Pause_threshold = 1.00

# set the options for saving the results ["Praat Info window", "Save as text file", "Table", "TextGrids(s) only"]
# "Table" : TextGrid(s) + Tables are saved to disk
# "Praat Info window" : TextGrid(s) + output in Praat Info window
# "TextGrids(s) only" : TextGrid(s) are saved to disk
# "Save as text file" : TextGrid(s) + csv file containing articulation rate etc.
# "Save all" : TextGrid(s) + csv file + Table if filled pauses saved to disk
data$ = "Save all"

# if textfile or praat info, do you want to overwrite or append? ["OverWriteData", "AppendData"]
dataCollectionType$ = "OverWriteData" 

# keep the objects in praat? boolean
keep_Objects = 0


clearinfo

####################################################
####################################################
##### Finding pitch floor and ceiling per wav ######
####   Extract pitch and intensity composites   ####
####################################################
####################################################

if procpint

    appendInfoLine: "finding pitch floor and ceiling"
    
    # read files
    fls = Create Strings as file list... list 'folder$'*.wav
    numberOfFiles = Get number of strings
    
    # initialise table with columns
    limits = Create Table with column names: "limits", numberOfFiles, "name floor ceiling"

    # loop through the files
    for ifile to numberOfFiles
        
        select Strings list
        fileName$ = Get string... ifile
        appendInfoLine: "Determining limits ", ifile, " of ", numberOfFiles, ": ", fileName$
        sound = Read from file... 'folder$''fileName$'
        
        # get individual pitch floor and ceiling
        To Pitch: pitch_step, pitch_floor, pitch_ceil
        q1 = Get quantile: 0, 0, 0.25, "Hertz"
        q3 = Get quantile: 0, 0, 0.75, "Hertz"
        Remove
        floor = q1*0.75
        ceiling = q3*2.5
        floor = floor(floor)
        ceiling = ceiling(ceiling)
        
        removeObject: sound
        
        # add results to the table
        selectObject: limits
        Set string value: ifile, "name", fileName$
        Set numeric value: ifile, "floor", floor
        Set numeric value: ifile, "ceiling", ceiling
        
    endfor
    
    removeObject: fls
    
    appendInfoLine: "Done determining limits."
    appendInfoLine: "-----------------------------------------"

    # get information about limits to calculate pitch time step
    selectObject: limits
    min_floor = Get minimum: "floor"
    pstep = 0.75 / min_floor

    clearinfo

    # print a single header line with column names
    writeFileLine: "'folder$''prefix$'_pitchIntensity.csv", "Name,Duration,PitchFloor,PitchCeiling,PitchMean,PitchSD,PitchMin,PitchMax,IntensityMean,IntensitySD,IntensityMin,IntensityMax"
    
    appendInfoLine: "Extracting pitch and intensity composites"
    
    for i to numberOfFiles
        
        # select the table
        selectObject: limits
        # extract the information
        p_floor = object[limits, i, "floor"]
        p_ceil = object[limits, i, "ceiling"]
        fileName$ = object$[limits, i, "name"]
        
        # read in the file again
        sound = Read from file: folder$ + fileName$

        # selecting the sound object
        selectObject: sound

        # interpret start and end
        tmax = Get end time
        tmin = Get start time
        dur = tmax - tmin

        # extract short name
        shortName$ = selected$("Sound")
        appendInfoLine: "Extracting ", i, " of ", numberOfFiles, ": ", shortName$, "(", dur, " seconds)"
    
        # extracting pitch
        pitch = To Pitch (ac): pstep, p_floor, 4, "no", 0.03, 0.45, 0.01, 0.35, 0.14, p_ceil
        meanp = Get mean: tmin, tmax, "Hertz"
        sdp = Get standard deviation: tmin, tmax, "Hertz"
        minp = Get minimum: tmin, tmax, "Hertz", "Parabolic"
        maxp = Get maximum: tmin, tmax, "Hertz", "Parabolic"
    
        # extracting intensity
        selectObject: sound
        intensity = To Intensity: int_pmin, int_step, "yes"
        meani = Get mean: tmin, tmax, "energy"
        sdi = Get standard deviation: tmin, tmax
        mini = Get minimum: tmin, tmax, "Parabolic"
        maxi = Get maximum: tmin, tmax, "Parabolic"
    
        appendFileLine: "'folder$''prefix$'_pitchIntensity.csv", "'shortName$','dur','p_floor','p_ceil','meanp','sdp','minp','maxp','meani','sdi','mini','maxi'"
    
        removeObject: sound, pitch, intensity
    
    endfor

    removeObject: limits

    appendInfoLine: "Done extracting pitch and intensity."
    appendInfoLine: "-----------------------------------------"

endif


###########################################################################
#                                                                         #
#  Praat Script Syllable Nuclei, version 3 (Syllable Detector)            #
#  Copyright (C) 2019  Nivja de Jong, Ton Wempe, Jos J A Pacilly          #
#                                                                         #
#    This program is free software: you can redistribute it and/or modify #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         #
#    This program is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with this program.  If not, see http://www.gnu.org/licenses/   #
#                                                                         #
###########################################################################

if procuhm

  appendInfoLine: "Using the uhm-o-meter."

  @processArgs
  
  for file to nrFiles
    appendInfoLine: "Processing of file ", file, " named ", filename$[file]
    idSnd#[file] = Read from file: directory$ + filename$[file]
    @findSyllableNuclei
    ext$ = right$(filename$[file], length(filename$[file])-rindex(filename$[file], ".")+1)
    if detect_Filled_Pauses
      selectObject: idSnd#[file], idTG#[file]
      runScript: "FilledPauses.praat", language$, filled_Pause_threshold, (data$ == "Table" or data$ == "Save all")
      if data$ == "Table" or data$ == "Save all"
        idTbl#[file] = selected("Table")
        selectObject: idTbl#[file]
        Save as tab-separated file: directory$ + replace$(filename$[file], ext$, ".Table", 1)
        endif
      @countFilledPauses: idTG#[file]
    else
      @terminateLines
      endif
    selectObject: idTG#[file]
    Save as text file: directory$ + replace$(filename$[file], ext$, ".TextGrid", 1)
    if not keep_Objects
      removeObject: idSnd#[file], idTG#[file]
      if detect_Filled_Pauses and (data$ == "Table" or data$ == "Save all")
        removeObject: idTbl#[file]
        endif
      endif
    endfor
  
  @coda
  
  
  procedure findSyllableNuclei
    name$ = selected$("Sound")
  
    if pre_processing$ == "None"
      idSnd = selected ("Sound")
    elif pre_processing$ == "Band pass (300..3300 Hz)"
      idSnd = Filter (pass Hann band): 300, 3300, 100
      Scale peak: 0.99
    elif pre_processing$ == "Reduce noise"
      idSnd = noprogress Reduce noise: 0, 0, 0.025, 80, 10000, 40, -20, "spectral-subtraction"
      Scale peak: 0.99
      endif
  
    tsSnd  = Get start time
    teSnd  = Get end time
    dur    = Get total duration
  
    idInt  = To Intensity: 50, 0, "yes"		; use intensity to get threshold
    dbMin  = Get minimum: 0, 0, "Parabolic"
    dbMax  = Get maximum: 0, 0, "Parabolic"
    dbQ99  = Get quantile: 0, 0, 0.99		; get .99 quantile to get maximum (without influence of non-speech sound bursts)
  
  # estimate Intensity threshold
    threshold  = dbQ99 + silence_threshold
    threshold2 = dbMax - dbQ99
    threshold3 = silence_threshold - threshold2
    if threshold < dbMin
      threshold = dbMin
      endif
  
  # get pauses (silences) and speakingtime
    idTG = To TextGrid (silences): threshold3, minimum_pause_duration, 0.1, "", "sound"
    Set tier name: 1, "Phrases"
    nrIntervals = Get number of intervals: 1
    nsounding   = 0
    speakingtot = 0
    for interval to nrIntervals
      lbl$ = Get label of interval: 1, interval
      if lbl$ <> ""
        ts = Get start time of interval: 1, interval
        te = Get end time of interval:   1, interval
        nsounding   += 1
        speakingtot += te - ts
        Set interval text: 1, interval, string$(nsounding)
        endif
      endfor
  
    selectObject: idInt
    idPeak = To IntensityTier (peaks)
    selectObject: idSnd
    idP = noprogress To Pitch (ac): 0.02, pitch_floor, 4, "no", 0.03, voicing_threshold, 0.01, 0.35, 0.25, 450
  
  # fill array with intensity values
    peakcount = 0
    selectObject: idPeak
    nrPeaks = Get number of points
    for peak to nrPeaks
      selectObject: idPeak
      time  = Get time from index: peak
      dbMax = Get value at index: peak
      selectObject: idP
      voiced = Get value at time: time, "Hertz", "Linear"
      if dbMax > threshold and (voiced <> undefined)
        peakcount      += 1
        t [peakcount*2] = time		; peaks at EVEN indices (base 2)
        db[peakcount*2] = dbMax
        endif
      endfor
  
  # get Minima between peaks		; minima at ODD indices (base 1) t[1, 3, 5..]
    t[0]               = Get start time
    t[2*(peakcount+1)] = Get end time
    selectObject: idInt
    for valley to peakcount+1
      t [2*valley-1] = Get time of minimum: t[2*(valley-1)], t[2*valley], "Parabolic"
      db[2*valley-1] = Get minimum:         t[2*(valley-1)], t[2*valley], "Parabolic"
      endfor
  
    selectObject: idTG
    Insert point tier: 1, "Nuclei"
  
    tierDebug = 0
    if tierDebug
      Insert point tier: tierDebug, "DEBUG"
      endif
  
  # fill array with the largest peaks *followed* by a dip > minimum_dip_near_peak (obsolete), OR
  # with the largest peaks *surrounded* by a dip > minimum_dip_near_peak (current default method)
  
    voicedcount = 0	; nrNuclei
    tp[voicedcount] = t[0]
    tRise           = t[0]
    tFall           = t[0]
    tMax            = t[0]
    tMin            = t[0]
    dbMax           = db[1]
    dbMin           = db[1]
    nrPoints        = 2*peakcount+1
    selectObject: idTG
  
    for point to nrPoints
      if tierDebug
        Insert point: tierDebug, t[point], fixed$(db[point], 2)
        endif
  
      if db[point] > dbMax
        tMax  =  t[point]
        dbMax = db[point]
        if db[point] - dbMin > minimum_dip_near_peak
          tRise =  t[point]
          dbMin = db[point]
          endif
  
      elif db[point] < dbMin
        tMin  =  t[point]
        dbMin = db[point]
        if dbMax - db[point] > minimum_dip_near_peak
          tFall =  t[point]
          dbMax = db[point]
          endif
        endif
  
  #   Insert voiced peaks in TextGrid (note that the code for the obsolete
  #   "peak-dip" parser is kept only for backward compatibility reasons)
  
      if parser$ ==     "peak-dip" and                   tRise < tFall and tFall <> t[0] or
  ...    parser$ == "dip-peak-dip" and tRise <> t[0] and tRise < tFall and tFall <> t[0]
  
        i  = Get interval at time: 2, tMax
        l$ = Get label of interval: 2, i
        if l$ <> ""
          voicedcount    += 1
          tp[voicedcount] = tMax
          Insert point: 1, tMax, string$(voicedcount)
          tMax  = t [point]
          tMin  = t [point]
          dbMin = db[point]
          dbMax = db[point]
          tRise = t [0]
          tFall = t [0]
          endif
        endif
      endfor
    tp[voicedcount+1] = t[2*(peakcount+1)]
  
  # verify that shift in time between various objects is no longer an issue
    tsTG = Get start time
    teTG = Get end time
    assert tsSnd == tsTG and teSnd == teTG
  
  # clean up before next sound file is opened
    removeObject: idInt, idPeak, idP
    if pre_processing$ <> "None"
      removeObject: idSnd
      endif
    idTG#[file] = idTG
  
  # summarize results in Info window
    speakingrate = voicedcount / dur
    articulationrate = voicedcount / speakingtot
    npause = nsounding - 1
    asd = speakingtot / voicedcount
  
    if data$ == "Praat Info window"
      appendInfo: "'name$', 'voicedcount', 'npause', 'dur:2', 'speakingtot:2', 'speakingrate:2', 'articulationrate:2', 'asd:3'"
    elif data$ == "Save all"
      appendFile: txtFile$, "'name$', 'voicedcount', 'npause', 'dur:2', 'speakingtot:2', 'speakingrate:2', 'articulationrate:2', 'asd:3'"
      appendFile: temporaryDirectory$ + "/SyllableNuclei.tmp", "'name$', 'voicedcount', 'npause', 'dur:2', 'speakingtot:2', 'speakingrate:2', 'articulationrate:2', 'asd:3'"
    elif data$ == "Save as text file"
      appendFile: txtFile$, "'name$', 'voicedcount', 'npause', 'dur:2', 'speakingtot:2', 'speakingrate:2', 'articulationrate:2', 'asd:3'"
    elif data$ == "Table"
      appendFile: temporaryDirectory$ + "/SyllableNuclei.tmp", "'name$', 'voicedcount', 'npause', 'dur:2', 'speakingtot:2', 'speakingrate:2', 'articulationrate:2', 'asd:3'"
      endif
    endproc
  
  procedure countFilledPauses: .id
    selectObject: .id
    .nrInt = Get number of intervals: 3
    .nrFP  = 0
    .tFP   = 0
    for .int to .nrInt
      .lbl$ = Get label of interval: 3, .int
      if .lbl$ == "fp"
        .nrFP += 1
        .ts = Get start time of interval: 3, .int
        .te = Get end time of interval: 3, .int
        .tFP += (.te - .ts)
        endif
      endfor
    if data$ == "Praat Info window"
      appendInfoLine: ", '.nrFP', '.tFP:3'"
    elif data$ == "Save all"
      appendFileLine: txtFile$, ", '.nrFP', '.tFP:3'"
      appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ", '.nrFP', '.tFP:3'"
    elif data$ == "Save as text file"
      appendFileLine: txtFile$, ", '.nrFP', '.tFP:3'"
    elif data$ == "Table"
      appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ", '.nrFP', '.tFP:3'"
      endif
    endproc
  
  procedure terminateLines
    if data$ == "Praat Info window"
      appendInfoLine: ""
    elif data$ == "Save all"
      appendFileLine: txtFile$, ""
      appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ""
    elif data$ == "Save as text file"
      appendFileLine: txtFile$, ""
    elif data$ == "Table"
      appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ""
      endif
    endproc
  
  procedure processArgs


    if fileSpec$ <> ""
      len        = length(fileSpec$)
      sep        = rindex_regex(fileSpec$, "[\\/]")
      directory$ =  left$(fileSpec$, sep)
      selection$ = right$(fileSpec$, len-sep)
      idStr      = Create Strings as file list: "fileList", directory$ + selection$
      nrFiles    = Get number of strings
      nrObjects  = 0
      #appendInfoLine: directory$ + selection$
    else
      exit Unsupported Input Selection
      endif
  
    if nrFiles
      idSnd# = zero#(nrFiles)
      idTG#  = zero#(nrFiles)
      idTbl# = zero#(nrFiles)
      endif
    for file to nrFiles
      filename$[file] = Get string: file
      endfor
      
    if data$ == "Save as text file" or data$ == "Save all"
      txtFile$ = directory$ + prefix$ + "_syllableNuclei.csv"
    endif
  
    header$ = "Name,nSyll,nPause,Duration,PhonationDuration,SpeechRate,ArticulationRate,SyllableDuration"
    if data$ == "Praat Info window" and dataCollectionType$ == "OverWriteData"
      # print a single header line with column names and units
      writeInfo: header$
    elif data$ == "Save all" and dataCollectionType$ == "OverWriteData"
      writeFile: txtFile$, header$
      writeFile: temporaryDirectory$ + "/SyllableNuclei.tmp", header$
    elif data$ == "Save as text file" and dataCollectionType$ == "OverWriteData"
      writeFile: txtFile$, header$
    elif data$ == "Table"
      writeFile: temporaryDirectory$ + "/SyllableNuclei.tmp", header$
      endif
  
    if detect_Filled_Pauses
      if data$ == "Praat Info window" and dataCollectionType$ == "OverWriteData"
        appendInfoLine: ", nrFP, tFP(s)"
      elif data$ == "Save all" and dataCollectionType$ == "OverWriteData"
        appendFileLine: txtFile$, ", nrFP, tFP(s)"
        appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ", nrFP, tFP(s)"
      elif data$ == "Save as text file" and dataCollectionType$ == "OverWriteData"
        appendFileLine: txtFile$, ", nrFP, tFP(s)"
      elif data$ == "Table"
        appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ", nrFP, tFP(s)"
        endif
    else
      if data$ == "Praat Info window" and dataCollectionType$ == "OverWriteData"
        appendInfoLine: ""
      elif data$ == "Save all" and dataCollectionType$ == "OverWriteData"
        appendFileLine: txtFile$, ""
        appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ""
      elif data$ == "Save as text file" and dataCollectionType$ == "OverWriteData"
        appendFileLine: txtFile$, ""
      elif data$ == "Table"
        appendFileLine: temporaryDirectory$ + "/SyllableNuclei.tmp", ""
        endif
      endif
    endproc
  
  procedure coda
    if fileSpec$ <> ""
      removeObject: idStr
      endif
    if data$ == "Table" or data$ == "Save all"
      tbl = Read Table from comma-separated file: temporaryDirectory$ + "/SyllableNuclei.tmp"
      deleteFile: temporaryDirectory$ + "/SyllableNuclei.tmp"
      if not keep_Objects
        removeObject: tbl
        endif
      endif
    if nrFiles and keep_Objects
      selectObject: idSnd#, idTG#
      if detect_Filled_Pauses and (data$ == "Table" or data$ == "Save all")
        plusObject: idTbl#
        endif
      endif
    endproc
    appendInfoLine: "Done with the uhm-o-meter."
endif