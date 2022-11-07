# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:13:39 2019

@author: afe02
"""

import numpy as np
import math
import csv
import switchDataForCPU
    
def timesAsArrays(executionTimes, releaseTimes, distLimit):
    etsSummed = [executionTimes[0]]
    releaseTimesResult=[releaseTimes[0]]
    for i in range(1, len(releaseTimes)):
        if releaseTimes[i] - releaseTimesResult[-1] < distLimit:
            etsSummed [-1] += executionTimes[i]
        else:
            etsSummed.append(executionTimes[i])
            releaseTimesResult.append(releaseTimes[i])
        
    etArray = np.zeros(len(etsSummed), dtype=np.int64)
    rtArray = np.zeros(len(etsSummed), dtype=np.int64)
    for i in range(len(etsSummed)):
        etArray[i] = etsSummed[i]
        rtArray[i] = releaseTimesResult[i]
    return etArray, rtArray
    
            
def exportTimingDataToCSV(filename, outputBaseFileName, processName):
    switchWakeupData = switchDataForCPU.getSwitchAndWakeupDataForCPU(filename, '[002]')

    releaseTimeDict, schedulingTimeDict, executionTimeDict, previousProcessList,\
    wakeupInLatencyProcessList, wakeupInExecutionProcessList = \
        switchDataForCPU.getTimeDicts(switchWakeupData, processName)

    executionTimes, releaseTimes = timesAsArrays(executionTimeDict['all'], releaseTimeDict["all"], 3e5)
    
    stateTimesFileName = outputBaseFileName + '_full.csv'
    with open(stateTimesFileName, 'w', newline='') as f:
        timeswriter = csv.writer(f, delimiter=',')
        timeswriter.writerow(("executionTime","releaseTime"))
        for j in range(len(executionTimes)):
            timeswriter.writerow((executionTimes[j],releaseTimes[j]))
    executionTimes = executionTimes[2000:len(executionTimes)-1]
    releaseTimes = releaseTimes[2000:len(releaseTimes)-1]
    stateTimesFileName = outputBaseFileName + '.csv'
    with open(stateTimesFileName, 'w', newline='') as f:
        timeswriter = csv.writer(f, delimiter=',')
        timeswriter.writerow(("executionTime","releaseTime"))
        for j in range(len(executionTimes)):
            timeswriter.writerow((executionTimes[j],releaseTimes[j]))
    stateTimesFileName = outputBaseFileName + '_et_log.csv'
    with open(stateTimesFileName, 'w', newline='') as f:
        timeswriter = csv.writer(f, delimiter=',')
        for j in range(len(executionTimes)):
            timeswriter.writerow((int(1000*math.log(executionTimes[j])),))
    

reportFileName = '../data/traces_reports_csv/controlReportDeadline'
executionTimesFileName = \
    '../data/traces_reports_csv/control_time'
exportTimingDataToCSV(reportFileName, executionTimesFileName, 'ripControl')

        
