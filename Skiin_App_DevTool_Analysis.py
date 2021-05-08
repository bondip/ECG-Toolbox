import os
import csv
import zlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import plotly.express as px
from plotly.offline import plot
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def pako_inflate(data):
    decompress = zlib.decompressobj(15)
    decompressed_data = decompress.decompress(data)
    decompressed_data += decompress.flush()
    return decompressed_data

def SNT_SigQuality(data, fs, row_num, col_num):
    
    snt_raw = []
    
    data = np.array(data)
    
    segment_size = 10       #this is the segment size in seconds
    
    num_segs = len(data)/(320*segment_size)
    
    for i in range(int(np.trunc(num_segs)-1)):
        data_seg = data[int(i*(len(data)/num_segs)):int((i+1)*(len(data)/num_segs))]
        fig.add_trace(go.Scatter(x=[int(i*(len(data)/num_segs)),int((i+1)*(len(data)/num_segs))],y=[np.max(data_seg),np.max(data_seg)], mode='lines',line=dict(color="magenta")), row=row_num, col=col_num)
        fig.add_trace(go.Scatter(x=[int(i*(len(data)/num_segs)),int((i+1)*(len(data)/num_segs))],y=[np.min(data_seg),np.min(data_seg)], mode='lines',line=dict(color="magenta")), row=row_num, col=col_num)
        # plt.plot([int(i*(len(data)/num_segs)),int((i+1)*(len(data)/num_segs))],[np.max(data_seg),np.max(data_seg)],'g')
        # plt.plot([int(i*(len(data)/num_segs)),int((i+1)*(len(data)/num_segs))],[np.min(data_seg),np.min(data_seg)],'g')
        snt_raw.append(np.max(data_seg) - np.min(data_seg))

    snt_avg = np.mean(snt_raw)
    
    return snt_raw, snt_avg

def count_sort_skiin_app_ecg(ecg_data):
    
    #### Note: This Counting Sort Method is not FIFO.
    
    ts = []
        
    for row in ecg_data:
        ts.append(int(row[30]))
        
    
    # Counting Sort is n+k complexity... so big numbers like unix timestamps
    # make it slow. So I make them smaller by subtrating the smallest timestamp
    start = min(ts)
    ts = [x - start for x in ts]
    
    # The output character array that will have sorted arr 
    output = [0 for i in range(len(ts))]
    data_output = [0 for i in range(len(ts))]
    
    # Create a count array to store count of inidividul characters and
    # initialize count array as 0 
    count = [0 for i in range(max(ts)+1)]
    
    # For storing the resulting answer since the string is immutable 
    sorted_ts = [0 for _ in ts] 
    sorted_data = [0 for _ in ts]
    
    # Store count of each character 
    for i in ts: 
        count[i] += 1

    # Change count[i] so that count[i] now contains actual position of this
    # character in output array 
    for i in range(1,len(count),1): 
        count[i] += count[i-1]
      
    # Build the output character array 
    for i in range(len(ts)): 
        output[count[ts[i]]-1] = ts[i]
        data_output[count[ts[i]]-1] = ecg_data[i]
        count[ts[i]] -= 1
      
    # Copy the output array to arr, so that arr now contains sorted characters 
    for i in range(len(ts)): 
        sorted_ts[i] = output[i]
        sorted_data[i] = data_output[i]
    
#    plt.figure(100)
#    plt.plot(sorted_ts)

    return sorted_data
        

def parse_ECG_data(file):
    
    data = dict([
        ('ECG_Ch1',[]),
        ('ECG_Ch2',[]),
        ('ECG_Ch3',[]),
        ('ECG_Ch1_Timestamp',[]),
        ('ECG_Ch2_Timestamp',[]),
        ('ECG_Ch3_Timestamp',[]),
        ('ECG_Ch1_LOD',[]),
        ('ECG_Ch1_LOD_Pos',[]),
        ('ECG_Ch1_LOD_Neg',[]),
        ('ECG_Ch2_LOD',[]),
        ('ECG_Ch2_LOD_Pos',[]),
        ('ECG_Ch2_LOD_Neg',[]),
        ('ECG_Ch3_LOD',[]),
        ('ECG_Ch3_LOD_Pos',[]),
        ('ECG_Ch3_LOD_Neg',[])
        ])
    
    unsorted_ch1_data = []
    unsorted_ch2_data = []
    unsorted_ch3_data = []
    
    ecg1_timestamp = []
    ecg1 = []
    lod_ch1 = []
    lod_ch1_pos = []
    lod_ch1_neg = []
    
    ecg2_timestamp = []
    ecg2 = []
    lod_ch2 = []
    lod_ch2_pos = []
    lod_ch2_neg = []
    
    ecg3_timestamp = []
    ecg3 = []
    lod_ch3 = []
    lod_ch3_pos = []
    lod_ch3_neg = []
    
    with open(file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if row[0] == '26656':
                unsorted_ch1_data.append(row)
            if row[0] == '26657':
                unsorted_ch2_data.append(row)
            if row[0] == '26658':
                unsorted_ch3_data.append(row)
    
    sorted_ch1_data = count_sort_skiin_app_ecg(unsorted_ch1_data)
    sorted_ch2_data = count_sort_skiin_app_ecg(unsorted_ch2_data)
    # sorted_ch3_data = count_sort_skiin_app_ecg(unsorted_ch3_data)
    
    for row in sorted_ch1_data:
        try:
            ecg1.extend([int(x) for x in row[1:25]])
            ecg1_timestamp.append(int(row[30])/1000)
            lod_ch1.append(int(row[26]))
            lod_ch1_neg.append(int(row[27]))
            lod_ch1_pos.append(int(row[28]))
        except:
            continue
    
    for row in sorted_ch2_data:
        try:
            ecg2.extend([int(x) for x in row[1:25]])
            ecg2_timestamp.append(int(row[30])/1000)
            lod_ch2.append(int(row[26]))
            lod_ch2_neg.append(int(row[27]))
            lod_ch2_pos.append(int(row[28]))
        except:
            continue
        
    # for row in sorted_ch3_data:
    #     if row[0] == '26658':
    #         try:
    #             ecg3.extend([int(x) for x in row[1:25]])
    #             ecg3_timestamp.append(int(row[25]))
    #             lod_ch3.append(int(row[26]))
    #             lod_ch3_neg.append(int(row[27]))
    #             lod_ch3_pos.append(int(row[28]))
    #         except:
    #             continue

    ecg1 = medfilt (ecg1)
    ecg2 = medfilt (ecg2)
    
    if len(ecg1) > 0:
        ecg1 = codes_to_volts(ecg1)
    if len(ecg2) > 0:
        ecg2 = codes_to_volts(ecg2)
    # if len(ecg3) > 0:
    #     ecg3 = codes_to_volts(ecg3)
    
    
    data['ECG_Ch1'] = ecg1 if len(ecg1) > 0 else 0
    data['ECG_Ch2'] = ecg2 if len(ecg2) > 0 else 0
    # data['ECG_Ch3'] = ecg3 if len(ecg3) > 0 else 0
    data['ECG_Ch1_Timestamp'] = ecg1_timestamp if len(ecg1_timestamp) > 0 else 0
    data['ECG_Ch2_Timestamp'] = ecg2_timestamp if len(ecg2_timestamp) > 0 else 0
    # data['ECG_Ch3_Timestamp'] = ecg3_timestamp if len(ecg3_timestamp) > 0 else 0
    data['ECG_Ch1_LOD'] = lod_ch1 if len(lod_ch1) > 0 else 0
    data['ECG_Ch2_LOD'] = lod_ch2 if len(lod_ch2) > 0 else 0
    # data['ECG_Ch3_LOD'] = lod_ch3 if len(lod_ch3) > 0 else 0
    data['ECG_Ch1_LOD_Pos'] = lod_ch1_pos if len(lod_ch1_pos) > 0 else 0
    data['ECG_Ch2_LOD_Pos'] = lod_ch2_pos if len(lod_ch2_pos) > 0 else 0
    # data['ECG_Ch3_LOD_Pos'] = lod_ch3_pos if len(lod_ch3_pos) > 0 else 0
    data['ECG_Ch1_LOD_Neg'] = lod_ch1_neg if len(lod_ch1_neg) > 0 else 0
    data['ECG_Ch2_LOD_Neg'] = lod_ch2_neg if len(lod_ch2_neg) > 0 else 0
    # data['ECG_Ch3_LOD_Neg'] = lod_ch3_neg if len(lod_ch3_neg) > 0 else 0
    
    return data


def parse_acc_data(file):
    
    data = dict([
        ('accelx',[]),
        ('accely',[]),
        ('accelz',[]),
        ('accel_Timestamp',[])
        ])
        
    accelx = []
    accely = []
    accelz = []
    accel_timestamp = []
    
    with open(file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if row[0] == '26659':
                try:
                    accelx.extend([int(x) for x in row[1:37:3]])
                    accely.extend([int(x) for x in row[2:37:3]])
                    accelz.extend([int(x) for x in row[3:37:3]])
                    accel_timestamp.append(int(row[38])/1000)
                except:
                    continue
    
    data['accelx'] = accelx
    data['accely'] = accely
    data['accelz'] = accelz
    data['accel_Timestamp'] = accel_timestamp
    
    return data


def parse_temp_data(file):
    
    data = dict([
        ('heat_flux',[]),
        ('pod_temp',[]),
        ('core_body_temp',[]),
        ('corr_heat_flux',[]),
        ('corr_pod_temp',[]),
        ('temp_timestamp',[])
        ])
    
    heat_flux = []
    pod_temp = []
    core_body_temp = []
    corr_pod_temp = []
    corr_heat_flux = []
    temp_timestamp = []
    
    with open(file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if row[0] == '26660':
                try:
                    heat_flux.append(int(row[1]))
                    pod_temp.append(float(row[2])/100)
                    core_body_temp.append(float(row[3])/100)
                    corr_heat_flux.append(float(row[4]))
                    corr_pod_temp.append(float(row[5]))
                    temp_timestamp.append(int(row[6]))
                except:
                    continue
    
    data['heat_flux'] = heat_flux
    data['pod_temp'] = pod_temp
    data['core_body_temp'] = core_body_temp
    data['corr_heat_flux'] = corr_heat_flux
    data['corr_pod_temp'] = corr_pod_temp
    data['temp_timestamp'] = temp_timestamp
    
    return data

def medfilt (sig):
    win = 10
    starts = len(sig)//win
    medFilt_sig = []

    for i in range(starts):
        start = i*10
        seg = sig[start:start+win]
        med = np.median(seg)
        seg = [x-med for x in seg]
        medFilt_sig[start:start+win] = seg
    
    return medFilt_sig

def codes_to_volts(sig):
    
    #This FCN returns uV for the underwear module EVT3.5 and earlier
#    sig_v = [x*0.3814697 for x in sig]
    
    #This FCN returns uV for the underwear module EVT3.6 with the ADS IC
    sig = np.array(sig)
    sig -= sig[0]
    sig_v = sig/9.3333
    #9333 codes/mV or 9.333 codes/uV or 107.15 nV/code
    #sig_v = sig
    
    return sig_v

##############################################################################

ecg_fs = 320


if __name__ == '__main__':
    
    path = r'G:\.shortcut-targets-by-id\0BwbMF_G6GNgNTW44d2lTUVFsZms\SKIIN\Product Development\1. Hardware Development\Data Collection\DVT1 Testing\BMX vs BMI Testing\SNT BMI Pod 04'
    
    ecg_question = input("Do you want to analyze ECG Data? (y=Yes n=No): ")
    acc_question = input("Do you want to analyze Accelerometer Data? (y=Yes n=No): ")
    temp_question = input("Do you want to analyze CBT Data? (y=Yes n=No): ")
    
    data={}

    files = os.listdir(path) #get all the files inside of the folder
    all_files = [name for name in files if name.endswith('.csv')]
    for name in all_files:
        if name.find('ecg') != -1 and ecg_question == 'y':
            ecg_file = name
            data['ecg_data'] = parse_ECG_data(path + '/' + ecg_file)
        if name.find('acc') != -1 and acc_question == 'y':
            acc_file = name
            data['accel_data'] = parse_acc_data(path + '/' + acc_file)
        if name.find('temp') != -1 and temp_question == 'y':
            temp_file = name
            data['temp_data'] = parse_temp_data(path + '/' + temp_file)
        
    data['ecg_data']['ECG_Ch1'] = data['ecg_data']['ECG_Ch1'][:]
    data['ecg_data']['ECG_Ch2'] = data['ecg_data']['ECG_Ch2'][:]
    # data['ecg_data']['ECG_Ch3'] = data['ecg_data']['ECG_Ch3'][25000:]
    
    plt.figure(1)
    plt.plot(data['ecg_data']['ECG_Ch1'])
    plt.title('ECG Ch1 SNT')
    plt.ylabel('ECG Amplitude [uV]')
    plt.xlabel('ECG Sample Number')
    plt.grid()
    
    plt.figure(2)
    plt.plot(data['ecg_data']['ECG_Ch2'])
    plt.title('ECG Ch2 SNT')
    plt.ylabel('ECG Amplitude [uV]')
    plt.xlabel('ECG Sample Number')
    plt.grid()
    
    
    fig = make_subplots(
        rows=2, cols=1,
        vertical_spacing=0.1,
        horizontal_spacing=0.07,
        subplot_titles=["ECG Ch1 SNT", "ECG_Ch2 SNT"])
    
    fig.add_trace(go.Scatter(y=data['ecg_data']['ECG_Ch1'], name='ECG Ch1', mode='lines',legendgroup='Raw'), row=1, col=1)
    fig.add_trace(go.Scatter(y=data['ecg_data']['ECG_Ch2'], name='ECG Ch2', mode='lines',legendgroup='Raw'), row=2, col=1)

    
    ch1_avg, ch1_amps = SNT_SigQuality(data['ecg_data']['ECG_Ch1'], ecg_fs, 1,1)
    ch2_avg, ch2_amps = SNT_SigQuality(data['ecg_data']['ECG_Ch2'], ecg_fs, 2,1)
    
    
    fig.update_xaxes(title_text="Samples", row=1, col=1)
    fig.update_yaxes(title_text="ECG Amplitude [uV]", row=1, col=1)
    fig['layout']['yaxis1'].update(matches="y1")
    plot(fig, filename= "SNT.html")
    # accdata = r'C:\Users\Muammar\Desktop\acc-50981400165-1609376486721-4.3.4' #'C:/Users/Muammar/Desktop/New folder (2)/'
    # fs_acc = 25
    
    # _, AccX, AccY, AccZ, Acc_Timestamp = parse_acc_data(accdata)
    
    # tAcc = np.arange(0,(len(AccX)/fs_acc), 1/fs_acc)
    
    # fig2=plt.figure(2)
    # ax1 = plt.subplot(311)
    # ax1.plot(tAcc,AccX)
    # plt.grid()   
    # #plt.xlabel('Time (s)', fontsize=16)
    # plt.ylabel('Acc - X', fontsize=16)
    
    # ax2=plt.subplot(312, sharex=ax1)
    # ax2.plot(tAcc,AccY)
    # plt.grid()   
    # #plt.xlabel('Time (s)', fontsize=16)
    # plt.ylabel('Acc - Y', fontsize=16)

    # ax3=plt.subplot(313, sharex=ax1)
    # ax3.plot(tAcc,AccZ)
    # plt.grid()   
    # plt.xlabel('Time (s)', fontsize=16)
    # plt.ylabel('Acc - Z', fontsize=16)
    