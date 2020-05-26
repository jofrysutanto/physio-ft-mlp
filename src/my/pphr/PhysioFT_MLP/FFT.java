package my.pphr.PhysioFT_MLP;
/*
 Copyright (C) 2001, 2006 by Simon Dixon 
  This program is free software; you can redistribute it and/or modify 
 it under the terms of the GNU General Public License as published by 
 the Free Software Foundation; either version 2 of the License, or 
 (at your option) any later version. 
  This program is distributed in the hope that it will be useful, 
 but WITHOUT ANY WARRANTY; without even the implied warranty of 
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 GNU General Public License for more details. 
  You should have received a copy of the GNU General Public License along 
 with this program (the file gpl.txt); if not, download it from 
 http://www.gnu.org/licenses/gpl.txt or write to the 
 Free Software Foundation, Inc., 
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 
 This file was modified in 2012 by Thomas Friedel
   adapted to PPG analysis by Lee Seldon 2017 */ 
  
//import java.util.Arrays; 
//import java.awt.Frame; // to open PPG file
//import java.awt.FileDialog;
//import java.io.File;
//import java.io.FileReader;
//import java.io.BufferedReader;
//import java.io.PrintWriter;
//import java.io.FileNotFoundException;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.Scanner;

/** Class for computing a windowed fast Fourier transform.
 *  Implements some of the window functions for the STFT from 
 *  Harris (1978), Proc. IEEE, 66, 1, 51-83. 
 */ 
public class FFT { 
 
     // used in FFT.fft(float[], float[], int) to specify forward or reverse (inverse) transform
     public static final int FORWARD = -1; 
     public static final int REVERSE = 1; 
     public static final float twoPI = (float)(2. * Math.PI); 
     // FFT ARRAYS
     public static float[] fftRe = null;      // FFT real
     public static float[] fftIm = null;      // FFT imaginary
	 public static float[] fftPow = null;      // FFT power
     public static float[] fftReFilt = null;  // filtered FFT real 
     public static float[] fftImFilt = null;  // filtered FFT imaginary
     public static float[] reRecon = null;    // reverse FFT real
     public static float[] imRecon = null;    // reverse FFT imaginary
     // input data arrays for PhysioCam and MIMIC
     public static float[] inputTime, inputTimeInterval;   // sample times
     public static float[] inputYampl;                 // Y amplitude from YUV
     public static float yAmplMin, yAmplMax, yAmplRange, yAmplAve, yAmplStD;
     // input arrays only for PhysioCam data
     public static float[] inputUampl, inputVampl; // YUV video
	 //public static float uAmplMin, uAmplMax, uAmplRange, uAmplAve, uAmplStD;
     //public static float vAmplMin, vAmplMax, vAmplRange, vAmplAve, vAmplStD;
	 public static int[] inputHeartBeats;
     
     public static int nLgt = 8192; // FFT.nLgt = FFT array size max
	 public static int startI = 0;
     public static float totalTime, hrFreqHz; // hrFreq = HR in Hz
     // fftRe[hrPeakLoc] = element of HR peak. Same for harmonics of HR
     public static int hrPeakLoc = 0, hrHm1Loc = 0, hrHm2Loc = 0;     
     public static int totalHeartBeats = 0;
	 
	 /* input Meta-data for each sample, if known
        mtD[0] = HR as BPM
        mtD[1] = systolic BP
        mtD[2] = diastolic BP
        mtD[3] = BP range = systolic - diastolic
        mtD[4] = input file ref number - on output this comes FIRST
        mtD[5] = offset of sample within input file - ignored on output
			   = numPixels used for YUV capture of phone PPG */
     public static int mtD[] = new int[6];
     // BP variables - frequency-related
     // Ratios provide sort of normalization
     public static float HR_Hz; // the HR in Hz
     public static float HR_Sec; // the peak-peak time in sec = 1/HR
     public static float Hm1_Hz; // the Harmonic1 in Hz
     public static float Hm1Rel; // = Hm1_Hz/HR_Hz
     public static float Hm1_Sec; // the Hm1 as sec = 1/Hm1_Hz
     public static float Hm2_Hz; // the Harmonic2 in Hz
     public static float Hm2Rel; // = Hm2_Hz/HR_Hz
     public static float Hm2_Sec; // Hm2 as sec = 1/Hm2_Hz
     public static float Resp_Hz; // the Respiration rate in Hz
     public static float RespRel; // = Resp_Hz/HR_Hz
     public static float Resp_Sec; // resp as sec - may not be useful
     // BP variables - Power of peaks and waves
     // HR and Harmonics divide power spectrum into 3 parts: HR-Hm1, Hm1-Hm2, Hm2-end
     public static float pow01HR; // power of single HR peak
     public static float pow01HRa; // ave power of HR peak, including sides
     public static float pow02Hm1; // power of Hm1 peak
     public static float pow03Hm2; // power of Hm2 peak
     public static float pow04Resp; // power of Resp peak, including sides
     public static float pow05HR1; // power from HR-Hm1 peak
     public static float pow06HR2; // power from HR-Hm2 peak
     public static float pow07HRe; // power from HR peak -end
     public static float pow08Hm12; // power from Hm1-Hm2 peak
     public static float pow09Hm1e; // power from Hm1-end
     public static float pow10Hm2e; // power from Hm2-end
     // POWER RATIOS
     //  of peaks - Rp
     public static float powRp01Hm1_HR; // = pow02Hm1/pow01HR
     public static float powRp02Hm2_HR; // = pow03Hm2/pow01HR
     public static float powRp03Resp_HR; // = pow04Resp/pow01HR
     //  of areas or parts of waves - Ra
     public static float powRa01HR1_HR2; // = pow05HR1/pow06HR2
     public static float powRa02HR1_HRe; // = pow05HR1/pow07HRe
     public static float powRa03Hm12_HR1; // = pow08Hm12/pow05HR1
     public static float powRa04Hm12_HR2; // = pow08Hm12/pow06HR2
     public static float powRa05Hm12_HRe; // = pow08Hm12/pow07HRe
     public static float powRa06Hm12_Hm1e; // = pow08Hm12/pow09Hm1e
     public static float powRa07Hm1e_HR1; // = pow09Hm1e/pow05HR1
     public static float powRa08Hm1e_HRe; // = pow09Hm1e/pow07HRe
     public static float powRa09Hm2e_HR1; // = pow10Hm2e/pow05HR1
     public static float powRa10Hm2e_HR2; // = pow10Hm2e/pow06HR2
     public static float powRa11Hm2e_HRe; // = pow10Hm2e/pow07HRe
     public static float powRa12Hm2e_Hm1e; // = pow10Hm2e/pow09Hm1e
     // others
     public static float spO2, HRV; // SpO2, HR Variability
     public static float waveAmpMx, waveAmpAve; // PPG ampl
     // Label string. 23 more important vars are first
     public static String outHeadr = "HR_Hz,HR_Sec,Hm1_Hz,Hm1Rel,Hm1_Sec,Hm2_Hz,Hm2Rel,Hm2_Sec,"
       + "Resp_Hz,RespRel,Resp_Sec,pow01HR,pow01HRa,pow02Hm1,pow03Hm2,pow04Resp,"
       + "powRp01Hm1_HR,powRp02Hm2_HR,powRp03Resp_HR,"
       + "powRa03Hm12_HR1,powRa05Hm12_HRe,powRa06Hm12_Hm1e,powRa07Hm1e_HR1,"
       + "powRa09Hm2e_HR1,powRa10Hm2e_HR2,powRa11Hm2e_HRe,HRV,waveAmpMx,waveAmpAve,"
	   + "pow05HR1,pow06HR2,pow07HRe,pow08Hm12,pow09Hm1e,pow10Hm2e," 
	   + "powRa01HR1_HR2,powRa02HR1_HRe,powRa04Hm12_HR2,powRa08Hm1e_HRe,powRa12Hm2e_Hm1e," 
	   + "spO2,HRbpm,SBP,DBP,BPrg,recNum,numPix";
	 
     // input file with PPG
	 //  public static Frame fileDialogFrame;
	 //  public static FileDialog fDialog;
     public static String inputLine, openFilename;
     public static String[] lineParts;
     // output file with FFT
     //public static String outFilename, outFilenamePrefix;
     //public static PrintWriter fftFileWriter;
     // filter choices
	 // public static Scanner scn;
     public static int filtLowFreq = 0;  // 1 -> set FFT freq < HR = 0
     public static int filtNoise = 0;  // 1 -> set FFT freq < ave = 0
     public static int filtShiftHR = 0;  // if >0, shift HR peaks -> filtShiftHR
     public static int filtDelHR = 0;  // if > 0, del HR peaks between shiftHR-delHR
 
 /** The FFT method. Calculation is inline, for complex data stored
  *  in 2 separate arrays. Length of input data must be a power of two. 
  *  @param re        the real part of the complex input and output data 
  *  @param im        the imaginary part of the complex input and output data 
  *  @param direction -direction of Fourier transform (FORWARD or REVERSE) 
  *  @throws IllegalArgumentException if the length of the input data is 
  *     not a power of 2 
  */ 
 public static void fft(float re[], float im[], int direction) { 
    int n = re.length;
	// if n not 2^x, get it from FFT.nLgt, which was set in getArrays...()
	//while ((FFT.nLgt = (FFT.nLgt >>> 1)) > n);
	if (n > FFT.nLgt) n = FFT.nLgt;
    int bits = (int)Math.rint(Math.log(n) / Math.log(2)); 
    //if (n != (1 << bits)) 
    //    throw new IllegalArgumentException("FFT data must be power of 2"); 
    int localN; 
    int j = 0, kn, nby2, id;
    float Wjk_r, Wjk_i, Wj_r, Wj_i, theta, temp, tempr, tempi, wtemp;
    for (int i = 0; i < n-1; i++) { 
        if (i < j) { 
            temp = re[j];  re[j] = re[i];  re[i] = temp; 
            temp = im[j];  im[j] = im[i];  im[i] = temp; 
        } 
        kn = n / 2; 
        while ((kn >= 1) &&  (kn - 1 < j)) { 
            j = j - kn; 
            kn = kn / 2; 
        } 
        j = j + kn; 
    } // End for i = 0 to n-1 
    for(int m = 1; m <= bits; m++) { 
        localN = 1 << m; 
        Wjk_r = 1; 
        Wjk_i = 0; 
        theta = twoPI / localN; 
        Wj_r = (float)Math.cos(theta); 
        Wj_i = direction * (float)Math.sin(theta); 
        nby2 = localN / 2; 
        for (j = 0; j < nby2; j++) { 
            for (int k = j; k < n; k += localN) { 
                id = k + nby2; 
                tempr = Wjk_r * re[id] - Wjk_i * im[id]; 
                tempi = Wjk_r * im[id] + Wjk_i * re[id]; 
                re[id] = re[k] - tempr; 
                im[id] = im[k] - tempi; 
                re[k] += tempr; 
                im[k] += tempi; 
            } // end for k = j to n 
            wtemp = Wjk_r; 
            Wjk_r = Wj_r * Wjk_r  - Wj_i * Wjk_i; 
            Wjk_i = Wj_r * Wjk_i  + Wj_i * wtemp; 
        } // end for j = 0 to nby2 
    } // end for m = 1 to bits 
 } // END fft() 
  
 /** public static void doFFT(boolean doFilt, boolean doRecon)
  *  Creates FFT real & imag arrays, copies input to them, runs fft()
     Copies output into FFT filter real & imag arrays, cuts freq < heart rate
     Copies output into FFT recon arrays, runs reverse fft()
        to reconstruct filtered input
     @param doFilt - do the filtering
     @param doRecon - do the reconstruction 
     ALSO note last vals in writeFFTline(..5, 2) says write 5 metadata
        and write only 1/2 of inRe to output
 */
 public static void doFFT(boolean doFilt, boolean doRecon) {
    //System.out.println("forward FFT ->fftRe"); 
    if (FFT.fftRe == null) {       // IF FFT arrays do NOT yet exist, create them
        FFT.fftRe = new float[FFT.nLgt]; 
        FFT.fftIm = new float[FFT.nLgt];
		for (int j = 0; j < FFT.nLgt; j++) { 
			FFT.fftRe[j] = FFT.inputYampl[j];  // YUV amplitudes are real part
			FFT.fftIm[j] = FFT.inputTimeInterval[j];  // maybe time or time intervals could be imaginary
		}  
	} // end creating and filling fftRe and fftIm
	// IF FFT arrays ALREADY EXIST, assume they have been filled outside of this method
	
    fft(FFT.fftRe, FFT.fftIm, FORWARD);
    FFT.fftRe[0] = FFT.fftIm[0] = 0.f;  // zero DC offset
	// create and fill FFT power array
    if (FFT.fftPow == null) FFT.fftPow = new float[FFT.nLgt];
    for (int i1 = 0; i1 < FFT.nLgt; i1++) {
        FFT.fftPow[i1] = FFT.fftRe[i1] * FFT.fftRe[i1]
                        + FFT.fftIm[i1] * FFT.fftIm[i1];
    }
    // Find physiol params
    findPhysioInFFT(FFT.fftPow, FFT.fftImFilt);
	
    if (doFilt) {
        //System.out.println("Find HR, filter FFT fftRe->fftReFilt");
        if (FFT.fftReFilt == null) {       // first time create arrays
            FFT.fftReFilt = new float[FFT.nLgt]; 
            FFT.fftImFilt = new float[FFT.nLgt]; }
        for (int j = 0; j < FFT.nLgt; j++) { 
            FFT.fftReFilt[j] = FFT.fftRe[j];
            FFT.fftImFilt[j] = FFT.fftIm[j]; 
        } 
        // Filter FFT. Shift down to 1-3, del HR peaks within 0-4-5, del high-freq
        filterFFT(FFT.fftReFilt, FFT.fftImFilt, FFT.filtLowFreq, FFT.filtNoise, FFT.filtShiftHR, FFT.filtDelHR); 
        // Write output line. First half of fftRe, 5 metadata
//        writeFFTline(FFT.fftReFilt, FFT.fftImFilt, 5, 2);
        
        if (doRecon) {
            //System.out.println("reverse fftReFilt ->reRecon");
            if (FFT.reRecon == null) {       // first time create arrays
                FFT.reRecon = new float[FFT.nLgt]; 
                FFT.imRecon = new float[FFT.nLgt]; }
            for (int j = 0; j < FFT.nLgt; j++) { 
                FFT.reRecon[j] = FFT.fftReFilt[j];
                FFT.imRecon[j] = FFT.fftImFilt[j]; 
            } 
            fft(FFT.reRecon, FFT.imRecon, REVERSE); 
        } // end if doRecon
    } // end if doFilt
 } // End doFFT()
   
 /** findPhysioInFFT(float[] re, float[] im)
  * Search FT real[] and/or imag[] arrays for physiol features.
  * In reality we will use the fftPow[] power instead of real[]
  * Looks for physiological parameters in the PPG or the FFT of PPG.
  *    1) Heart Rate HR is the largest peaks in the FFT between 0.5-3 Hz (=30-180 BPM).
  *       The harmonics of HR might indicate HR variability, or the width of the PPG wave,
  *       or the dichrotic notch, or ...
  *       Number of freq in HR peak cluster indicates HR variability.
  *    2) The FFT peak(s) between 0.05-0.5 Hz might be respiratory rate (3-30 breaths/min).
  *    3) The magnitude of the range in U (red-green component of YUV) might be related to SpO2.
  *       If we assume that the absolute minimum O2 blood level is 70% (venous return), 
  *       and that corresponds to an extreme U value, then 100% O2 blood level should correspond
  *       to the theoretically maximal opposite extreme U value. Not yet implemented 2017-11-21
  * 1) For base HR search from 0.5 Hz (HR 30bpm fft[index]~totalTime/2)to 3 Hz (HR 180 bpm)
       HOW TO DEFINE HR PEAK?
           Look for greatest re[] amplitude between 0.5-3 Hz. Elem re[i1]
           IF re[i1]--re[i1-1]--re[i1-2] ALL large, then HR peak is re[i1-1]
           ELSE IF only re[i1]--re[i1-1] large, then HR peak is the greater of the 2
           ELSE re[i1] is the HR peak
       PPG diagrams indicate that the HR peak itself is ~1/3 the length of a PPG wave,
         so made of frequencies ~ 3*HR.
       Dichrotic notch is ~halfway through a PPG wave,
         so made of frequencies ~ 2*HR.
       Dichrotic notch width is ~1/10 of PPG wave,
         so made of frequencies ~ 10*HR.
       So filter should find frequencies from ~HR - ~ 10*HR. (The upper bound not found in FFT.)
   FFT over totalTime(sec) ->
       re[0] = constant offset ~signal average ampl - delete this in all cases
       re[1] = 1 wave / totalTime = 1/totalTime Hz
       re[2] = 2 wave / totalTime = 2/totalTime Hz 
       ...re[1/2 totalTime] = 1/2 totalTime wave / totalTime = 1/2 Hz....
	   ...re[totalTime] = 1*totalTime waves / totalTime = 1 Hz....
  *  @param re - FFT real array
     @param im - FFT imaginary array
 */
 protected static void findPhysioInFFT(float[] re, float[] im) {
    int lgt = re.length;    // length of re[]
    int lgtHalf = lgt / 2;  // halfway point in re[]
    int iEnd = lgt - 1;     // last index in re[]

    int point5Hz = (int)(FFT.totalTime / 2000.);      // HR 30/min
	int point7Hz = (7 * point5Hz) / 5;      // HR 40/min
    int to3Hz = (int)(3. * FFT.totalTime / 1000.);  // HR 180/min
    int point05Hz = (int)(FFT.totalTime / 20000.);  // breath rate 3/min
    float aPeakHRre, aPeakHm1, aPeakHm2, aPeakBreath; // peak amplitudes
    int iPeakHRre, iPeakHm1, iPeakHm2, iPeakBreath; // peak indices
    
    if (point05Hz < 3) point05Hz = 3;
    
    aPeakHRre = aPeakHm1 = aPeakHm2 = aPeakBreath = 0.f;
    iPeakHRre = iPeakHm1 = iPeakHm2 = iPeakBreath = -1;
    // Look for HR peak in 0.5 or 0.7 -3 Hz
    for (int i = point7Hz; i < to3Hz; i++) {
        if ((Math.abs(re[i]) + Math.abs(re[i-1])) > aPeakHRre) { // find highest peak in range
            iPeakHRre = i;
            aPeakHRre = Math.abs(re[i]) + Math.abs(re[i-1]);
            if (Math.abs(re[i]) < Math.abs(re[i-1])) {
                iPeakHRre = i - 1;
            }
        } 
    } // end for i in range to find HR freq
	if (iPeakHRre < 1) iPeakHRre = 1;
    FFT.hrPeakLoc = iPeakHRre; // save location of peak
    //FFT.HR_Hz = 1000. * ((float)iPeakHRre) / FFT.totalTime;
    // Look for first Harmonic above HR peak 
    for (int i = iPeakHRre + (iPeakHRre/3); i < lgtHalf; i++) {
        if ((Math.abs(re[i]) + Math.abs(re[i-1])) > aPeakHm1) { // find highest peak in range
            iPeakHm1 = i;
            aPeakHm1 = Math.abs(re[i]) + Math.abs(re[i-1]);
            if (Math.abs(re[i]) < Math.abs(re[i-1])) {
                iPeakHm1 = i - 1;
            }
        } 
    } // end for i in range to find iPeakHm1
	if (iPeakHm1 < iPeakHRre) iPeakHm1 = iPeakHRre;
    FFT.hrHm1Loc = iPeakHm1;
    // Look for 2nd Harmonic above first Harmonic
    for (int i = iPeakHm1 + (iPeakHRre/3); i < lgtHalf; i++) {
        if ((Math.abs(re[i]) + Math.abs(re[i-1])) > aPeakHm2) { // find highest peak in range
            iPeakHm2 = i;
            aPeakHm2 = Math.abs(re[i]) + Math.abs(re[i-1]);
            if (Math.abs(re[i]) < Math.abs(re[i-1])) {
                iPeakHm2 = i - 1;
            }
        } 
    } // end for i in range to find aPeakHm2
	if (iPeakHm2 < iPeakHm1) iPeakHm2 = iPeakHm1;
    FFT.hrHm2Loc = iPeakHm2;
    // Look for respiration peak
    for (int i = point05Hz + 1; i < point5Hz; i++) {
        if ((Math.abs(re[i]) + Math.abs(re[i-1])) > aPeakBreath) { // find highest peak in range
            iPeakBreath = i;
            aPeakBreath = Math.abs(re[i]) + Math.abs(re[i-1]);
            if (Math.abs(re[i]) < Math.abs(re[i-1])) {
                iPeakBreath = i - 1;
            }
        } 
    } // end for i in range to find aPeakBreath
    
    // CALC all the Physio Var
    FFT.HR_Hz = 1000.f * ((float)iPeakHRre) / FFT.totalTime; // the HR in Hz
    FFT.HR_Sec = 1.f / FFT.HR_Hz; // the peak-peak time in sec = 1/HR
    FFT.Hm1_Hz = 1000.f * ((float)iPeakHm1) / FFT.totalTime; // the Harmonic1 in Hz
    FFT.Hm1Rel = FFT.Hm1_Hz / FFT.HR_Hz; // = Hm1_Hz/HR_Hz
    FFT.Hm1_Sec = 1.f / FFT.Hm1_Hz; // the Hm1 as sec = 1/Hm1_Hz
    FFT.Hm2_Hz = 1000.f * ((float)iPeakHm2) / FFT.totalTime; // the Harmonic2 in Hz
    FFT.Hm2Rel = FFT.Hm2_Hz / FFT.HR_Hz; // = Hm2_Hz/HR_Hz
    FFT.Hm2_Sec = 1.f / FFT.Hm2_Hz; // Hm2 as sec = 1/Hm2_Hz
    FFT.Resp_Hz = 1000.f * (iPeakBreath) / FFT.totalTime; // the Respiration rate in Hz
    FFT.RespRel = FFT.Resp_Hz / FFT.HR_Hz; // = Resp_Hz/HR_Hz
    FFT.Resp_Sec = 1.f / FFT.Resp_Hz; // resp as sec - may not be useful
    // BP variables - Power of peaks and waves
    // HR and Harmonics divide power spectrum into 3 parts: HR-Hm1, Hm1-Hm2, Hm2-end
    FFT.pow01HR = FFT.fftPow[iPeakHRre]; // power of HR peak
    FFT.pow01HRa = (FFT.fftPow[iPeakHRre] + FFT.fftPow[iPeakHRre -1]
                + FFT.fftPow[iPeakHRre +1]) / 3.f; // AVE power of HR peak, including sides
    FFT.pow02Hm1 = FFT.fftPow[iPeakHm1]; // power of Hm1 peak
    FFT.pow03Hm2 = FFT.fftPow[iPeakHm2]; // power of Hm2 peak
    FFT.pow04Resp = FFT.fftPow[iPeakBreath]; // power of Resp peak, including sides
    FFT.pow05HR1 = 0.f; // power from HR-Hm1 peak
    for (int i2 = iPeakHRre; i2 < iPeakHm1; i2++) {
        FFT.pow05HR1 += FFT.fftPow[i2]; }
    FFT.pow08Hm12 = 0.f; // power from Hm1-Hm2 peak
    for (int i2 = iPeakHm1; i2 < iPeakHm2; i2++) {
        FFT.pow08Hm12 += FFT.fftPow[i2]; }
    FFT.pow10Hm2e = 0.f; // power from Hm2-end
    for (int i2 = iPeakHm2; i2 < lgtHalf; i2++) {
        FFT.pow10Hm2e += FFT.fftPow[i2]; }
    FFT.pow06HR2 = FFT.pow05HR1 + FFT.pow08Hm12; // power from HR-Hm2 peak
    FFT.pow07HRe = FFT.pow06HR2 + FFT.pow10Hm2e; // power from HR peak -end
    FFT.pow09Hm1e = FFT.pow08Hm12 + FFT.pow10Hm2e; // power from Hm1-end
    
    // POWER RATIOS
    //  of peaks - Rp
    FFT.powRp01Hm1_HR = FFT.pow02Hm1/FFT.pow01HR; // = pow02Hm1/pow01HR
    FFT.powRp02Hm2_HR = FFT.pow03Hm2/FFT.pow01HR; // = pow03Hm2/pow01HR
    FFT.powRp03Resp_HR = FFT.pow04Resp/FFT.pow01HR; // = pow04Resp/pow01HR
    //  of areas or parts of waves - Ra
    FFT.powRa01HR1_HR2 = FFT.pow05HR1/FFT.pow06HR2; // = pow05HR1/pow06HR2
    FFT.powRa02HR1_HRe = FFT.pow05HR1/FFT.pow07HRe; // = pow05HR1/pow07HRe
    FFT.powRa03Hm12_HR1 = FFT.pow08Hm12/FFT.pow05HR1; // = pow08Hm12/pow05HR1
    FFT.powRa04Hm12_HR2 = FFT.pow08Hm12/FFT.pow06HR2; // = pow08Hm12/pow06HR2
    FFT.powRa05Hm12_HRe = FFT.pow08Hm12/FFT.pow07HRe; // = pow08Hm12/pow07HRe
    FFT.powRa06Hm12_Hm1e = FFT.pow08Hm12/FFT.pow09Hm1e; // = pow08Hm12/pow09Hm1e
    FFT.powRa07Hm1e_HR1 = FFT.pow09Hm1e/FFT.pow05HR1; // = pow09Hm1e/pow05HR1
    FFT.powRa08Hm1e_HRe = FFT.pow09Hm1e/FFT.pow07HRe; // = pow09Hm1e/pow07HRe
    FFT.powRa09Hm2e_HR1 = FFT.pow10Hm2e/FFT.pow05HR1; // = pow10Hm2e/pow05HR1
    FFT.powRa10Hm2e_HR2 = FFT.pow10Hm2e/FFT.pow06HR2; // = pow10Hm2e/pow06HR2
    FFT.powRa11Hm2e_HRe = FFT.pow10Hm2e/FFT.pow07HRe; // = pow10Hm2e/pow07HRe
    FFT.powRa12Hm2e_Hm1e = FFT.pow10Hm2e/FFT.pow09Hm1e; // = pow10Hm2e/pow09Hm1e
    
    // SpO2 is calculated in PhysioCam.PreviewCallback() 
    if ((FFT.spO2 < 0.3f) || (FFT.spO2 > 1.f)) FFT.spO2 = 0.99f;
    // estimate HRV as ratio of powers ((HR peak freq - 1) + (HR peak freq + 1)) / 2*(HR peak freq)
	// e.g.  at HR freq power 0 4 0 => (0 + 0)/(2*4) = 0 means HR has single frequency
	// e.g.  at HR freq power 4 4 4 => (4 + 4)/(2*4) = 1 means HR has 3 equal frequencies
	//float hrvPowVariabil = ( ((re[iPeakHRre - 1] * re[iPeakHRre - 1]) + (im[iPeakHRre - 1] * im[iPeakHRre - 1]))
	//		+ ((re[iPeakHRre + 1] * re[iPeakHRre + 1]) + (im[iPeakHRre + 1] * im[iPeakHRre + 1])) )
	//		/ (2. * ((re[iPeakHRre] * re[iPeakHRre]) + (im[iPeakHRre] * im[iPeakHRre])) );
    FFT.HRV = (float)((Math.abs((double)re[iPeakHRre - 1]) + Math.abs((double)re[iPeakHRre + 1])) / (2. * Math.abs((double)re[iPeakHRre])));
	
    // ppg Y ampl - see readPhysio or readMIMIC methods to set it.
    
    // store HR as first item in metadata mtD[]
    FFT.mtD[0] = ((int)((60. * FFT.HR_Hz)));  // + 2.5)/5.)) * 5;
    
 } // End protected static void findPhysioInFFT(float[] re, float[] im)

 /** filterFFT(float[] re, float[] im, int filtLow, int filtNoise, int shiftHR, int delHR)
   PPG diagrams indicate that the HR peak itself is ~1/3 the length of a PPG wave,
     so made of frequencies ~ 3*HR.
   Dichrotic notch is ~halfway through a PPG wave,
     so made of frequencies ~ 2*HR.
   Dichrotic notch width is ~1/10 of PPG wave,
     so made of frequencies ~ 10*HR.
   So filter should find frequencies from ~HR - ~ 10*HR.
   FFT over totalTime(sec) ->
       re[0] = constant offset ~signal average ampl - delete this in all cases
       re[1] = 1 wave / totalTime = 1/totalTime Hz
       re[2] = 2 wave / totalTime = 2/totalTime Hz 
       ...re[1/2 totalTime] = 1/2 totalTime wave / totalTime = 1/2 Hz....
  *  @param re - FFT real array
     @param im - FFT imaginary array
     @param filtLow - filter low freq FFT, up to HR = 0. 0 - no filt
     @param filtNoise - filter noise (ampl < average/filtNoise) = 0. 0 - no filt
     @param shiftHR - shift FFT down so HR freq at 'shiftHR' element. 0 - no shift
     @param delHR - delete the large-amplitude HR peaks in FFT. 0 - no del
 */
 public static void filterFFT(float[] re, float[] im, int filtLow, int filtNoise, int shiftHR, int delHR) {
    int lgt = re.length;    // length of re[]
    int lgtHalf = lgt / 2;  // halfway point in re[]
    int iEnd = lgt - 1;     // last index in re[]
    int lowCutoff, midNumCut;
    float aveHRReLow = 0.f, aveHRReHigh = 0.f, aveHRImLow = 0.f, aveHRImHigh = 0.f;
    float aveRe;  // average of non-zero re[] elements
    int aveReN;    // number of non-zero re[] elements
    // FILTER NOISE FREQ, those below aveRe / noise cutoff level
    aveRe = 0.f;
    aveReN = 0;
    for (int i = 0; i < lgt; i++) {
        if (i == 0) re[0] = im[0] = 0.f; // del DC offset
        if (re[i] != 0.f) {        // find ave of re[]
            aveRe += Math.abs(re[i]);
            aveReN ++;
        }
    } // end for i loop thru arrays to find averages and HR freq
    if (aveReN > 0) aveRe /= (float)aveReN;
    /* del center-(high)-freq elements, down to first big high-freq peak?
     * Problem is: if HR = 3 Hz, then 10*HR = 30 Hz, so need 60 Hz frame
     * sampling to reach that. Even for HR = 1 Hz need 20 Hz frame sampling.
     * So it is debatable whether we should do anything to our high-freq
     * re[], which at frame rate 16 Hz would only be ~8 Hz.
    */
    if (filtNoise > 0) {
        aveRe /= (float)filtNoise;
        for (int i = 0; i < lgt; i++) {
            if ((float)Math.abs((double)re[i]) < aveRe) re[i] = im[i] = 0.f;
        }
    }  
    
    /* LOW-FREQ CUTOFF 
       IF FFT.hrPeakLoc = iPeakHRre < 0 that means did not find a HR peak
       THEN FFT.hrFreqHz = FFT.mtD[0] <= 0 
       Use that as flag to NOT write output line for this FFT */
    // FILTER LOW FREQ elements at ends of re[]
    // FFT: 0=DC offset, 1==255, 2==254...127==129, 128 is unique
    lowCutoff = FFT.hrPeakLoc - 1;  // cut off up max peak-1
	// The re[hrPeakLoc - 1] is either similar to re[hrPeakLoc] - means part of HR variability
    //  OR is much smaller - means part of low-freq background
    //  IF the latter, then remove it
    if (Math.abs(re[lowCutoff]) < Math.abs(re[FFT.hrPeakLoc])/2.) lowCutoff++;
    if (filtLow > 0) {
        for (int i = 0; i < lowCutoff; i++) {
            re[i] = im[i] = 0.f;
            if (i > 0) re[lgt - i] = im[lgt - i] = 0.f;
        }
    } // end filter low end
    
    /* IF SHIFT HR is ON, shift half arrays from FFT.hrPeakLoc down to shiftHR
       and upper half up to lgt-shiftHR-1
       THEN zero the middle of the arrays.
       Note: in 256 FFT 0=DC offset, 1-127 & 129-255 = FFT, 128 = unique midpoint
       Also find averages of the 4 values - large HR peaks around shiftHR */
    if ((shiftHR > 0) && (shiftHR < lowCutoff)) {
        for (int i = lowCutoff, sh = shiftHR; i <= lgtHalf; i++, sh++) {
            re[sh] = re[i];
            im[sh] = im[i];
            re[lgt - sh] = re[lgt - i];
            im[lgt - sh] = im[lgt - i];
            if (sh < shiftHR + delHR) {
                aveHRReLow += Math.abs(re[sh]);
                aveHRImLow += Math.abs(im[sh]);
                aveHRReHigh += Math.abs(re[lgt - sh]);
                aveHRImHigh += Math.abs(im[lgt - sh]);
            }
            re[i] = im[i] = 0.f;
            re[lgt - i] = im[lgt - i] = 0.f;
        } // end for i loop to shift values
        // Also reset the peak location in mtDextra[]
        //FFT.mtDextra[0] += shiftHR - FFT.hrPeakLoc;
    } // end if do shift
    // FILTER DEL HR PEAKS
    if (delHR > 0) {
        float dblDelHr = (float)delHR;
        aveHRReLow /= dblDelHr;
        aveHRImLow /= dblDelHr;
        aveHRReHigh /= dblDelHr;
        aveHRImHigh /= dblDelHr;
    }
    // if delHR = delete HR peaks > 0, del those peaks > aveHRRe or aveHRIm
    if ((shiftHR > 0) && (delHR > 0)) {
        for (int i = shiftHR; i < shiftHR + delHR; i++) {
            if ((float)Math.abs((double)re[i]) > aveHRReLow) re[i] = 0.f;
            if ((float)Math.abs((double)im[i]) > aveHRImLow) im[i] = 0.f;
            if ((float)Math.abs((double)re[iEnd - i]) > aveHRReHigh) re[iEnd - i] = 0.f;
            if ((float)Math.abs((double)im[iEnd - i]) > aveHRImHigh) re[iEnd - i] = 0.f;
        }
    } 
 } // End public static void filterFFT(float[] re, float[] im)

    
 /** public static boolean getPhysioCamArrays(long[] tm, int[] yAmpl, int numSampl)
    Reads PhysioCam arrays for sample time and Y amplitude.
  *    Could replace Y with smoothed Y or U or V
  * Amplitudes go into FFT.fftRe[]. The amplitudes are sums of Y-U-V over numPixels (e.g. 240)
  * Sample time differences go into fftIm[]. Time array shows that sample times are not regular,
  *    so the variations in sample times are treated like phases
  * fftRe[] and fftIm[] are normalized to 0-1
  * Y amplitude range (per pixel) is stored in FFT.waveAmpMx. Div by 50 empirical value to make range around 0-2
    Sample input:
        Time: 0,32,102,153,220,291,335,395,453,518,584,641,702,...
        Y: 18042,18134,18273,18288,18356,18394,18486,18581,18666,...
        Ysm: 18042,18134,18273,18288,18356,18394,18486,18581,18666,...
        U: 12305,12272,12275,12233,12235,12231,12245,12112,12184,...
        V: 24503,24508,24479,24467,24425,24419,24378,24404,24379,24282,...
  *	@param tm - time array, long int
  *	@param yAmpl - Y amplitude array 
  * @param numSampl - actual number of samples (frame data) in the arrays, which may not be full
  *        BECAUSE FFT requires array length = 2^x, calculate FFT.nLgt = highest power of 2 <= numSampl
  *        - take only FFT.nLgt values from input
  *        IF nLgt < numSampl, and BECAUSE the start of a recording is sometimes noisier than the rest,
  *        - take the last nLgt points (not the first nLgt)
  */
 public static void getPhysioCamArrays(long[] tm, int[] yAmpl, int numSampl) {
	//long tmMx = 0, tmMn = 100000;
	long tmIntvl = 0, tmIntvlMx = 0, tmIntvlMn = 100000;
	int yAmpMx = 0, yAmpMn = 100000;
	float tmRng, tmIntvlRng, yAmpRng;
	float timeIntvlMin = 0.f, timeIntvlMax = 0.f, timeIntvlRange = 0.f;
    // for SpO2 calcs
    //float Ysum = 0., YsumSq = 0., Usum = 0., UsumSq = 0., Vsum = 0., VsumSq = 0.;
	// find largest 2^exp in input list. Use for FFT
	if (FFT.nLgt > numSampl) {
		//FFT.nLgt = 8192;                             // reset nLgt
		while ((FFT.nLgt = (FFT.nLgt >>> 1)) > numSampl); // divide nLgt/2 until <= numSamples
		if (FFT.nLgt < 1) FFT.nLgt = 1;
	}
	/*QUESTION: IF nLgt < numSampl
	    THEN we can use the START of arrays tm[0,nLgt] yAmpl[0,nLgt]
	                 OR the END of arrays tm[numSampl-nLgt,numSampl] yAmpl[numSampl-nLgt,numSampl]
	    IF the frame-grabbing is slow to start, it would be better to use the 2nd. */
	FFT.startI = numSampl - FFT.nLgt;  // use nLgt datapoints from startI to numSampl-1
	//FFT.inputTime = new float[FFT.nLgt];
	FFT.fftIm = new float[FFT.nLgt]; // time intervals between samples like phases
	FFT.fftRe = new float[FFT.nLgt]; // for amplitudes
	//FFT.inputUampl = new float[FFT.nLgt];
	//FFT.inputVampl = new float[FFT.nLgt];
	
	// FFT is ONLY on range [startI] - [numSamples-1]
	FFT.totalTime = (float)(tm[numSampl - 1] - tm[FFT.startI]); // or (tm[FFT.nLgt - 1] - tm[0]);
	if (FFT.totalTime <= 0.f) FFT.totalTime = 30000.f;
	// find range of Y amplitudes and time intervals (time between one frame and next varies)
	for (int i1 = FFT.startI; i1 < numSampl; i1++) {
		//if (tm[i1] > tmMx) tmMx = tm[i1];
		//if (tm[i1] < tmMn) tmMn = tm[i1];
		if (yAmpl[i1] > yAmpMx) yAmpMx = yAmpl[i1];
		if (yAmpl[i1] < yAmpMn) yAmpMn = yAmpl[i1];
		if (i1 > 0) {
			tmIntvl = tm[i1] - tm[i1-1];
			if (tmIntvl > tmIntvlMx) tmIntvlMx = tmIntvl;
			if (tmIntvl < tmIntvlMn) tmIntvlMn = tmIntvl;
		}
		//if (FFT.inputUampl[i1] > FFT.uAmplMax) FFT.uAmplMax = FFT.inputUampl[i1];
		//if (FFT.inputUampl[i1] < FFT.uAmplMin) FFT.uAmplMin = FFT.inputUampl[i1];
		//Usum += FFT.inputUampl[i1];
		//UsumSq += FFT.inputUampl[i1] * FFT.inputUampl[i1];
		//if (FFT.inputVampl[i1] > FFT.vAmplMax) FFT.vAmplMax = FFT.inputVampl[i1];
		//if (FFT.inputVampl[i1] < FFT.vAmplMin) FFT.vAmplMin = FFT.inputVampl[i1];
		//Vsum += FFT.inputVampl[i1];
		//VsumSq += FFT.inputVampl[i1] * FFT.inputVampl[i1];
	} // end for i1 loop to find max min
	//tmRng = (float)(tmMx - tmMn);
	tmIntvlRng = (float)(tmIntvlMx - tmIntvlMn);
	yAmpRng = (float)(yAmpMx - yAmpMn);
	// Y range is FEATURE - ADD yAmpRng to waveAmpMx
	//if (PhysioCam.numPixelsUsed > 0) FFT.waveAmpMx = yAmpRng / PhysioCam.numPixelsUsed;
	//else FFT.waveAmpMx = yAmpRng / 240.;  // pick default numPixelsUsed=240
	//FFT.waveAmpMx /= 50.;
	
	// copy LAST nLgt ELEMENTS of tm[] and yAmpl[] into FFT.fftRe[] and fftIm[]
	// normalize fftRe[] (amplitude) and fftIm[] (time intervals)
	for (int i2 = startI + 1; i2 < numSampl; i2++) {
		//FFT.inputTime[i2] = (tm[i2] - tmMn) / tmRng;
		FFT.fftIm[i2 - startI] = (tm[i2] - tm[i2-1] - tmIntvlMn) / tmIntvlRng;
		FFT.fftRe[i2 - startI] = (yAmpl[i2] - yAmpMn) / yAmpRng;
	}
	FFT.fftRe[0] = (yAmpl[startI] - yAmpMn) / yAmpRng; //FFT.fftRe[1];
	FFT.fftIm[0] = FFT.fftIm[1];
 } // End public static void getPhysioCamArrays()

} // class FFT