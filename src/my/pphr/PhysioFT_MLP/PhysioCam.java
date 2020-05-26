package my.pphr.PhysioFT_MLP;
/** PhysioCam class
 *  Derived from HRmonitorTSY, Tan ShiYang 2012, and PhysioCam
 *  Intended to collect sensor (including camera YUV) data for N seconds into arrays.
 *  Then use FFT on the arrays to seek periodic functions like Heart Rate HR, SpO2, 
 *		Blood Pressure SBP DBP
 *  Does NOT attempt real-time processing, as that slows the sampling rate.
 * @author HL Seldon 2017
 * Edited: HL Seldon 2017-
 */
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import android.hardware.Camera;
import android.hardware.Camera.PreviewCallback;
import android.os.Bundle;
import android.os.Environment;
import android.app.Activity;
//import android.app.Dialog;
import android.util.Log;
import android.view.SurfaceHolder;
import android.view.SurfaceView;
import android.view.View;
import android.view.View.OnClickListener;
import android.widget.Button;
import android.widget.TextView;
import android.widget.EditText;
import android.widget.ToggleButton;
import java.util.Date;
import java.text.SimpleDateFormat;
import android.text.method.ScrollingMovementMethod;

public class PhysioCam extends Activity {

    public TextView viewAvgHR; //, viewCurrentHR, viewMaxHR, viewMinHR;
    public TextView viewCurrentSpO; //, viewAvgSpO, viewMaxSpO, viewMinSpO;
	public TextView viewAvgResp;
    public TextView viewCurrentY, viewAvgY, viewMaxY, viewMinY;
    //public TextView viewCurrentU, viewAvgU, viewMaxU, viewMinU;
    //public TextView viewCurrentV, viewAvgV, viewMaxV, viewMinV;
    public TextView viewTime, viewNumFrm, viewFrmSec;
	public TextView viewSBPguess, viewDBPguess, viewRBPguess;
	public TextView messageTxt;
	private EditText getBPval, getHRval;
    private boolean experimental = true; // to control display of Y U V params
    public String bpTxt = "111/11", sBP = "0", dBP = "0", bpRng = "0", hrTxt = "11";
	
    private Camera camera; //Camera object in this activity
    private SurfaceView surfacepreview;
    private SurfaceHolder surfaceholder;
    
    /* Arrays to store Y U V
     * index of arrays is numFrames - i.e. one entry for each video Frame.
     * Wraps around at end, so actual index = numFrames % arraySize 
     * Arrays should be saved for offline analysis at end of session. */
    final int arraySize = 1300; //Size of arrays to store Y U V and system time. 2.5min at 16 frame/sec
    public long[] frameTime; //to store time of each frame
    public byte[] frameIsPeak; // flags to identify frames with peak intensity etc
    public int[] frameY; // Y values from frame = intensity
    public int[] frameYsmooth; // fast phones have noisy Y values; this stores smoothed Ys
    public int[] frameU;
    public int[] frameV;
    public int smoothOver = 1; // average Y values between this many frames
    
    int minus1 = -1; // default value for frameY, timePrevPeak
    long minus1L = -1L;
    long timeStart = minus1L; // start measurements. Delay 2 seconds before real reading
    long timeCurrentPeak; // Time to show instant heart rate
    long timeFirstPeak = minus1L, timeLastPeak; // Time to show overall heart rate
    //long time1PeakAgo = minus1L, time2PeakAgo = minus1L, time3PeakAgo = minus1L;
    
    int totalBeat; //Total beat in whole duration, used to show overall heart rate
    float HRCurrent = 0.f, HRmin = 999.f, HRmax = -999.f;
    float HRSum = 0.0f;
    float aveHR = 0.f, aveHRbeats = 0.f; // aveHR is highest peak in real FT. aveHRbeats = total beats / time
    
    float SpOCurrent = 0.f, SpOmin = 999.f, SpOmax = -999.f;
    float SpOSum = 0.f;
    
    // Frame sample params
    static int exclRim = 15; // default 15% of frame rim excluded from sample
    static int heightStart = -1, heightEnd = 0, widthStart = -1, widthEnd = 0; 
    static int rowStep = 1, colStep = 1;  // increments through rows, columns
    public static int numPixelsUsed = 0;
    public static int numFrames; //Total frames processed. 
    public static float framesPerSec = 1.f; // frames processed per second
    
    // Statistics
    static int YFrameSum = 0, YFrameSumMin = 100000, YFrameSumMax = 0;
	static int YFrameAve = 0, YFrameAveMin = 100000, YFrameAveMax = 0;
    static float YAve = 0.f;
    static long YframeAveSum = 0L; //sum over frames of ave Y per frame
    static int UFrameSum = 0, UFrameSumMin = 100000, UFrameSumMax = 0;
	static int UFrameAve = 0, UFrameAveMin = 100000, UFrameAveMax = 0;
    static float UAve = 0.f, UStdDev = 0.f;
    static long UframeAveSum = 0L, UAveSumSqr = 0L;
    static int VFrameSum = 0, VFrameSumMin = 100000, VFrameSumMax = 0;
	static int VFrameAve = 0, VFrameAveMin = 100000, VFrameAveMax = 0;
    static float VAve = 0.f, VStdDev = 0.f;
    static long VframeAveSum = 0L, VAveSumSqr = 0L;
    
    static int byteMask = 0x000000FF; // to convert byte to int
    
	public static File physFile; // pointer to PhysFT_Feat41.csv or PhysFT_TmYUV.csv file
	public static FileWriter physFileWrite;
	public static FileReader physFileRead;
	public BufferedReader inputBuf; // wrapper for physFileRead
	
	// Parameters to LOAD Feat41.csv file for MLP
	public static char driverMode, modeFindMLP = 'f', modeGuessBP = 'g', modeMore = 'm';
	ArrayList<String> fileCnt;
	ArrayList<String> sampleDateTime;
	String inputLine;
	String[] inputHeadrLn = null;
	int numInputLines = 0, numInputVars = 0, numInputVarsUse = 35;
	int iHRbpm = 0;
	public float[][] inputArr; // all input vects. numInputLines x numInputVars
	public int[] inputSBP, inputDBP, inputRBP; // save input SBP DBP RBP vals
	public int vExpect = 0, iExpect = 0;

	// OUTPUT or EXPECTED values = expected SBP, DBP or BPrange
	public int numOutputVals;  // output dimension = num of possible values
	// the measured (expected) values of training set
	public int sbpMin = 80, sbpMax = 200, sbpNumVals; // SBP range of values
	public int dbpMin = 50, dbpMax = 100, dbpNumVals; // DBP
	public int rbpMin = 20, rbpMax = 100, rbpNumVals; // BPrange

	public float[][] outputArr;  // numInputLines x numOutputVals
	// expect... arrays have SBP etc values as 1 in arrays of 0, steps of 2
	// e.g. DBP=62 -> (62-50)/2 = 6, so expectDBParr[i] = {0,0,0,0,0,0,1,0,..0}
	public float[][] expectSBParr; // numInputLines x 60
	public float[][] expectDBParr; // numInputLines x 25
	public float[][] expectRBPArr; // numInputLines x 40
	public static float[] dsrTotalErr = { 0.f, 0.f, 0.f };  // calculated total network error for DBP SBP RBP
	public static float[] overallErrPrev = {999.f, 999.f, 999.f}; // total error of previous iteration
	// MLP Structure params
    // layer 0 = input, 1 = hidden1, 2 = hidden2, 3 = output
    public int nodesPerLayLgt = 4;
    public int[] nodesPerLay = new int[nodesPerLayLgt];
    public float learnRate = 0.15f;   // e.g. 0.15
    public float momentum = 0.01f;
    public float minErr = 0.0001f;    // target ave error/sample e.g. 0.001
    public int maxIter = 20000;     // max iterations
    public float convergeSpeedCutoff = 0.999f;
    // number of MLPs to test for best fit
    public int numMLPtoTrain = 30;

    // Params to make MLPs for BP
	private MLP[] nnMLP = new MLP[3]; // MLP[0]=DBP, [1]=SBP, [2]=RBPrange
    private int dsrTyp = 0;   // 0=DBP, 1=SBP, 2=BPrange
	// Open weights file in Driver, pass to MLP[]
    public String[] wgtsFilename = new String[nnMLP.length]; 

    /** static void decodeYUV420SP(byte[] yuv420sp, int width, int height) 
     * For each video frame (numFrames)
     *    find average Y U V
     *    from middle 50x50 pix of camera preview frame
     * Called from other class PreviewCallBack, so must be static
     * Uses static parameters in PhysioCam.
     * @param yuv420sp preview frame in YUV format
     * @param width width of preview frame
     * @param height height of preview frame
     * return null. Sets values of static variables YUVFrameSum, YUVFrameAve, numFrames
     *      Also sets YUVAveSum, YUVAveSumSqr for SD calcs
     */
    public static void decodeYUV420SP(byte[] yuv420sp, int width, int height) 
    {
        if (yuv420sp == null) return;
        int frameSize = width * height;
        if (frameSize <= 0) return;   
        int y = 0, uvp = 0, u = 0, v = 0;
            
        // find area in image to scan
        if ((heightStart < 0) && (widthEnd < 1)) {
		// default exclRim = 15, so exclude 15% rim around edge of image
		if (exclRim < 0) exclRim = 0;
		if (exclRim > 45) exclRim = 45; // cannot exclude more than 45% on each edge
		heightStart = (height * exclRim) / 100;
		heightEnd = height - heightStart;
		// default sample only 14 rows, and 14 pixels in each row
		rowStep = (heightEnd - heightStart) / 14; // want to sample 14 rows
		widthStart = (width * exclRim) / 100;
		widthEnd = width - widthStart;
		colStep = (widthEnd - widthStart) / 14; // want to sample 14 cols
        } // end first time - find sample window
            
        numPixelsUsed = 0; // counter for actual number of pixels processed
        YFrameSum = UFrameSum = VFrameSum = 0;
        /* SEE ColorSpaceNotes.txt in /docs folder 
         * YUV: First frameSize =height*width bytes are Y. yp = index to these
         * Next frameSize/2 values are UV pairs. uvp = index to these
         * LOOP thru the sample area in middle of frame*/
        // BUG? yp = heightStart * width; at start?? Orig version no width
        for (int row = heightStart, yp = heightStart * width; row < heightEnd; row += rowStep, yp += rowStep * width) {
            uvp = frameSize + (row >> 1) * width;
            for (int i = widthStart; i < widthEnd; i += colStep, yp+=colStep) {
                // Y studio range: 16..236  U&V: 16..240
                y = (PhysioCam.byteMask & ((int) yuv420sp[yp])); // - 16; if convert RGB
                if (y < 0) y = 0;
                YFrameSum += y;
                if (((i - widthStart) & 1) == 0) { // if i is even number from start
                    v = (PhysioCam.byteMask & (int)yuv420sp[uvp++]); // - 128; if convert RGB
                    u = (PhysioCam.byteMask & (int)yuv420sp[uvp++]); // - 128; if convert RGB
                    UFrameSum += u;
                    VFrameSum += v;
                }
                else u = v = 0;
                numPixelsUsed ++; // count number of pix sampled
            } // end for i loop thru pixels of each row 
        } // end loop thru rows. End of one Frame sampling.
          
        // Find max min for frame sums
		if (PhysioCam.YFrameSumMin > PhysioCam.YFrameSum) PhysioCam.YFrameSumMin = PhysioCam.YFrameSum;
		if (PhysioCam.YFrameSumMax < PhysioCam.YFrameSum) PhysioCam.YFrameSumMax = PhysioCam.YFrameSum;
		if (PhysioCam.UFrameSumMin > PhysioCam.UFrameSum) PhysioCam.UFrameSumMin = PhysioCam.UFrameSum;
		if (PhysioCam.UFrameSumMax < PhysioCam.UFrameSum) PhysioCam.UFrameSumMax = PhysioCam.UFrameSum;
		if (PhysioCam.VFrameSumMin > PhysioCam.VFrameSum) PhysioCam.VFrameSumMin = PhysioCam.VFrameSum;
		if (PhysioCam.VFrameSumMax < PhysioCam.VFrameSum) PhysioCam.VFrameSumMax = PhysioCam.VFrameSum;
		// FOR EACH FRAME calculate YUV FrameAve, add to FrameSum, SD
        YFrameAve = (YFrameSum + numPixelsUsed/2) / numPixelsUsed;
        YframeAveSum += YFrameAve;
        //YAveSumSqr += YFrameAve * YFrameAve;
        UFrameAve = (UFrameSum + numPixelsUsed/2) / numPixelsUsed;
        UframeAveSum += UFrameAve;
        UAveSumSqr += UFrameAve * UFrameAve;
        VFrameAve = (VFrameSum + numPixelsUsed/2) / numPixelsUsed;
        VframeAveSum += VFrameAve;
        VAveSumSqr += VFrameAve * VFrameAve;
        // add 1 to number of frames processed
        numFrames++;
    } // END static void decodeYUV420SP(byte[] yuv420sp, 

    /** public void resetVals()
     * resets the parameters, arrays, etc */
    public void resetVals() {
        for(int i = 0; i < arraySize; i++) { 
            frameYsmooth[i] = frameY[i] = frameU[i] = frameV[i] = minus1;
            frameTime[i] = minus1L;
			frameIsPeak[i] = 0;
        }
		timeStart = timeFirstPeak = minus1L;
        HRCurrent = HRSum = 0.0f;
        HRmin = SpOmin = 999.f;
        HRmax = -HRmin;
		totalBeat = 0;
        SpOCurrent = SpOSum = SpOmax = 0.0f;
        PhysioCam.YFrameAve = PhysioCam.UFrameAve = PhysioCam.VFrameAve = 0;
		PhysioCam.YframeAveSum = PhysioCam.UframeAveSum = PhysioCam.VframeAveSum = 0;
        PhysioCam.YFrameSumMin = PhysioCam.UFrameSumMin = PhysioCam.VFrameSumMin = 100000;
        PhysioCam.YFrameSumMax = PhysioCam.UFrameSumMax = PhysioCam.VFrameSumMax = 0;
		PhysioCam.YframeAveSum = PhysioCam.UframeAveSum = PhysioCam.UAveSumSqr = PhysioCam.VframeAveSum = PhysioCam.VAveSumSqr = 0L;
        //YMin = UMin = VMin = 100000;
        //YMax = UMax = VMax = 0;
		PhysioCam.UStdDev = PhysioCam.VStdDev = 0.0f;
        PhysioCam.numFrames = 0;
    } // END public void resetVals()
        
    /** void onCreate(Bundle savedInstanceState)
     * Called when the activity is started by user. 
     * Loads the GUI,
     * show dialog box that tells user how to use the app and
     * assigns necessary textbox and toggle button to variables. 
     * Starts surfaceCB to preview camera' view. 
     * @param savedInstanceState - the current saved state of this instance
     */
    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.physiocam_main);// Load the layout of GUI
        
        // CREATE ARRAYS for SAMPLE TIME, INTENSITY, ...
        frameTime = new long[arraySize];
        frameY = new int[arraySize];
        frameYsmooth = new int[arraySize];
        frameU = new int[arraySize];
        frameV = new int[arraySize];
        frameIsPeak = new byte[arraySize];
		
		// Add Calc button to original PhysioCam file 
		// This passes frameY(smooth)[] to FFT, calculates FFT and metadata FEATURES
        final Button calcBtn = (Button)findViewById(R.id.calc); 
        calcBtn.setOnClickListener(new OnClickListener() {
            public void onClick(View viewParam) {
				//String bpTxt = "111/11", sBP = "0", dBP = "0", bpRng = "0", hrTxt = "11";
				// to count heartbeats, skip first 3 sec to bypass startup noise
				int startPeakCountFrm = (int)(3. * PhysioCam.framesPerSec);
			if (experimental) {
				// get user HR text value, to check if =default 11
				hrTxt = getHRval.getText().toString();
				
				// Choose number samples to process arrays
				PhysioCam.numFrames = PhysioCam.numFrames % arraySize;
				// Choose length of FFT arrays. Must be power of 2
				FFT.nLgt = 8192;                             // reset nLgt
				if (FFT.nLgt > PhysioCam.numFrames) {
					while ((FFT.nLgt = (FFT.nLgt >>> 1)) > PhysioCam.numFrames); // divide nLgt/2 until <= numFrames
					if (FFT.nLgt < 1) FFT.nLgt = 1;
				}
				FFT.startI = PhysioCam.numFrames - FFT.nLgt;  // use nLgt datapoints from startI to numFrames-1
				
				// For fast noisy cameras calculate the average frameYsmooth over "smoothOver" frames
				// weighted ave: a[i] = (3*a[i] + a[i-1] + a[i+1])/5 (or 4 if at start or end of a[])
				smoothOver = (int)(PhysioCam.framesPerSec / 5.);
				if (smoothOver < 1) smoothOver = 1;
				if (smoothOver > 3) smoothOver = 3;
				int smoothI, smoothWgt;
				totalBeat = 0;
				timeFirstPeak = minus1L;
				for (int frI = 0; frI < PhysioCam.numFrames; frI++) {
					frameYsmooth[frI] = 3 * frameY[frI];
					smoothWgt = 3;
					switch(smoothOver) {
						case 3: smoothI = (frI + 1);
								if (smoothI < PhysioCam.numFrames) {
									frameYsmooth[frI] += frameY[smoothI];  smoothWgt++;
								} 
						case 2: smoothI = (frI - 1);
								if (smoothI >= 0) {
									frameYsmooth[frI] += frameY[smoothI]; smoothWgt++;
								} 
						case 1: // frameYsmooth[frI] += 2 * frameY[frI];  smoothWgt += 2;
					}
					frameYsmooth[frI] = (frameYsmooth[frI] + smoothWgt/2) / smoothWgt;
					
					// to count heartbeats, skip first 3 sec to bypass startup noise
					if (frI < startPeakCountFrm) continue;
					/* COUNT PEAKS = totalBeats in frameYsmooth[]
						PEAK is 2 or 3 increasing Y values [0-2] followed by one decreasing [3] */
					if ((frameYsmooth[frI - 4] < frameYsmooth[frI - 3])
							&& (frameYsmooth[frI - 3] < frameYsmooth[frI - 2])
							&& (frameYsmooth[frI - 2] < frameYsmooth[frI - 1])
							&& (frameYsmooth[frI - 1] > frameYsmooth[frI]))
					{
						frameIsPeak[frI - 1] = 1;
						timeCurrentPeak = frameTime[frI - 1];
						//timeFirstPeak and timeCurrentPeak used to show average heart rate
						if (timeFirstPeak == minus1L) timeFirstPeak = timeCurrentPeak; 
						totalBeat++; // Increase total heart beat by 1
					}
					else {
						frameIsPeak[frI - 1] = 0;
					}
				} // end frI loop to make frameYsmooth[] and count peaks
				if (totalBeat > 4) { 
                /* Find ave HR as total beats / total time since start, */
                    aveHRbeats = 60000.f * (float)(totalBeat - 1) / (float)(timeCurrentPeak - timeFirstPeak);
                    viewAvgHR.setText(String.format("%.1f", aveHRbeats));
                } // end if totalBeat>4
				
				/* Y range is FEATURE - FFTp.yAmpl. Normalize per pixel.
				   CALC Y amplitude = range for individual waves in Y[]. 
					The waveAmpMx is the 2nd-highest wave amplitude to avoid outliers.
					Also calc waveAmpAve as average of individual wave amplitudes, ignoring top
					  and bottom values as outliers. */
				int[] waveAmpL = new int[FFT.nLgt];
				int dfY;
				// Calc diff between each Y and previous one, starting at startI
				waveAmpL[0] = 0;
				for (int wv = 1; wv < FFT.nLgt; wv++) {
					dfY = frameYsmooth[wv + FFT.startI] - frameYsmooth[wv + FFT.startI - 1];
					// IF dfY >= 0, ie positive slope, add to previous cumulative positive wave
					if (dfY >= 0) {
						waveAmpL[wv] = waveAmpL[wv-1] + dfY;
						waveAmpL[wv-1] = 0;
					}
					else waveAmpL[wv] = dfY; // do not sum negative slope
				} // end for wv loop to fill waveAmp with cumulative pos diffs in Y
				/* Peak is defined as positive value with >= 2 zeroes (prev pos values)
				   before it and a neg or 0 val after it. 
				   e.g. -10 -23 0 0 0 0 44 -21 0 44 -12 -> first 44 is peak, 2nd is not */
				for (int wv = 2; wv < FFT.nLgt - 1; wv++) {
					if (waveAmpL[wv] > 0) {
						if ((waveAmpL[wv-1] == 0) && (waveAmpL[wv-2] == 0)
							&& (waveAmpL[wv+1] <= 0) ) continue;
						else waveAmpL[wv] = -1; // not a real peak
					} // end if have positive num
				} // end for wv loop to filter too-short peaks
				Arrays.sort(waveAmpL); // sort peak list from low to high
				// choose 2nd-highest peak for Y amplitude. Highest may be outlier
				if (numPixelsUsed > 0) FFT.waveAmpMx = (float)waveAmpL[FFT.nLgt - 2] / numPixelsUsed;
				else FFT.waveAmpMx = (float)waveAmpL[FFT.nLgt - 2] / 240.f;
				FFT.waveAmpMx /= 2.;  // 2 =  empirical value to get range 0-2
				// next get the average ampl
				float waveAmpSum = 0.f;
				int waveAmpCnt = 0;
				for (int wv = 0; wv < FFT.nLgt - 1; wv++) {  // skip top waveAmpL[] as outlier
					if (waveAmpL[wv] <= 0) continue; // only peaks, ie pos numbers
					waveAmpSum += (float)waveAmpL[wv];
					waveAmpCnt ++;
				} // end for wv loop to calc average
				if (waveAmpCnt > 0) {
					if (numPixelsUsed > 0) FFT.waveAmpAve = waveAmpSum / (2.f * waveAmpCnt * numPixelsUsed);
					else FFT.waveAmpAve = waveAmpSum / (2.f * waveAmpCnt * 240.f);
				}
				else FFT.waveAmpAve = 0.f;
				
				// QUESTION: use frameY[] or frameYsmooth[] arrays? Peak detection better in smoothed
				FFT.getPhysioCamArrays(frameTime, frameYsmooth, PhysioCam.numFrames);
				FFT.doFFT(false, false); // do FFT, find physio features in FFT, NO filter, NO inverse FFT
				
				// Calc aveHR from FFT, display it in the editable HRbpm box
				aveHR = (60.f * FFT.HR_Hz);  // HR calculated from FFT peaks
				//viewAvgHR.setText(String.format("%d", FFT.mtD[0])); //= ((int)((60. * FFT.HR_Hz)))
				getHRval.setText(String.format("%d", (int)(aveHR + 0.5)));
				// Calc ave respiratory rate. For now, do NOT display SpO2 current
				viewAvgResp.setText(String.format("%.1f", FFT.Resp_Hz * 60.));
				//viewCurrentSpO.setText(String.format("%.1f", SpOCurrent));
				messageTxt.setText("HR FT Hz=" + String.format("%.2f",FFT.HR_Hz) + " bpm="
						+ String.format("%.1f",aveHR) + "\nHR count Beats=" + String.format("%d",totalBeat) 
						+ " in sec=" + String.format("%.0f", (timeCurrentPeak - timeFirstPeak) / 1000.));
				//resetVals();
            } // end if experimental
            } // End calcBtn onClick()
        }); // end calc button Listener Object - 
    	
		/* Add Save button to save T Y U V, FFT FEATURES to file
		   saveTMYUVtoFile() is SET TO NOT SAVE T Y U V (false), and TO SAVE FFT Features (true) */
		final Button saveHR = (Button)findViewById(R.id.save); 
        saveHR.setOnClickListener(new OnClickListener() {
            public void onClick(View viewParam) {
				/* Prompt for blood pressure */
				bpTxt = getBPval.getText().toString();
				int bpTxtDiv = bpTxt.indexOf("/");
				if (bpTxtDiv > 0) {
					sBP = bpTxt.substring(0, bpTxtDiv);
					dBP = bpTxt.substring(bpTxtDiv + 1);
					try { bpRng = String.valueOf(Integer.parseInt(sBP) - Integer.parseInt(dBP));
					} catch(NumberFormatException nfe) {
						bpRng = "NaN";	}
				}
				hrTxt = getHRval.getText().toString();
				String saveHR_BP = hrTxt /*String.valueOf((int)aveHR)*/ + "," + sBP + "," + dBP + "," + bpRng + ",0," + String.valueOf(PhysioCam.numPixelsUsed);
				saveTmYUV_Feat_File(saveHR_BP, false, true); //} // false = do NOT save YUV, true = DO save Features
				resetVals();
			 } // End onClick()
        }); // end saveHR button Listener Object - save instantHR data into pphr_android
		
		// Add Train button to train MultiLayer Perceptron MLP
        final Button trainBtn = (Button)findViewById(R.id.train); 
        trainBtn.setOnClickListener(new OnClickListener() {
            public void onClick(View viewParam) {
				PhysioCam.driverMode = PhysioCam.modeFindMLP;
				loadFeat41csv();
				makeMLPs();
				runMLPs();
			} // End onClick()
        }); // end Train button Listener
	
		// Add Get BP btn
		final Button guessBPBtn = (Button)findViewById(R.id.getBP);
		guessBPBtn.setOnClickListener(new OnClickListener() {
			public void onClick(View viewParam) {
				PhysioCam.driverMode = PhysioCam.modeGuessBP;
				loadFeat41csv();
				makeMLPs();
				runMLPs();
			} // end onClick()
		} ); // end guessBPBtn listener
		
		// Add Quit button 
        final Button quitBtn = (Button)findViewById(R.id.quit); 
        quitBtn.setOnClickListener(new OnClickListener() {
            public void onClick(View viewParam) {
                finish();
            } // End onClick()
        }); // end quit button Listener
		
		// VIEW RESULTS
        viewTime = (TextView)findViewById(R.id.tvTimeValue);
		viewNumFrm = (TextView)findViewById(R.id.tvNpxValue);
        viewFrmSec = (TextView)findViewById(R.id.tvFrmSecValue);
		viewAvgHR = (TextView)findViewById(R.id.tvHRAvgValue);
		viewAvgResp = (TextView)findViewById(R.id.tvRespAvgValue);
		viewCurrentSpO = (TextView)findViewById(R.id.tvSpOValue);
		//viewCurrentHR = (TextView)findViewById(R.id.tvHRValue);
        //viewMaxHR = (TextView)findViewById(R.id.tvHRMaxValue);
        //viewMinHR = (TextView)findViewById(R.id.tvHRMinValue);
        //viewAvgSpO = (TextView)findViewById(R.id.tvSpOAvgValue);
        //viewMaxSpO = (TextView)findViewById(R.id.tvSpOMaxValue);
        //viewMinSpO = (TextView)findViewById(R.id.tvSpOMinValue);
        
        if (experimental) {
            viewCurrentY = (TextView)findViewById(R.id.tvYValue);
            viewAvgY = (TextView)findViewById(R.id.tvYAvgValue);
            viewMaxY = (TextView)findViewById(R.id.tvYFrameAveMaxValue);
            viewMinY = (TextView)findViewById(R.id.tvYFrameAveMinValue);
			viewSBPguess = (TextView)findViewById(R.id.tvSBPvalue);
			viewDBPguess = (TextView)findViewById(R.id.tvDBPvalue);
			viewRBPguess = (TextView)findViewById(R.id.tvRBPvalue);
			getHRval = (EditText)findViewById(R.id.hrValue);
			getBPval = (EditText)findViewById(R.id.bpValue);
			messageTxt = (TextView)findViewById(R.id.msgTxt);
			messageTxt.setMovementMethod(new ScrollingMovementMethod());
            //viewCurrentU = (TextView)findViewById(R.id.tvUValue);
            //viewAvgU = (TextView)findViewById(R.id.tvUAvgValue);
            //viewMaxU = (TextView)findViewById(R.id.tvUFrameAveMaxValue);
            //viewMinU = (TextView)findViewById(R.id.tvUFrameAveMinValue);
            //viewCurrentV = (TextView)findViewById(R.id.tvVValue);
            //viewAvgV = (TextView)findViewById(R.id.tvVAvgValue);
            //viewMaxV = (TextView)findViewById(R.id.tvVFrameAveMaxValue);
            //viewMinV = (TextView)findViewById(R.id.tvVFrameAveMinValue);
        } // end if experimental, show Y U V
                
        surfacepreview = (SurfaceView)findViewById(R.id.cameraPreview);
        surfaceholder = surfacepreview.getHolder();
        surfaceholder.addCallback(surfaceCB);
        surfaceholder.setType(SurfaceHolder.SURFACE_TYPE_PUSH_BUFFERS);
    }  // END onCreate()
    

    /** SurfaceHolder.Callback surfaceCB
     * Called to hold preview of camera. 
     * It opens the camera and set the effect to none, 
     * set the preview size to the smallest
     *  and start the preview of camera. 
     */
    public SurfaceHolder.Callback surfaceCB = new SurfaceHolder.Callback ()
    {
        public void surfaceCreated(SurfaceHolder holder)
        {
            camera = Camera.open(); //To access camera   
            Camera.Parameters param = camera.getParameters();   
            param.setColorEffect(Camera.Parameters.EFFECT_NONE); //Set effect of preview to NONE
            
            List<Camera.Size> camerasize = param.getSupportedPreviewSizes(); //Get supported preview size
            //Set lowest preview size
            param.setPreviewSize(camerasize.get(camerasize.size()-1).width,
                    camerasize.get(camerasize.size()-1).height);
            camera.setDisplayOrientation(90); // set vertical orientation
            camera.setParameters(param); //Make necessary changes to camera parameters
            camera.startPreview(); //Start to preview           
            try {
                camera.setPreviewDisplay(surfaceholder); //Set camera preview in specify surface holder
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        } // End surfaceCreated(...)
        
        public void surfaceChanged(SurfaceHolder holder, int format, int width, int height) {
        }
        
        public void surfaceDestroyed(SurfaceHolder holder) {           
        }
    }; // END SurfaceHolder.Callback object

    /** void onStartStopClicked(View v)
     * Listener for toggleStartStop ToggleButton in layout xml:
     *     android:id="@+id/toggleStartStop"
     * To Start: Reset variable values,  turn on camera video preview. 
     * To Stop: delete camera preview. 
     * @param v - the parent View (ToggleButton in layout main.xml)
     */    
    public void onStartStopClicked(View v)  {
        if (((ToggleButton) v).isChecked()) { // START to detect heart rate
            resetVals();    
            camera.setPreviewCallback(previewCB);
        } // end toggle to START
        else  { // STOP detecting heart rate
            camera.setPreviewCallback(null); //Stop processing preview frame
        }
    } // END onStartStopClicked()

    /** PreviewCallback preview
     * Gets preview frame from camera and process it. 
     * Extracts Y-U-V values from sample of numPix pixels. 
	 * Stores them in frameTime[] frameY[] frameU[] frameV[]
	 * Calculates averages, standard deviations
	 * Uses empirical linear regression to calculate SpO2 from from the Y-U-V stats.
     */
    private PreviewCallback previewCB = new PreviewCallback ()
    {
        public void onPreviewFrame(byte[] data, Camera camera)
        { 
            long currentInterval = 1;
            
            /* decodeYUV... see description above
             * numFrames++ = number of video frames processed */
            PhysioCam.decodeYUV420SP(data, camera.getParameters().getPreviewSize().width, camera.getParameters().getPreviewSize().height);
            int frameI = (numFrames - 1) % arraySize; // index to Frame arrays
            frameTime[frameI] = System.currentTimeMillis();
            if (timeStart == minus1L) timeStart = frameTime[frameI]; // initialize start time
            // save Y U V
			frameY[frameI] = PhysioCam.YFrameSum;
			frameU[frameI] = PhysioCam.UFrameSum;
            frameV[frameI] = PhysioCam.VFrameSum;
			
			// For fast cameras calculate the weighted average YFrameSum over "smoothOver" frames
            // Do that under the 'Save' button because do not need it in real time
			
            viewTime.setText(String.format("%d", (frameTime[frameI] - frameTime[0])/1000));
			viewNumFrm.setText(String.format("%d", PhysioCam.numFrames));
            // after a startup period, calculate framesPerSecond and smoothOver (number of frames to average)
            if (PhysioCam.numFrames > 60) {
                PhysioCam.framesPerSec = (1000.f * PhysioCam.numFrames) / (float)(frameTime[frameI] - timeStart);
				viewFrmSec.setText(String.format("%.1f", PhysioCam.framesPerSec));
            }
			
            // YUV min max ave are calculated and displayed
            PhysioCam.YFrameAveMin = (PhysioCam.YFrameSumMin + numPixelsUsed/2) / numPixelsUsed;
			PhysioCam.YFrameAveMax = (PhysioCam.YFrameSumMax + numPixelsUsed/2) / numPixelsUsed;
			PhysioCam.UFrameAveMin = (PhysioCam.UFrameSumMin + numPixelsUsed/2) / numPixelsUsed;
			PhysioCam.UFrameAveMax = (PhysioCam.UFrameSumMax + numPixelsUsed/2) / numPixelsUsed;
			PhysioCam.VFrameAveMin = (PhysioCam.VFrameSumMin + numPixelsUsed/2) / numPixelsUsed;
			PhysioCam.VFrameAveMax = (PhysioCam.VFrameSumMax + numPixelsUsed/2) / numPixelsUsed;
            if (experimental) {
                viewCurrentY.setText(String.format("%d", PhysioCam.YFrameAve));
                viewMinY.setText(String.format("%d", PhysioCam.YFrameAveMin));
                viewMaxY.setText(String.format("%d", PhysioCam.YFrameAveMax));
                //viewCurrentU.setText(String.format("%d", PhysioCam.UFrameAve));
                //viewMinU.setText(String.format("%d", PhysioCam.UFrameAveMin));
                //viewMaxU.setText(String.format("%d", PhysioCam.UFrameAveMax));
                //viewCurrentV.setText(String.format("%d", PhysioCam.VFrameAve));
                //viewMinV.setText(String.format("%d", PhysioCam.VFrameAveMin));
                //viewMaxV.setText(String.format("%d", PhysioCam.VFrameAveMax));
            } // end if experimental, show Y U V
            // Calc YUV overall averages. Display
            if (PhysioCam.numFrames > 0) {
                PhysioCam.YAve = (float)PhysioCam.YframeAveSum / PhysioCam.numFrames;
                PhysioCam.UAve = (float)PhysioCam.UframeAveSum / PhysioCam.numFrames;
                PhysioCam.VAve = (float)PhysioCam.VframeAveSum / PhysioCam.numFrames;
                if (experimental) {
                    viewAvgY.setText(String.format("%.1f", PhysioCam.YAve));
                    //viewAvgU.setText(String.format("%.1f", PhysioCam.UAve));
                    //viewAvgV.setText(String.format("%.1f", PhysioCam.VAve));
                } // end if experimental, show Y U V
            }  // end if numFrames > 0
            // if numFrames >1, calc overall YUV StdDev. Display.
            if (PhysioCam.numFrames > 1) {
                //PhysioCam.YStdDev = Math.sqrt((PhysioCam.YAveSumSqr - PhysioCam.YframeAveSum*PhysioCam.YframeAveSum/(double)PhysioCam.numFrames) / (PhysioCam.numFrames-1));
                PhysioCam.UStdDev = (float)Math.sqrt((PhysioCam.UAveSumSqr - PhysioCam.UframeAveSum*PhysioCam.UframeAveSum/(double)PhysioCam.numFrames) / (PhysioCam.numFrames-1));
                PhysioCam.VStdDev = (float)Math.sqrt((PhysioCam.VAveSumSqr - PhysioCam.VframeAveSum*PhysioCam.VframeAveSum/(double)PhysioCam.numFrames) / (PhysioCam.numFrames-1));
                // Convert hue to SpO2 with (R = 0.67, Rsquare = 0.45)
                SpOCurrent = 108.31f     // 108.362 - 0.0013 * 40. altitude (fixed at 40 m)
                    + (1.04f * PhysioCam.UAve - 1.62f * PhysioCam.VAve) / PhysioCam.YAve
                    + (14.62f * PhysioCam.UFrameAveMax - 34.15f * PhysioCam.UStdDev) / PhysioCam.UAve
                    + (40.8f * PhysioCam.VStdDev - 22.01f * PhysioCam.VFrameAveMax) / PhysioCam.VAve;
                FFT.spO2 = SpOCurrent / 100.f; // save in range 0-1 for normalization
				viewCurrentSpO.setText(String.format("%.0f", SpOCurrent));
			}  // end if numFrames > 1
        } //end of onPreviewFrame method
    };  // END of new PreviewCallback object
    
    /** public void saveTmYUV_Feat_File(String hrBPstr, boolean saveYUVflag, boolean saveFeatFlag)
	 *  FIRST writes time, Y U V array to Download/Phys_YUV.csv
	 *  Saves ONLY the FFT.nLgt elements at ends of arrays - discard beginnings.
	 *  SECOND writes features from FFT analysis, and user HR BP input to Download/Phys_Feat41.csv
	 *  @param hrBPstr - String with user-given HR BP values, numPixels used in video frame
	 *  @param saveYUVflag - boolean to save the Time Y U V values to storage
	 *  @param saveFeatFlag - boolean to save the Features derived from FFT to storage */
	public void saveTmYUV_Feat_File(String hrBPstr, boolean saveYUVflag, boolean saveFeatFlag) {
		String commaS = ",";
		Date dNow;
		SimpleDateFormat dNowFt;
		dNowFt = new SimpleDateFormat("yyyyMMdd_HHmmss");
		dNow = new Date( );
        String currentTimeStr = dNowFt.format(dNow) + commaS;
		//String currentTimeStr = String.valueOf(System.currentTimeMillis()) + commaS;
		int saveUpTo = PhysioCam.numFrames % arraySize; // if wrap-around in arrays, only save newest frames
		String saveTm = "Time" + currentTimeStr;
		String saveY = "Y" + currentTimeStr;
		String saveYsmooth = "Ysm" + currentTimeStr;
		String saveU = "U" + currentTimeStr;
		String saveV = "V" + currentTimeStr;
		//String saveFFTre = "YFFTre" + currentTimeStr;
		//String saveFFTim = "YFFTim" + currentTimeStr;
		String saveFFTfeat = "FTfeatVal" + currentTimeStr;
		// CREATE Time Y U V Strings to SAVE
		if (saveYUVflag) {
			// SAVE ONLY the part used for FFT, ie from FFT.startI. Discard early values.
			for (int i = FFT.startI; i < saveUpTo; i++) { 
				if (frameTime[i] <= 0L) break;
				saveTm += String.valueOf(frameTime[i] - frameTime[FFT.startI]) + commaS;
				saveY += String.valueOf(frameY[i]) + commaS; // - PhysioCam.YFrameSumMin
				saveYsmooth += String.valueOf(frameYsmooth[i]) + commaS; // - PhysioCam.YFrameSumMin
				saveU += String.valueOf(frameU[i]) + commaS; // - PhysioCam.UFrameSumMin
				saveV += String.valueOf(frameV[i]) + commaS; // - PhysioCam.VFrameSumMin
			}
		} // end if saveYUVflag
		
		// CREATE FEATURE STRINGS to SAVE
		if (saveFeatFlag) {
			saveFFTfeat += String.format("%.7f,", FFT.HR_Hz); 
			saveFFTfeat += String.format("%.7f,", FFT.HR_Sec); 
			saveFFTfeat += String.format("%.7f,", FFT.Hm1_Hz); 
			saveFFTfeat += String.format("%.7f,", FFT.Hm1Rel); 
			saveFFTfeat += String.format("%.7f,", FFT.Hm1_Sec); 
			saveFFTfeat += String.format("%.7f,", FFT.Hm2_Hz); 
			saveFFTfeat += String.format("%.7f,", FFT.Hm2Rel); 
			saveFFTfeat += String.format("%.7f,", FFT.Hm2_Sec); 
			saveFFTfeat += String.format("%.7f,", FFT.Resp_Hz); 
			saveFFTfeat += String.format("%.7f,", FFT.RespRel); 
			saveFFTfeat += String.format("%.7f,", FFT.Resp_Sec); 
			saveFFTfeat += String.format("%.7f,", FFT.pow01HR); 
			saveFFTfeat += String.format("%.7f,", FFT.pow01HRa); 
			saveFFTfeat += String.format("%.7f,", FFT.pow02Hm1); 
			saveFFTfeat += String.format("%.7f,", FFT.pow03Hm2); 
			saveFFTfeat += String.format("%.7f,", FFT.pow04Resp); 
			saveFFTfeat += String.format("%.7f,", FFT.powRp01Hm1_HR); 
			saveFFTfeat += String.format("%.7f,", FFT.powRp02Hm2_HR); 
			saveFFTfeat += String.format("%.7f,", FFT.powRp03Resp_HR); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa03Hm12_HR1); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa05Hm12_HRe); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa06Hm12_Hm1e); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa07Hm1e_HR1); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa09Hm2e_HR1); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa10Hm2e_HR2); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa11Hm2e_HRe); 
			saveFFTfeat += String.format("%.7f,", FFT.HRV); 
			saveFFTfeat += String.format("%.7f,", FFT.waveAmpMx);
			saveFFTfeat += String.format("%.7f,", FFT.waveAmpAve);
			saveFFTfeat += String.format("%.7f,", FFT.pow05HR1); // lin regr says following ones dubious
			saveFFTfeat += String.format("%.7f,", FFT.pow06HR2); 
			saveFFTfeat += String.format("%.7f,", FFT.pow07HRe); 
			saveFFTfeat += String.format("%.7f,", FFT.pow08Hm12); 
			saveFFTfeat += String.format("%.7f,", FFT.pow09Hm1e); 
			saveFFTfeat += String.format("%.7f,", FFT.pow10Hm2e); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa01HR1_HR2); // lin regr says following NA
			saveFFTfeat += String.format("%.7f,", FFT.powRa02HR1_HRe); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa04Hm12_HR2); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa08Hm1e_HRe); 
			saveFFTfeat += String.format("%.7f,", FFT.powRa12Hm2e_Hm1e); 
			saveFFTfeat += String.format("%.7f,", FFT.spO2); 

			// add HR, SBP, DBP, BPrng, RecNum to saveFFTfeat
			saveFFTfeat += hrBPstr;  //String.valueOf((int)aveHR) + commaS;
			//saveFFTfeat += sBP + commaS + dBP + commaS + bpRng + commaS + "0," + String.valueOf(PhysioCam.numPixelsUsed);
		} // end if saveFeatFlag
		
		try {
			// WRITE Time YUV to FILE
			if (saveYUVflag) {
				PhysioCam.physFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS), 
				"PhysFT_YUV.csv"); // + String.valueOf(System.currentTimeMillis()) + ".csv");
				PhysioCam.physFileWrite = new FileWriter(PhysioCam.physFile, true);
				PhysioCam.physFileWrite.write(saveTm +"\n" + saveYsmooth +"\n" + saveU +"\n" + saveV +"\n"); // + saveY +"\n"
				PhysioCam.physFileWrite.close();
				saveTm = saveY = saveYsmooth = saveU = saveV = null;
			}
			
			// WRITE FEATURES TO FILE
			if (saveFeatFlag) {
				PhysioCam.physFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS), 
				"PhysFT_Feat41.csv"); // + String.valueOf(System.currentTimeMillis()) + ".csv");
				PhysioCam.physFileWrite = new FileWriter(PhysioCam.physFile, true);
				if (PhysioCam.physFile.length() < 50L)  // could also try ! .exists()
					PhysioCam.physFileWrite.write("FTfeatLbl," + FFT.outHeadr + "\n"); //write header first time 
				PhysioCam.physFileWrite.write(saveFFTfeat +"\n"); // write data every time
				PhysioCam.physFileWrite.close(); // end writing to .csv file
				messageTxt.setText("Save PhysFT_Feat41.csv\n" + saveFFTfeat.substring(0, 25) + "\n" + hrBPstr);
			}
		} // end try to output file
		catch (NullPointerException pe) {
			Log.d("OutFile", "Error creating file: " + pe.getMessage());
		} catch (FileNotFoundException ne) {
			Log.d("OutFile", "File not found: " + ne.getMessage());
		} catch (IOException oe) {
			Log.d("OutFile", "Error accessing file: " + oe.getMessage());
		} 
	} // end saveTmYUV_Feat_File()
	
	/** public void loadFeat41csv()
     *  Input sample lines look like:
     *      FTfeatLbl,HR_Hz,HR_Sec,Hm1_Hz,Hm1Rel,... or
     *      FTfeatVal20190323_114953,1.83848108,0.543927273,3.342692873,...
     *  Loads all input line data into ArrayList fileCnt,
	 *      and the date-time, e.g.20190323_114953, into ArrayList sampleDateTime.
     *  Loops through fileCnt, 
	 *      extracts header lines, then removes them from fileCnt
	 *      looks for default SBP-DBP values ",111,11," (means BP UNKNOWN)
	 *        if 'findMLP' mode, removes those lines
	 *        if 'guessBP' mode, keeps those lines and removes all others
     *  Remaining fileCnt is number of input sample (data) lines
     *  Create inputArr[][] for number of input samples and total number of features (41)
     *      Vars are sorted so that only the first 'numInputVarsUse' e.g. 35 ones are used for calcs
     *  Create expectNNNarr[][]... arrays for number of expected output values
     *     SBP 60 (=200-80 in steps of 2), DBP 25 (=100-50), RBP 40 (=100-20)
     *  Loops through ArrayList fileCnt
     *     Reads input line i
     *     splits into parts
     *        gets known (training) HR SBP DBP RBP from end of line, i.e. end of split list
     *        puts first (35) feature vals into inputArr[i][],
     *           and expected results into expectNNNarr[i][] as 0 0 .. 1 0 .. 0,
     *           and DBP after (35) feature vals
     *  When run MLP, outputArr[][] = corresponding expectNNNarr[][]
     */
	public void loadFeat41csv() {
		this.fileCnt = new ArrayList<String>(); 
		this.sampleDateTime = new ArrayList<String>();
		this.sbpNumVals = (this.sbpMax - this.sbpMin)/2;
		this.dbpNumVals = (this.dbpMax - this.dbpMin)/2;
		this.rbpNumVals = (this.rbpMax - this.rbpMin)/2;
		String[] lineParts;
		// OPEN Downloads/PhysFT_Feat41.csv, load lines into ArrayList
		int i1 = 0;
		try {
			PhysioCam.physFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS), 
			"PhysFT_Feat41.csv");
			PhysioCam.physFileRead = new FileReader(PhysioCam.physFile);
			inputBuf = new BufferedReader(PhysioCam.physFileRead);
			while ((inputLine = inputBuf.readLine()) != null) {
				// Check data lines 'FTfeatVals20190506_123456'
				if (this.inputLine.startsWith("FTfeatVal")) {
                    // lines with unknown BP (SBP 111 DBP 11 RBP 100) NOT used for training (modeFind, modeTrainmore)
                    if ((PhysioCam.driverMode != PhysioCam.modeGuessBP) 
                        && (this.inputLine.indexOf(",111,11,100") >= 0) ) { 
                        }
                    // else save line in ArrayList fileCnt
                    else {
                        i1 = this.inputLine.indexOf(",");
                        if (i1 > 9) sampleDateTime.add(this.inputLine.substring(9, i1));
                        //else sampleDateTime.add(" ");
                        fileCnt.add(this.inputLine.substring(i1 + 1)); // data after FTfeatVal20190506_123456,
                    }
                } // end FtfeatVal lines
				// Check header lines 'FTfeatLbl'
                else if (this.inputLine.startsWith("FTfeatLbl")) {
                    this.inputLine = this.inputLine.substring(this.inputLine.indexOf(",") + 1); // trim FTfeatLbl,
                    if (this.inputHeadrLn == null) {       // process first header line
                        this.inputHeadrLn = this.inputLine.split(",");
                        for (int r2 = 0; r2 < this.inputHeadrLn.length; r2++) {
                            if (this.inputHeadrLn[r2].equals("HRbpm")) this.iHRbpm = r2;
                        }
                    } // end process first header line, ignore other header lines
                } // end if inputLine starts with FTfeatLbl
			} // end while reading inputLine
			inputBuf.close();
		} // end try
		catch (FileNotFoundException fne) {
			messageTxt.setText("PhysFT_Feat41 ERR not found");
		}
		catch (IOException ioe) {
			messageTxt.setText("PhysFT_Feat41 ERR open or read input");
		}
		
		this.numInputLines = fileCnt.size(); // num samples = arraylist size after filtering
		// IF <100 data lines, exit with "<100 measurements"
		if ((PhysioCam.driverMode != PhysioCam.modeGuessBP) && (this.numInputLines < 100)) {
			messageTxt.setText("numInputLines <100 measurements. Not enough");
			return;
		}
		/* CONVERT feature data to inputArr[num cases][41 features]
		   AND SBP DBP RBP data to expectedOutput[num cases][num possible values] */
		// user input. Use FIRST numVarsUse. Add 1 place for DBP value for SBP calcs
		this.numInputVars = this.iHRbpm;  //this.numVarsUse; // + 1; 
		// Make array to hold ALL the input feature variables
		this.inputArr = new float[this.numInputLines][this.numInputVars];
		this.inputSBP = new int[this.numInputLines];
		this.inputDBP = new int[this.numInputLines];
		this.inputRBP = new int[this.numInputLines];
		// Make array to hold the expected number of output values for SBP DBP RBP
		expectSBParr = new float[numInputLines][this.sbpNumVals];  // SBP expect 80-200 / 2
		expectDBParr = new float[numInputLines][this.dbpNumVals];  // DBP expect 50-100 / 2
		expectRBPArr = new float[numInputLines][this.rbpNumVals]; // range expect 20-100 / 2
		// outputArr will correspond to expect[SBP DBP RBP]arr 
		// read input feature lines, save feature values in inputArr
		// save expected output values in expect[SBP DBP RBP]arr
		for (int iL = 0; iL < this.numInputLines; iL++) {
			// zero expected output arrays for input line i
			for (int i2 = 0; i2 < this.sbpNumVals; i2++) {
				if (i2 < this.dbpNumVals) expectDBParr[iL][i2] = expectSBParr[iL][i2] = expectRBPArr[iL][i2] = 0.f;
				else if (i2 < this.rbpNumVals) expectSBParr[iL][i2] = expectRBPArr[iL][i2] = 0.f;
				else expectSBParr[iL][i2] = 0.f;
			}
			lineParts = (fileCnt.get(iL)).split(",");
			/* use first line to print param values
			if (i == 0) {
				System.out.printf("    Total line lgt=%d, numVarsUse=%d, labelHR_SBP_DBP_RBP..=%d\n",
					lineParts.length, this.numVarsUse, lineParts.length - this.iHRbpm);
			} */
			// convert lineParts[] into feature values
			try {
				for (int iV = 0; iV < this.numInputVars; iV++) { // ONLY 35 relevant. Skip [0]
					inputArr[iL][iV] = Float.parseFloat(lineParts[iV]);
				} // end for ff loop getting feature vals from one line
				// next get SBP DBP RBP vals from end of line. Fill expectArrays
				this.inputSBP[iL] = vExpect = Integer.parseInt(lineParts[iHRbpm + 1]); // SBP val
				iExpect = (vExpect - this.sbpMin) / 2;  // element in expectSBParr
				if (iExpect < 0) iExpect = 0;
				if (iExpect >= expectSBParr[0].length) iExpect = expectSBParr[0].length - 1;
				expectSBParr[iL][iExpect] = 1.f; // set element = 1
				//System.out.printf("Ln %3d SBP %3d -> %2d .", i, vExpect, iExpect);

				this.inputDBP[iL] = vExpect = Integer.parseInt(lineParts[iHRbpm + 2]); // DBP
				iExpect = (vExpect - this.dbpMin) / 2;
				if (iExpect < 0) iExpect = 0;
				if (iExpect >= expectDBParr[0].length) iExpect = expectDBParr[0].length - 1;
				expectDBParr[iL][iExpect] = 1.f;
				//System.out.printf(" DBP %3d -> %2d .", vExpect, iExpect);

				this.inputRBP[iL] = vExpect = Integer.parseInt(lineParts[iHRbpm + 3]); // RBP
				iExpect = (vExpect - this.rbpMin) / 2;
				if (iExpect < 0) iExpect = 0;
				if (iExpect >= expectRBPArr[0].length) iExpect = expectRBPArr[0].length - 1;
				expectRBPArr[iL][iExpect] = 1.f;
				//System.out.printf(" RBP %3d -> %2d\n", vExpect, iExpect);
			} // end try
			catch(NumberFormatException nne) {
				messageTxt.setText("ERR: input vars not num -" + nne.getMessage());
			}
		} // end for i loop reading input 
		/* NORMALIZE input columns (features)
           - find mean m, SD sd for each column
           - replace each xk = (xk - m) / sd */
        float fMn, fSD, fSum, fSum2, fVar;
        for (int nV = 0; nV < this.numInputVars; nV++) {  // LOOP thru features
            fSum = fSum2 = 0.f;
            for (int nL = 0; nL < this.numInputLines; nL++) { // for each feature, loop entries
                fSum += inputArr[nL][nV];
                fSum2 += inputArr[nL][nV] * inputArr[nL][nV];
            } // end nL loop adding entries
            fMn = fSum / this.numInputLines;  // MEAN
            fSum2 = fSum2 - (fSum * fSum) / this.numInputLines; // Sum of squares
            fVar = fSum2 / (this.numInputLines - 1);  // Variance
            if (fVar <= 0.0f) fVar = 1.f;
            fSD = (float)Math.sqrt(fVar);
            for (int nL = 0; nL < this.numInputLines; nL++) { // replace vals with norms
                inputArr[nL][nV] = (inputArr[nL][nV] - fMn) / fSD;
            } // end nL loop to normalize feature values
        } // end nV loop thru features
		messageTxt.setText(String.valueOf(PhysioCam.driverMode) + " MLP: Num Samples=" + this.numInputLines 
			+ "\n      Num Feat=" + this.numInputVars + " -use=" + this.numInputVarsUse);
	} // End loadFeat41csv()
	
	/** public void makeMLPs() */
	public void makeMLPs() {
		// BIG LOOP to train or test DBP SBP RBP - dsrTyp
		for (dsrTyp = 0; dsrTyp < nnMLP.length; dsrTyp++) {
			/* select corresponding array of expected outputs
			   0-DBP 1-SBP 2-RBP */
			if (dsrTyp == 0) {  outputArr = expectDBParr;  }
			else if (dsrTyp == 1) {  outputArr = expectSBParr; }
			else if (dsrTyp == 2) {  outputArr = expectRBPArr; }
			numOutputVals = outputArr[0].length;
			//System.out.printf("MODE %d TYPE %s numInputVars %d numOutputVals %d wgtsFilename %s\n",
			//	this.driverMode, MLP.makeBPtyp[dsrTyp], this.numInputVars, 
			//	this.numOutputVals, this.wgtsFilename[dsrTyp]);

			/* IF 'findmlp' mode 'f', i.e. first search for best MLP,
			   THEN create layers and nodes from inputArr dimensions. */
			if (PhysioCam.driverMode == PhysioCam.modeFindMLP) {
				// Create list of 3-4 layers with number nodes in each layer 
				nodesPerLay = new int[nodesPerLayLgt]; 
				// Layer 0 = input. num input dim = num nodes in layer 0
				nodesPerLay[0] = numInputVarsUse; // subset of numInputVars
				nodesPerLay[1] = 3 * nodesPerLay[0]; // or 2
				nodesPerLay[2] = numOutputVals + 3; 
				nodesPerLay[nodesPerLay.length - 1] = numOutputVals;
				convergeSpeedCutoff = 0.999f;
				wgtsFilename[dsrTyp] = null;
			} // end if (PhysioCam.driverMode == 0)
			/* ELSE IF 'guessBP' mode, open saved MLP weights file OR set wgtsFilename to default*/
			else if (PhysioCam.driverMode == PhysioCam.modeGuessBP) {
				/*messageTxt.setText(String.format("MODE %c CHOOSE SAVED %s WEIGHTS FILE", 
					PhysioCam.driverMode, MLP.makeBPtyp[dsrTyp]) );
				//this.wgtsFilename[dsrTyp] = this.openWgtsFile(); // gets this.wgtsFilename[dsrTyp] */
				wgtsFilename[dsrTyp] = MLP.mlpWgtsFileDefault + MLP.makeBPtyp[dsrTyp] + ".csv";
				convergeSpeedCutoff = 0.999f;
			} // end if (this.driverMode == 2)
			/* ELSE IF 'trainmore' mode 'm', i.e. add new samples to existing best MLP,
			   THEN create layers and nodes from saved weights files. */
			else if (PhysioCam.driverMode == PhysioCam.modeMore) {
				// For each dsrTyp (DBP SBP RBP) open saved MLP weights file
				/*messageTxt.setText(String.format("MODE %c CHOOSE SAVED %s WEIGHTS FILE", 
					PhysioCam.driverMode, MLP.makeBPtyp[dsrTyp]) );
				//this.wgtsFilename[dsrTyp] = this.openWgtsFile(); // gets this.wgtsFilename[dsrTyp] */
				wgtsFilename[dsrTyp] = MLP.mlpWgtsFileDefault + MLP.makeBPtyp[dsrTyp] + ".csv";
				convergeSpeedCutoff = 0.999f;
			} // end if (this.driverMode == 2)

			/* For this.numMLPtoTrain, CREATE a MLP structure.
			 *  Tell MLP how many MLPs to set up and train,
			 *  and max iterations to train each one.
			 * Each MLP starts with random weights.
			 * also saves best MLP to PhysMLPwgtNNNseriesBest.csv  */
			nnMLP[dsrTyp] = new MLP(this, dsrTyp, PhysioCam.driverMode, 
					numMLPtoTrain, maxIter, nodesPerLay, 
					inputArr, numInputVarsUse, outputArr, 
					learnRate, momentum, minErr, 
					wgtsFilename[dsrTyp], convergeSpeedCutoff); //, testRes); 
		} // end LOOP thru dsrTyp
	} // End public void makeMLPs()
	
	public void runMLPs() {	
		nnMLP[0].start(); //start all 3 Threads
		nnMLP[1].start();
		nnMLP[2].start();
		try {
			nnMLP[1].join(); //wait for SBP
		} // wait for the MLP to die 
		catch(InterruptedException iex) { } 
		messageTxt.setText(MLP.testResStr.toString() + String.format("\nMLP Err DBP=%.6f", PhysioCam.overallErrPrev[0]) 
				+ String.format("\n SBP=%.6f", PhysioCam.overallErrPrev[1]) 
				+ String.format(" RBP=%.6f", PhysioCam.overallErrPrev[2]) );
		viewSBPguess.setText(MLP.guessSBPvalStr);
		viewDBPguess.setText(MLP.guessDBPvalStr);
		viewRBPguess.setText(MLP.guessRBPvalStr);
	} // End public void runMLPs()
	
	/** onPause()
     * Called when user pauses the app. 
     * Stop camera preview and release the camera. 
     */
    @Override
    public void onPause()  {
        super.onPause();
        try {
            camera.stopPreview();
            camera.setPreviewCallback(null);
            camera.release();
        } catch (Exception e) {  // ignore: tried to stop non-existent preview
            Log.d("Camera", "Error stopping camera preview: " + e.getMessage());
        }
    }  // End onPause()
    
    /** onResume()
     * Called when the activity is resumed by user. 
     * Reset necessary variables and start new calculation. 
     */
    @Override
    public void onResume()  {
        super.onResume();
        resetVals();
    } // End onResume()
    
    /* onStop - should save data here? */
    @Override
    public void onStop() {
        super.onStop();
        //Log.e(TAG, "-- ON STOP --");
    }  // End onStop()

    @Override
    public void onDestroy() {
        super.onDestroy();
        //Log.e(TAG, "--- ON DESTROY ---");
    }  // End onDestroy()

} // END class PhysioCam
