package my.pphr.PhysioFT_MLP;
/** public class MLP extends Thread 
 * is based on a multi-layer Back Propagation Neural Network
 * Created by Anthony J. Papagelis & Dong Soo Kim
 *  DateCreated:    15 September, 2001
 *  Last Update:    24 October, 2001
               2017-12-14 HL Seldon 2018-02-13 2018-10-26
 * MLP has 2 major functions (methods):
 *    - trainNetwork() - to run iterations of MLP to optimize weights. 
 *          This is to find best MLP. (PhysioCam trains many MLPs, selects best)
 *    - testNetwork(float[] input) - given 1 new input sample, predicts output SBP DBP RBP
 * FORWARD CALCS:
 * Loop thru inputArr[i] of samples. Within sample features are columns j
 * inputArr[i][j] 
 *     -> Layer[0].nodesOutput[j] x Layer[1].nodesWgts[j][k]
 *     -> Layer[1].nodesInputSum[k] 
 *     -> (via outputFcn) Layer[1].nodesOutput[k] x Layer[2].nodesWgts[k][l] 
 *     -> Layer[2].nodesInputSum[l] 
 *     -> (via outputFcn) Layer[2].nodesOutput[l] x Layer[3].nodesWgts[l][m]
 *     -> Layer[3].nodesInputSum[m] 
 *     -> (via outputFcn) Layer[3].nodesOutput[m] 
 *     -> collect Layer[3].nodesOutput[m] in actualOutput[i][m]
 *  End Loop
 *  and compare actualOutput[i][m] with expectedOutput[i][m]
 *  nodesWgts[w][n] arrays have nodes in columns[n]
 *      and weights for each node in rows[w] under the column.
 *      For matrix multiplication.
 *  
 *  LAYER class is included as inner class. 2018-10-26
 *      each LAYER contains 1D arrays of node inputs & outputs, 2D arrays of nodes-weights etc.
*/
import android.os.Environment;
import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
//import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.io.FileNotFoundException;

public class MLP extends Thread {
	public PhysioCam physioParent;
    public static String[] makeBPtyp = {"DBP", "SBP", "RBP"}; // DBP SBP RBP
    public int dsrTypMLP;  // index to makeBPtyp[], cf Driver dsrTyp
    public String bpTyp = "Typ";
    /* mlpMode = mode of operation.  PhysioCam.modeFindMLP - generate random MLPs and find best one.
       PhysioCam.modeGuessBP - test a MLP, load weights from file.
       PhysioCam.modeMore - train an existing MLP more, load weights from file. */
    public char mlpMode = PhysioCam.modeFindMLP;
    public int numMLPtoTry = 1; // number of MLPs to try for 'findmlp' Mode 1
    
    // The user-defined input pattern for a set of samples
    private float inputArr[][]; //points to InputSamples[][] from PhysioCam
    // Use ONLY the first numInVarsToUse variables
    //  should = numNodesPerLayer[0].length as layer[0] contains the input vars to use
    private int numInVarsToUse = 35; // make input parameter
    // User defined learning rate - used for updating the network weights
    private float LearningRate;
    // Users defined momentum - used for updating the network weights
    private float Momentum;
    // The user-defined expected output pattern for a set of input samples
    private float expectedOutput[][]; // points to expectInputSamples[][] from Driver
    public float  actualOutput[][]; // one row per input sample. cols = 0 or 1. Matches expectedOutput[][]
    // Error Function variable that is calculated using the calcOverallError() function
    private float OverallError;
    public float totalErr;
    // The max error per input sample targeted by findmlp or trainmore. e.g. 0.0001
    private float targetError;
    // min error found over series of MLP trials (in one MLP object)
    public float minMLPseriesErr = 10.f;
    
    // Number of layers in the network - includes the input, output and hidden layers
    private int numLayers;
    public int[] numNodesPerLayer; // array of nodes per layer
    public int outLayer;        // pointer to last (output) layer
    public float[] output_nodes; // shortcut to (Layer[Layer.length - 1]).nodesOutput
    public LAYER Layer[];         // the layers
    
    // Number of training sets
    public int numInputSamples;
    // Current training set/sample that is used to train network
    private int sampleN;
    private int sampleNstart, sampleNstep; // sampleN start value, step (forwards, backwards...)
    // Maximum number of Epochs before the traing stops training - user defined
    private long maxNumIteration = 20000;
    /* breakThreshold - control descent. Each 1000th iteration must have 
     * error < breakThreshold * previous error, else no good so break iterations */
    public float breakThreshold = 0.999f;
    public static StringBuffer testResStr; // String to show guess results for each sample to guess
	public String sampleN_DTstr;
	public int sampleN_DTstrIdx, sampleN_DTstrLgt;
	public static String guessSBPvalStr, guessDBPvalStr, guessRBPvalStr; // to pass values from synchronized train() to PhysioCam

    // Saved files of weights
    public String mlpWgtsFileName; // name of file with stored weights
    public static String mlpWgtsFileDefault = /*PhysioCam.dataDir +*/ "PhysFT_Wgt";
    public File mlpWgtsFile;
	public boolean mlpWgtsFileExists = false;
    /* this.mlpWgtsFile - first line
               MLP SN [0]
               layer0 numNodes [1] --this.numNodesPerLayer[0]
               layer1 numNodes [2]
               layer2 numNodes [3]
               layer3 numNodes [4] .. maybe more
               learning rate [5]
               momentum [6]
               num input samples used [7]
               min error achieved by the series [8]
        Following lines are thresholds and weights of nodes in layers*/
    public BufferedReader mlpWgtsFileRead;
    public PrintWriter mlpWgtsFileWrite;
    public int saveFileSN = 0; // SN of best-weights file, allows updating file
    public String wgtsLine;
    public int mlpWgtsFileNumSamples = 0;
    public int mlpWgtsFileSN = 0;
    public float mlpWgtsFileErr = 0.f;
    
    
    /** 
	 * @param physioPar - pointer to parent PhysioCam
     * @param bpTypIdx - index to type of MLP: 0-"DBP", 1-"SBP", 2-"RBP"
     * @param mode - mode of operation. PhysioCam.modeFindMLP-find best MLP. PhysioCam.modeGuessBP-test or guess. PhysioCam.modeMore-train more.
     * @param howManyMLP - number of MLPs to try - see public int numMLPtoTry
     * @param MaxIter - max num iterations to train or test each MLP
     * @param numNodesPerLayer[] - array of num nodes in each layer
     * @param InputSamples[][] - each row is one input vector,
     *                         within row, cols are input dimensions
     * @param numVarsUse - num of input variables / features to use, e.g. first 35 out of 41 vars
     * @param expectOutputSamples[][] - each row matches one input vector
     *                         col for expected val = 1, all others = 0
     * @param LearnRate - learning rate 0-1
     * @param Moment - momentum
     * @param MinError - target total error
     * @param wgtsFileNm - name of a saved weights file, if any, else null.
     *          If there is a weights file, the parameters are taken from line 1
     *          then the weights are loaded from the file. If no weights file,
     *          params are those in the constructor.
     * @param breakOnErr - value passed to nnMLP[this.dsrTypMLP].trainNetwork() to check if err decreasing
     *                      break if err > val * previous err, e.g. 0.999 * previous err
     * @param dTestRes[][] - pointer to testRes[][] in PhysioCam
    */
    public MLP(PhysioCam physioPar, int bpTypIdx, char mode, int howManyMLP, long MaxIter, 
                    int inpNodesPerLayer[], float InputSamples[][], int numVarsUse,
                    float expectOutputSamples[][], float LearnRate, 
                    float Moment, float MinError, String wgtsFileNm, 
                    float breakOnErr) { //, int dTestRes[][]) {
        // Initialize variables
		this.physioParent = physioPar;
        this.dsrTypMLP = bpTypIdx; // index for type of MLP (0-DBP, 1-SBP, 2-RBP)
        this.bpTyp = MLP.makeBPtyp[this.dsrTypMLP]; // label for type of MLP
        this.mlpMode = mode; // mode of operation 0-findmlp 1-test 2-trainMore
        this.numMLPtoTry = howManyMLP; // number of MLPs for 'findmlp' mode
        if (this.numMLPtoTry < 1) this.numMLPtoTry = 1;
        this.maxNumIteration = MaxIter;
        if (this.mlpMode == PhysioCam.modeGuessBP) this.maxNumIteration = 1; // IF test mode, only 1 iteration
        this.numInputSamples = InputSamples.length;
        this.inputArr = InputSamples; // pointer to PhysioCam.inputArr[][]
        this.numInVarsToUse = numVarsUse; // num input vars to use for MLP
        this.expectedOutput = expectOutputSamples; // pointer to PhysioCam.expectOutput[][]
        this.breakThreshold = breakOnErr;
        if ((this.breakThreshold > 1.2f) || (this.breakThreshold < 0.99f)) 
                this.breakThreshold = 0.999f;
        this.minMLPseriesErr = 10.f;
        //this.driverTestRes = dTestRes; // pointer to PhysioCam.testRes[][]
        //PhysioCam.mlpLogWrite.printf("inputArr rows %d cols %d\tNum Input Samples %d\tNum features used %d\n", 
        //    this.inputArr.length, this.inputArr[0].length, this.numInputSamples, this.numInVarsToUse);
        //System.out.printf("inputArr rows %d cols %d\tNum Input Samples %d\tNum features used %d\n", 
        //    this.inputArr.length, this.inputArr[0].length, this.numInputSamples, this.numInVarsToUse);
        
        /* 2. To MAKE LAYERS (this.numNodesPerLayer)
         *   use either the input inpNodesPerLayer
         *   OR, if there is a wgtsFile, get values from its first line. */
        /* IF test mode 1 or train further mode 2
         * AND IF there exists this.mlpWgtsFileName 
               THEN read first line - this.loadMLPwgts(this.mlpWgtsFileName, 1)
                  and get 
               MLP SN [0]
               layer0 numNodes [1] --this.numNodesPerLayer[0]
               layer1 numNodes [2]
               layer2 numNodes [3]
               layer3 numNodes [4] .. maybe more
               learning rate [5]
               momentum [6]
               num input samples used [7]
               min error achieved by the series [8]
           and calc this.numLayers */
        // Check IF there is a mlpWgtsFile
		this.mlpWgtsFileName = MLP.mlpWgtsFileDefault + this.bpTyp + ".csv";
		try {
            this.mlpWgtsFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS), this.mlpWgtsFileName);
            this.mlpWgtsFileExists = this.mlpWgtsFile.exists();
			//mlpWgtsFile = new File(this.mlpWgtsFileName);
        } catch(NullPointerException ioe) {
            this.physioParent.messageTxt.setText("ERR creating mlpWgtsFile " + ioe.getMessage());
        }
		// IF mlpWgtsFile does NOT EXIST, set mlpWgtsFileName = null for safety
		//if (!this.mlpWgtsFile.exists()) this.mlpWgtsFileName = null;
		//if (wgtsFileNm != null) {
        //    this.mlpWgtsFileName = new String(wgtsFileNm); } // name of saved weights file
        //else this.mlpWgtsFileName = null;
		// IF 'test' or 'trainmore' mode, load params from first line of mlpWgtsFile
        if ((this.mlpMode != PhysioCam.modeFindMLP) && (this.mlpWgtsFileExists)) { 
            if ( !this.loadMLPwgts(this.mlpWgtsFileName, 1)) { // IF ERROR loading params from weights file
				this.physioParent.messageTxt.setText("ERR loadMLPWgts line 1");
				this.targetError = MinError;
				this.LearningRate = LearnRate;
				this.Momentum = Moment;
				this.numNodesPerLayer = inpNodesPerLayer;
				this.numLayers = this.numNodesPerLayer.length;
			}
        } // end if want to load params from weights file
        // ELSE (mlpMode = 'findmlp' or there is no weights file), get input params
        else {
            this.targetError = MinError;
            this.LearningRate = LearnRate;
            this.Momentum = Moment;
            this.numNodesPerLayer = inpNodesPerLayer;
            this.numLayers = this.numNodesPerLayer.length;
        }
        this.outLayer = this.numLayers -1; // index of output layer
        
        //PhysioCam.mlpLogWrite.printf("%s Model Layers/Nodes\n", this.bpTyp);
        // Create network layers from array this.numNodesPerLayer
        Layer = new LAYER[this.numLayers];
        /* LAYER constructor will create array of numNodesPerLayer[i] nodes
           Each node will have numNodesPerLayer[i-1] weights. */
        /* Layer[0] = input, so num Nodes = num input vect dimensions
         *  Each node has 'num input vect dimensions' weights */
        Layer[0] = new LAYER(this, this.numNodesPerLayer[0],this.numNodesPerLayer[0]);
        numInVarsToUse = this.numNodesPerLayer[0];
        //PhysioCam.mlpLogWrite.printf("layer\t0\tnumNodes\t%d\tnumWgt/node\t%d\n",
        //    this.numNodesPerLayer[0], this.numNodesPerLayer[0]);
        //System.out.printf("layer\t0\tnumNodes\t%d\tnumWgt/node\t%d\n",
        //    this.numNodesPerLayer[0], this.numNodesPerLayer[0]);
        // Create remaining layers and node arrays and nodes
        for (int i = 1; i < this.numLayers; i++) {
            Layer[i] = new LAYER(this, this.numNodesPerLayer[i],this.numNodesPerLayer[i-1]);
            //PhysioCam.mlpLogWrite.printf("layer\t%d\tnumNodes\t%d\tnumWgt/node\t%d\n", i,
            //    numNodesPerLayer[i], numNodesPerLayer[i-1]);
            //System.out.printf("layer\t%d\tnumNodes\t%d\tnumWgt/node\t%d\n", i,
            //    numNodesPerLayer[i], numNodesPerLayer[i-1]);
        }
        
        /* 3. To make OUTPUT arrays, create actualOutput[][] arrays
         * Layer[last] = final output, so num nodes = num distinct output values */
        actualOutput = new float[numInputSamples][Layer[this.outLayer].numNodes];
		// init testResStr to display test results
		MLP.testResStr = new StringBuffer("");
    } // End MLP constructor

    /** public void trainNetwork()
     * This actually includes 'test' option. Choice is via the this.mlpMode parameter:
     *    f - findmlp - test numMLPtoTry and select one with lowest error
     *    g - test MLP using selected saved weights files. Set numMLPtoTry = 1, maxIterations = 1
     *    m - trainmore - select saved weights files and train for more iterations
     * FLOW:
     * LOOP thru numMLPtoTry, number MLPs to make and train
     *   LOOP thru iterations for each MLP
     *      LOOP thru samples for each iteration
               Load input into Layer[0]
               Calc output
               if 'test' mode, find best match 
               else
               Load output from Layer[this.outLayer] into actualOutput
               Calc output errors
               Do backpropagation
            END sample LOOP
            if 'test' mode, exit loops
            Calc total error
            On every 1000th iteration, break if Error shows insufficient change 
         END iteration LOOP
         IF total error < previous total, save MLP weights to file
       END LOOP thru MLPs to try
     */
    public synchronized void trainNetwork() {
        int j, lyr, nd, ndp1, iterN = 0, n1, n2, backL, backN, backW;
        float layerNodeOutputErr = 0.f;
        int testWinner = 0;
        this.output_nodes = (Layer[this.numLayers - 1]).nodesOutput;
        
        // LOOP thru ALL MLPs to search for best one
        for (int mp = 0; mp < this.numMLPtoTry; mp++) {
            if (PhysioCam.driverMode == PhysioCam.modeFindMLP)
					PhysioCam.overallErrPrev[this.dsrTypMLP] = 999.f; // init previous total error for 'findMLP' mode
            //System.out.printf("%s %s MLP: %d of %d\n", 
            //   MLP.makeBPtyp[this.dsrTypMLP], PhysioCam.driverMode[this.mlpMode], mp, this.numMLPtoTry);
            //PhysioCam.mlpLogWrite.println(MLP.makeBPtyp[this.dsrTypMLP] + " MLP: " + mp + " of " + this.numMLPtoTry);
                
            /* INITIALIZE NODE WEIGHTS (LOAD or GENERATE RANDOM)
             * IF train mode ('findmlp' 0), generate random weights
             * IF test mode 1 or train further mode 2
             * AND IF there is weights file, load node weights from it */
            if ((this.mlpMode != PhysioCam.modeFindMLP) && (this.mlpWgtsFileExists)) { 
                if ( !this.loadMLPwgts(this.mlpWgtsFileName, 2)) { // IF ERROR loading weights from weights file
					this.physioParent.messageTxt.setText("ERR loadMLPWgts weights"); }
				/* IF 'train further mode m' generate random weights around 
                 * existing best weights, but with narrower range. */
                if (this.mlpMode == PhysioCam.modeMore) {
                    this.initRandomWeights(1, 0.25f);
                    //physioParent.messageTxt.setText("Wgts loaded for " + this.numNodesPerLayer.length + " layers");
                } // end if (this.mlpMode == m)
            } // end if mlpMode = g or m
            else { /* ELSE for initial training mode f, no existing wgt file
                      ONLY need wgts for layers 1-outLayer. Layer 0 is input */
                this.initRandomWeights(0, 1.f);
            } // End loading node weights
            
            /*PhysioCam.mlpLogWrite.printf("MLP %s Mode %s Layers=%d inSamp=%d inDim=%d outVals=%d learnR=%.2f momentum=%.2f maxIter=%d\n",
                MLP.makeBPtyp[this.dsrTypMLP], PhysioCam.driverMode[this.mlpMode], this.numLayers, 
                this.numInputSamples, this.numInVarsToUse, this.expectedOutput[0].length,
                this.LearningRate, this.Momentum, this.maxNumIteration);
            System.out.printf("MLP %s Mode %s Layers=%d inSamp=%d inDim=%d outVals=%d learnR=%.2f momentum=%.2f maxIter=%d\n",
                MLP.makeBPtyp[this.dsrTypMLP], PhysioCam.driverMode[this.mlpMode], this.numLayers, 
                this.numInputSamples, this.numInVarsToUse, this.expectedOutput[0].length,
                this.LearningRate, this.Momentum, this.maxNumIteration);
            //PhysioCam.mlpLogWrite.println("iterat\tAveError"); */
               
            // LOOP through ALL ITERATIONS
            for (iterN = 0; iterN < this.maxNumIteration; iterN++) {
                // For each iteration, LOOP through ALL or some SAMPLES (sampleN)
                // CHOOSE DIRECTION of loop, step size
                if (iterN % 1000 == 0) { // for every 1000th, do all samples from start
                    this.sampleNstart = 0;
                    this.sampleNstep = 1;
                } // end if even 1000
                else { // if NOT multiple of 1000
                    if (iterN % 2 == 0) { // even numbers - 0 to end
                        this.sampleNstart = 0;
                        this.sampleNstep = 2; // 1 -all samples, 2 -even samples
                    }
                    else { // odd numbers - reverse - last odd num to 1
                        if ((this.numInputSamples - 1) % 2 == 0)
                            this.sampleNstart = this.numInputSamples - 2;
                        else 
                            this.sampleNstart = this.numInputSamples - 1;
                        this.sampleNstep = -2; // -1 -all samples, -2 -odd samples
                    }
                } // end if NOT even 1000
                for (sampleN = this.sampleNstart; (sampleN >= 0) && (sampleN < this.numInputSamples); sampleN += this.sampleNstep) {
					// For each sample, point Layer[0].layerInput[] to that Input line
                    Layer[0].layerInput = inputArr[sampleN];
                    
                    //this.calcOutput(); // Calc output
                    /* public void calcOutput() {  Calculate the node activations */
                    Layer[0].nodesOutput = Layer[0].layerInput; // Layer[0] is the input
                    for (int calcOutL = 0; calcOutL < this.numLayers; calcOutL++) {
                        if (calcOutL > 0) {
                            //Layer[calcOutL].calcLayerOutput(); // already did layer0
                            /* calcLayerOutput()
                             * For each node[i] in the layer
                             *     Sums the weighted inputs = inputVect[j] * nodeWeight[j][i]
                             *     Calc the Output as a activationFcn function. */
                            //public void calcLayerOutput() 
                            for (int nD = 0; nD < Layer[calcOutL].numNodes; nD++) {
                                Layer[calcOutL].nodesInputSums[nD] = Layer[calcOutL].nodesThreshold[nD];
                                for (int iIN = 0; iIN < Layer[calcOutL].numInputs; iIN++)
                                    Layer[calcOutL].nodesInputSums[nD] += Layer[calcOutL].layerInput[iIN] * Layer[calcOutL].nodesWgts[iIN][nD];
                                Layer[calcOutL].nodesOutput[nD] =  1.f / (1.f + (float)Math.exp(-Layer[calcOutL].nodesInputSums[nD]));
                            }
                        } // End calcOutput() for layers >0
                        // Unless we have reached the last layer,
                        //  nodesOutput from Layer[i] = InputVect of Layer[i+1]
                        if (calcOutL < this.outLayer)
                            Layer[calcOutL+1].layerInput = Layer[calcOutL].nodesOutput;
                    } // end for calcOutL loop thru layers
                    //} // End calcOutput()
                    
                    if (this.mlpMode == PhysioCam.modeGuessBP) { // IF 'test' mode
                    /* get output layer of nodes and look for the node with largest output */
                        testWinner = 0;
                        for (int k=0; k < this.output_nodes.length; k++) {
                            if (this.output_nodes[testWinner] < this.output_nodes[k]) {
                                testWinner = k;
                            } // if
                        } // for k loop through this.output_nodes
                        /* IF BP unknown, i.e. DBP SBP = default values
							Construct sampleN line in MLP.testResStr to show DBP SBP RBP results  */
						if ((this.physioParent.inputDBP[sampleN] == 11) && (this.physioParent.inputSBP[sampleN] == 111)) {
							this.sampleN_DTstr = this.physioParent.sampleDateTime.get(sampleN); // get datetime as string
							this.sampleN_DTstrLgt = this.sampleN_DTstr.length();
							if (MLP.testResStr.indexOf(this.sampleN_DTstr) < 0)
								MLP.testResStr.append("\n" + this.sampleN_DTstr);
							this.sampleN_DTstrIdx = MLP.testResStr.indexOf(this.sampleN_DTstr) + this.sampleN_DTstrLgt;
							if (this.dsrTypMLP == 0) {
								MLP.guessDBPvalStr = String.valueOf(this.physioParent.dbpMin + 2 * testWinner);
								MLP.testResStr.insert(this.sampleN_DTstrIdx, " D=" + MLP.guessDBPvalStr);
							}
							else if (this.dsrTypMLP == 1) { 
								MLP.guessSBPvalStr = String.valueOf(this.physioParent.sbpMin + 2 * testWinner);
								MLP.testResStr.insert(this.sampleN_DTstrIdx, " S=" + MLP.guessSBPvalStr); 
							}
							else if (this.dsrTypMLP == 2) { 
								MLP.guessRBPvalStr = String.valueOf(this.physioParent.rbpMin + 2 * testWinner);
								MLP.testResStr.insert(this.sampleN_DTstrIdx, " R=" + MLP.guessRBPvalStr);
							}
						} // end if DBP = 11 & SBP = 111 default values
                    } // end this.mlpMode == g - testing
                    
                    else {  // IF findmlp or trainmore mode
                        /* For each sample[sampleN], Copy calculated output 
                         * from Layer[last].nodesOutput[] to actualOutput[sampleN][] */
                        for (j = 0; j < Layer[this.outLayer].numNodes; j++)
                            actualOutput[sampleN][j] = Layer[this.outLayer].nodesOutput[j];
                        
                        /* private void calcOutputErrors() 
                        Calculate output error in each output layer node[nd]
                           = (expected[for sampleN] - actual) * actual * (1 - actual)
                        private void calcOutputErrors() {  */
                        // FIRST for output Layer
                        for (nd = 0; nd < Layer[this.outLayer].numNodes; nd++) {
                            Layer[this.outLayer].nodesOutputError[nd] 
                                = (expectedOutput[sampleN][nd]
                                        - Layer[this.outLayer].nodesOutput[nd])
                                  * Layer[this.outLayer].nodesOutput[nd] 
                                  * (1.f - Layer[this.outLayer].nodesOutput[nd]);
                        }
                        // SECOND Calculate signal error for all nodes in the hidden layers
                        // Loop lyr from next-last layer back to first hidden layer
                        for (lyr = this.numLayers-2; lyr > 0; lyr--) {
                            // Loop nd thru nodes in Layer lyr
                            for (nd = 0; nd < Layer[lyr].numNodes; nd++) {
                                layerNodeOutputErr = 0.f;
                                // Loop ndp1 thru nodes in Layer lyr+1
                                // sum node[ndp1].weight[nd] * node[ndp1].error
                                for (ndp1 = 0; ndp1 < Layer[lyr+1].numNodes; ndp1++) {
                                    layerNodeOutputErr = layerNodeOutputErr + Layer[lyr+1].nodesWgts[nd][ndp1] * 
                                        Layer[lyr+1].nodesOutputError[ndp1]; }
                                // error in Layer[i].node[j] = node[j]output * (1 - node[j]output) * sum from layer i+1
                                Layer[lyr].nodesOutputError[nd]
                                    = Layer[lyr].nodesOutput[nd] * (1.f - Layer[lyr].nodesOutput[nd]) * layerNodeOutputErr;
                            } // end for nd loop thru nodes in Layer lyr
                        } // end for lyr loop backwards thru layers
                        //}  // End private void calcOutputErrors()
                        
                        /* private void BackPropagateWgt() {  */
                        // Update Weights
                        for (backL = this.outLayer; backL > 0; backL--) {
                            for (backN = 0; backN < Layer[backL].numNodes; backN++) {
                                // Calculate Bias weight difference to node backN
                                Layer[backL].nodesThresholdDiff[backN]
                                    = LearningRate * Layer[backL].nodesOutputError[backN]
                                      + Momentum * Layer[backL].nodesThresholdDiff[backN];
                
                                // Update Bias weight (threshold) to node backN
                                Layer[backL].nodesThreshold[backN] += Layer[backL].nodesThresholdDiff[backN];
                
                                // Update Weights
                                for (backW = 0; backW < Layer[backL].numInputs; backW++) {
                                    // Calculate weight difference between node backN and backW
                                    // Layer[backL] inputs = Layer[backL-1].nodesOutput[]
                                    Layer[backL].nodesWgtDiff[backW][backN] = 
                                        LearningRate * 
                                        Layer[backL].nodesOutputError[backN] * Layer[backL-1].nodesOutput[backW]
                                        + Momentum * Layer[backL].nodesWgtDiff[backW][backN];
                
                                    // Update weight between node backN and backW
                                    Layer[backL].nodesWgts[backW][backN] += Layer[backL].nodesWgtDiff[backW][backN];
                                } // end for backW loop through weights in each node
                            } // end for backN loop through nodes in each layer
                        } // end for backL loop through layers
                        //} // End private void BackPropagateWgt()
                    } // End else findmlp or trainmore mode
    
                } // end for sampleN loop
				
				// if 'test' mode, show test results, then finish
                if (this.mlpMode == PhysioCam.modeGuessBP) {
					this.physioParent.messageTxt.setText(MLP.testResStr);
					break;
				} 
                /* private void calcOverallError()
                 * LOOP through all input samples
                 *    for each input sample LOOP through the nodes in the output layer
                 *       for each node, add (expectedOutput - actualOutput)^2 / 2.
                 * Divide total error by num input samples to get average error per sample.
                private void calcOverallError() { */
                this.OverallError = 0.f;
                for (sampleN = 0; sampleN < numInputSamples; sampleN++) {
                    for (nd = 0; nd < Layer[this.outLayer].numNodes; nd++) {
                        this.OverallError += (float)(0.5 * ( Math.pow(expectedOutput[sampleN][nd] - actualOutput[sampleN][nd], 2.) ));
                    }
                } // end for sampleN loop
                if (numInputSamples > 0) this.OverallError /= numInputSamples;
                //} // end calcOverallError()
                
                /* On every 1000th iteration, check if OverallError decreasing fast enough.
                   Use iter 1000 as baseline. Generous empirical cutoff thresholds.
                   IF not decreasing enough, BREAK out of iteration loop */
                if (iterN % 1000 == 0) {
                    /*PhysioCam.mlpLogWrite.printf("%s\t%s\t%3d of %3d\t%6d\t%f\t%f\n", 
                       MLP.makeBPtyp[this.dsrTypMLP], PhysioCam.driverMode[this.mlpMode], mp, this.numMLPtoTry, 
                        iterN, this.OverallError, PhysioCam.overallErrPrev[this.dsrTypMLP]); */
                    /*this.physioParent.messageTxt.setText(String.format("%s\t%3d of %3d\t%6d\t%f\t%f\n", 
                       MLP.makeBPtyp[this.dsrTypMLP], mp, this.numMLPtoTry, 
                        iterN, this.OverallError, PhysioCam.overallErrPrev[this.dsrTypMLP]) ); */
                    if (this.OverallError > this.breakThreshold * PhysioCam.overallErrPrev[this.dsrTypMLP]) break;
                    PhysioCam.overallErrPrev[this.dsrTypMLP] = this.OverallError;
                } // end if iterN % 1000 = 0
            } // end iterN loop through iterations
            
            if (this.mlpMode == PhysioCam.modeGuessBP) break; // if 'test' mode, finished
            // if total error < best up to now, save these MLP weights
            if (this.OverallError < this.minMLPseriesErr) {
                this.saveMLPwgtsFile();
                this.minMLPseriesErr = this.OverallError;
            }
			if (this.OverallError < PhysioCam.overallErrPrev[this.dsrTypMLP]) {
                PhysioCam.overallErrPrev[this.dsrTypMLP] = this.OverallError; }
            if (this.OverallError < this.targetError) break; // may stop if reached < 0.0001
        } // End mp LOOP thru numMLPtoTry
        
    } // End public void trainNetwork()

    
    /** void saveMLPwgtsFile()
     * dumps nodes/layer,learn rate,momentum to a fixed csv file
     * dumps the weights to the csv file. */
    public void saveMLPwgtsFile() { 
        int i, j, w;
		// the Constructor creates mlpWgtsFileName and mlpWgtsFile
        //this.mlpWgtsFileName = MLP.mlpWgtsFileDefault + this.bpTyp + ".csv";
                    //       + "_S" + String.format("%4d", this.numInputSamples) 
                    //       + "E" + String.format("%1.6f", this.OverallError)
        this.physioParent.messageTxt.setText(String.format("%s OverallErr=%.6f this.minMLPseriesErr=%.6f -saveNetwork\n",
                    MLP.makeBPtyp[this.dsrTypMLP], this.OverallError, this.minMLPseriesErr) );  
        try {
            //mlpWgtsFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS), this.mlpWgtsFileName);
            //mlpWgtsFile = new File(this.mlpWgtsFileName);
            mlpWgtsFileWrite = new PrintWriter((mlpWgtsFile));
        } catch(IOException ioe) {
            this.physioParent.messageTxt.setText("ERR saveNetToFile " + ioe.getMessage());
            return;
        }
        // FIRST LINE = Serial Num, num Nodes in each layer
        mlpWgtsFileWrite.printf("%d,", this.saveFileSN + 1);
        for (i = 0; i < numLayers; i++) {
            mlpWgtsFileWrite.printf("%d,", Layer[i].numNodes);
        }
        // also learn rate, momentum, average network error per sample
        mlpWgtsFileWrite.printf("%.7f,%.7f,%d,%.7f\n", this.LearningRate,
                this.Momentum, this.numInputSamples, this.OverallError);
        // following lines are the Layers, Nodes in Layer, and Weights in each Node
        // Layer[0] weights are NOT USED, so skip them and start with Layer[1]
        for (i = 1; i < numLayers; i++) {
            for (j = 0; j < Layer[i].numNodes; j++) {
                // write layer and node num, node Threshold at start of line
                mlpWgtsFileWrite.printf("%d,%d,%.7f,", i, j, Layer[i].nodesThreshold[j]);
                // write weights[w] for node [j] in Layer[i] on rest of line
                for (w = 0; w < Layer[i].numInputs; w++) {
                    mlpWgtsFileWrite.printf("%.7f,", Layer[i].nodesWgts[w][j]);
                } // end for w loop through Weights in node j in layer i
                mlpWgtsFileWrite.println();
            } // end for j loop through nodes in a layer
        } // end for i loop through layers
        mlpWgtsFileWrite.flush();
        mlpWgtsFileWrite.close(); 
    } // End void saveMLPwgtsFile()
    
    /** public boolean loadMLPwgts(String fromFile, int processLn)
     *  Reads a text file of saved MLP weights.
     *  First line: nodes_Lay0,nodes_Lay1,nodes_Lay2,nodes_Lay3,learnRate,momentum,totalError
     *  Subsequent lines:
     *    LayerNum,NodeNum,NodeThreshold,Wgt1,Wgt2,....,WgtN
     *  Loops through ArrayList
     *     Reads input line
     *     splits into parts
     *     puts NodeThreshold into net[this.dsrTypMLP].Layer[LayerNum].nodeArray[nodeNum].Threshold
     *     puts Wgts into net[this.dsrTypMLP].Layer[LayerNum].nodeArray[nodeNum].Weight[]
     *  @param processLn - 1=process only first line, 2=process all other lines
     *  @return boolean - false = IO error, true = OK
     */
    public boolean loadMLPwgts(String fromFile, int processLn) {
        // the Constructor creates mlpWgtsFileName and mlpWgtsFile
		//this.mlpWgtsFileName = MLP.mlpWgtsFileDefault + this.bpTyp + ".csv"; // GET REAL NAME OF WGTS FILE in
                     //      + "_S" + String.format("%4d", this.numInputSamples) 
                     //      + "E" + String.format("%1.6f", this.OverallError)
		String[] lineParts;
        int lineCnt = 0, layN, nodeN;
        // open file and process first line
        if (processLn == 1) {
            try {
				//this.mlpWgtsFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS), 
				//this.mlpWgtsFileName); //new File(fromFile);
                mlpWgtsFileRead = new BufferedReader(new FileReader(this.mlpWgtsFile));
                this.wgtsLine = mlpWgtsFileRead.readLine(); // READ FIRST LINE
                //System.out.println(this.wgtsLine);
                lineParts = this.wgtsLine.split(",");
                /* First line is
                 * SN,numNode0,numNode1,numNode2,numNode3,LearnR,Momentum,numSamples,Error
                 */
                this.mlpWgtsFileErr = Float.parseFloat(lineParts[lineParts.length - 1]);
                PhysioCam.dsrTotalErr[this.dsrTypMLP] = this.totalErr = this.minMLPseriesErr = this.mlpWgtsFileErr;
				if (PhysioCam.driverMode != PhysioCam.modeFindMLP) { // use existing error for guessBP and trainMore modes
					PhysioCam.overallErrPrev[this.dsrTypMLP] = this.mlpWgtsFileErr;	}
                this.mlpWgtsFileNumSamples = Integer.parseInt(lineParts[lineParts.length - 2]);
                this.Momentum = Float.parseFloat(lineParts[lineParts.length - 3]);
                this.LearningRate = Float.parseFloat(lineParts[lineParts.length - 4]);
                this.mlpWgtsFileSN = Integer.parseInt(lineParts[0]);
                this.numLayers = lineParts.length - 4 - 1; // should = 9-4-1
                this.numNodesPerLayer = new int[this.numLayers];
                for (int i = 0; i < this.numLayers; i++) {
                    this.numNodesPerLayer[i] = Integer.parseInt(lineParts[i + 1]);
                    //System.out.printf(" Layer=%d Nodes=%d", i, this.numNodesPerLayer[i]);
                }
                mlpWgtsFileRead.close();
                return true;
            } // end try
            catch(FileNotFoundException fne) {
                this.physioParent.messageTxt.setText(fne.getMessage());
                return false; }
            catch(NumberFormatException nne) {
                this.physioParent.messageTxt.setText("ERR: first line not num -" + nne.getMessage());
                return false;   }
            catch (IOException ioe) {
                this.physioParent.messageTxt.setText("ERR reading MLPweights");
                return false;  }  
        } // end process first line, return true or false
        
        // NOT FIRST LINE. process lines with weights
        try {
            mlpWgtsFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DOWNLOADS), 
				this.mlpWgtsFileName); //new File(fromFile);
            mlpWgtsFileRead = new BufferedReader(new FileReader(mlpWgtsFile));
            this.wgtsLine = mlpWgtsFileRead.readLine(); // READ FIRST LINE. skip it
            while ((this.wgtsLine = mlpWgtsFileRead.readLine()) != null) {
                lineParts = this.wgtsLine.split(",");
                layN = Integer.parseInt(lineParts[0]);
                nodeN = Integer.parseInt(lineParts[1]);
                this.Layer[layN].nodesThreshold[nodeN] = Float.parseFloat(lineParts[2]);
                for (int ff = 3; ff < lineParts.length; ff++) {
                    this.Layer[layN].nodesWgts[ff-3][nodeN] = Float.parseFloat(lineParts[ff]);
                } // end for ff loop parsing one line
            } // end while read input lines
            mlpWgtsFileRead.close();
        } // end try
        catch(FileNotFoundException fne) {
            this.physioParent.messageTxt.setText(fne.getMessage());
            return false; }
        catch(NumberFormatException nne) {
            this.physioParent.messageTxt.setText("ERR: wgts line not num -" + nne.getMessage());
            return false;  }
        catch (IOException ioe) {
            this.physioParent.messageTxt.setText("ERR reading MLPweights");
            return false;  }  
        return true;
    } // End loadMLPwgts(..)
    
    /** initRandomWeights(int newMore, float rng)
     * LOOPs thru layers, from first hidden to outLayer
     *    LOOPs thru nodes in each layer
     *       LOOPs thru weights in each node in each layer
     *          creates a randomly generated number range of 2*rng,
     *          adds it to the existing Threshold of the current node.
	 *			IF new weight, sets weight = random num
	 *			ELSE there is already weight, adds the random num to old weight
                Sets Threshold Diff =  WeightDiff[] = 0 
                   so that the Momentum term can work during the first iteration
     * @param newMore - indicate whether new weights (0) or add to existing weights (1)
     * @param rng - to adjust the range of generated weights. 
     *              1. = range 2 (-1 to 1); 0.5 = range 1 (-0.5 to 0.5); etc
     *              rangeF = rng
     *              OR rangeF = Xavier linear distrib in range +-sqrt(6/(numIn [+ numOut])) */
    public void initRandomWeights(int newMore, float rng) {
        int inL, inN, inW; // indices for Layer Node Wgt
        float rangeF, tmpWgt;
        for (inL = 1; inL < this.numLayers; inL++) {
            // rangeF = Math.sqrt(6. / (Layer[inL].numInputs + Layer[inL].numNodes)); linear version of Xavier
            rangeF = rng; // * rangeF
            for (inN = 0; inN < Layer[inL].numNodes; inN++) {  // Loop thru nodes columns in layer
                Layer[inL].nodesThreshold[inN] += 2. * (Math.random() - 0.5) * rangeF;           
                Layer[inL].nodesThresholdDiff[inN] = 0.f;     
                for (inW = 0; inW < Layer[inL].numInputs; inW++) { // loop thru input rows in node
                    tmpWgt = 2.f * (float)(Math.random() - 0.5f) * rangeF; // range of wgts
                    if (newMore == 0) Layer[inL].nodesWgts[inW][inN] = tmpWgt; // new wgt
                    else Layer[inL].nodesWgts[inW][inN] += tmpWgt; // add to existing wgt
                    Layer[inL].nodesWgtDiff[inW][inN] = 0.f;
                }
            } // end inN loop thru nodes in layer inL
        } // End inL loop thru layers
    }  // End private void initRandomWeights()
    
    /** public void run()
     * needed to implement threading. */
    public void run() {
       this.trainNetwork(); 
    }  // End run() 

    // to notify the network to stop training.
    //public void kill() { die = true; }
    
/** The following java code is based on a multi-layer 
 * Back Propagation Neural Network Class (BackPropagation.class)
 * Created by Anthony J. Papagelis & Dong Soo Kim
 *  DateCreated:    15 September, 2001
 *  Last Update:    14 October, 2001 
                 2017-12-15, 2018-10-24, 2019-01-31 HL Seldon
*/
class LAYER {
    private MLP parentMLP = null;
    public int numNodes, numInputs;
    // Vector of inputs signals from previous layer to the current layer
    public float layerInput[]; // pointer to output of previous layer
    // 2D array. Inputs as rows, nodes as columns. Multiply times layerInput 
    public float nodesWgts[][]; 
    // nodesInputSums[]. dimension numNodes. Results of layerInput[]*nodesWgts[]
    public float nodesInputSums[];  
    // nodesOutput[]. dimension numNodes. Result of outputFcn on nodesInputSums
    public float nodesOutput[];
    // nodesOutputError[]. dimension numNodes. Diff between expected - actual output
    public float nodesOutputError[];
    // nodesThreshold[]. dimension numNodes
    public float nodesThreshold[];
    // nodesThresholdDiff[]. dimension numNodes. Diff from one iteration to next.
    public float nodesThresholdDiff[];
    // nodesWgtDiff[][]. Diff in weights from one iteration to next
    public float nodesWgtDiff[][];
    
    /** Constructor 
     *  @param parent - pointer to parent MLP
     *  @param numNodes2 - number nodes in this layer
     *  @param numInputs2 - number of input vector dimensions,
     *        or nodes in next-lower layer */
    public LAYER (MLP parent, int numNodes2, int numInputs2) {
        this.parentMLP = parent;
        this.numNodes = numNodes2;
        this.numInputs = numInputs2;
        // create arrays
        layerInput = new float[numInputs]; // OR point to prev layer output
        nodesWgts = new float[numInputs][numNodes];
        nodesInputSums = new float[numNodes];
        nodesOutput = new float[numNodes];
        nodesOutputError = new float[numNodes];
        nodesThreshold = new float[numNodes];
        nodesThresholdDiff = new float[numNodes];
        nodesWgtDiff = new float[numInputs][numNodes];
    }  // End constructor
} // END class LAYER

} // END class MLP
