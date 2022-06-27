/*
This macro was writte by Damon Poburko in the Department of Biomedical Physiology and Kinesiology, Simon Fraser University, Canada ca. 2019. 

Please feel free to use and adapt this code to suit your needs. If you use it to contribute to future publications, please cite:

KE Simpson, S Faizi, R Venkateshappa, M Yip, R Johal, D Poburko, YM Cheng, D Hunter, E Lin, GF Tibbits, Thomas W. Claydon.(2022) CRISPR-Cas9-Mediated Precise Knock-in Edits in Zebrafish Hearts. Journal of Visualized Experiments

*/


getDimensions(width, height, channels, slices, frames);

framesLast = frames;
var doFolder = true;
var	doPlot = true;
var batchMode = true;
var plotROIs =  "0"; //comma separated list of ROIs to plot
var	pkFitMaxTolerance = 0.03;
var frameRate = 8;   // frames per second
var	maxFFTPeriod = 15; //(frames?)
var	FFTPeakToleranceFraction = 0.01;
var	returnArrayLabels = newArray("nPeaks","tPk2PkMean", "tPk2PkSD", "pkHeightMean", "pkHeightSD","baselineMean", "baselineSD","fwhmMean", "fwhmSD", "decayConstMean", "decayConstSD", "gaussr2Mean", "gaussr2SD","expr2Mmean", "expr2SD","preferredModelMean", "perMain", "amplMain","per2", "ampl2","tPeak0","tPeak1"); 
var	nParameters = returnArrayLabels.length;
var	returnArray = newArray(nParameters);
var img0 = "";
var imgAnalyzed = "";
var totalROIs = 0;
var nFrames = 0;
var numNN = 1;
var minPksAnalyzed = 1;
 var framesFirst = 1;
var traceHolder = newArray(1);
var traceStarts = newArray(1);
var verticiesHolder = newArray(1);
var verticiesStarts = newArray(1);
var roiName = "";
var smoothFactor = 1;

//Exit if ROIs not selected
if (isOpen("ROI Manager")==false) exit("please select an ROI on the heart, /n and add it to the ROI manager (T) /n to measure heart rate");
if (roiManager("count")<1) exit("please select an ROI on the heart, /n and add it to the ROI manager (T) /n to measure heart rate");


/*
 * *Zebra fish HR analysis:*
Make dialog
- DONE frame rate (default with over ride) (assume 18 fps)
    - look for time info in file as default FR
    - opton to restrict to a range of frames
- DONE option to adjust peak height for detection
- PS analysis - try analyzing normalized to % max/min
- DONE span to analyze (frames or time with toggle)
- save ROIs, save csv output
    - report name, #frames analyzed, and maybe range used
    - look into options to save outputs (trace values and peaks)
- note # of intervals that are >1.5x mean period
    - create flag or report
 */

//User Dialog box
Dialog.create("zFish Heart Analyzer");
Dialog.addNumber("cameras frame rate (FPS)", parseInt(call("ij.Prefs.get", "dialogDefaults.frameRate", "8")),0,5," Hz");
Dialog.addNumber("Analysis Range: Start frame", framesFirst);
Dialog.addNumber("Analysis Range: Last frame", framesLast);
Dialog.addNumber("Peak detection: peak height tolerance", parseFloat(call("ij.Prefs.get", "dialogDefaults.pkFitMaxTolerance", "0.03")),3,5," dF/Fo");
Dialog.addCheckbox("analyze normalized power spectrum", call("ij.Prefs.get", "dialogDefaults.normPS", "1"));
Dialog.addNumber("Power Spectrum Analysis: min HR (BPM)", parseInt(call("ij.Prefs.get", "dialogDefaults.minHR", 60)));
Dialog.addNumber("Power Spectrum Analysis: peak height tolerance", parseFloat(call("ij.Prefs.get", "dialogDefaults.FFTPeakToleranceFraction", 0.01)));
Dialog.addNumber("P-S smooth factor", parseFloat(call("ij.Prefs.get", "dialogDefaults.smoothFactor", 1)));
Dialog.addCheckbox("show trace plots", call("ij.Prefs.get", "dialogDefaults.doPlot", "1"));
Dialog.show();

frameRate= Dialog.getNumber();	call("ij.Prefs.set", "dialogDefaults.frameRate", frameRate);
framesFirst = Dialog.getNumber();	call("ij.Prefs.set", "dialogDefaults.framesFirst", framesFirst);
framesLast = Dialog.getNumber();	call("ij.Prefs.set", "dialogDefaults.framesLast", framesLast);
pkFitMaxTolerance= Dialog.getNumber();	call("ij.Prefs.set", "dialogDefaults.pkFitMaxTolerance", pkFitMaxTolerance);
normPS = Dialog.getCheckbox();	call("ij.Prefs.set", "dialogDefaults.normPS", normPS);
minHR = Dialog.getNumber();	call("ij.Prefs.set", "dialogDefaults.minHR", minHR);
FFTPeakToleranceFraction = Dialog.getNumber();	call("ij.Prefs.set", "dialogDefaults.FFTPeakToleranceFraction", FFTPeakToleranceFraction);
smoothFactor = Dialog.getNumber();	call("ij.Prefs.set", "dialogDefaults.smoothFactor", smoothFactor);
doPlot = Dialog.getCheckbox();	call("ij.Prefs.set", "dialogDefaults.doPlot", doPlot);

//maxFFTPeriod = 1/(minHR*60);     // 1/BPM * 1 min/60 sec

t0 = getTime;

//If subregion of trace used, duplicate that region and define a new name for that image. Otherwise, analyze the original image
img0 = getTitle;
imgAnalyzed = img0;
duplicated = false;
if (( framesFirst!=1) || (framesLast!=frames)) {
	if (lastIndexOf(img0,".")>0) imgAnalyzed = substring(img0,1,lastIndexOf(img0,"."))+"_"+framesFirst+"-"+framesLast;
	if (lastIndexOf(img0,".")<0) imgAnalyzed = img0+"_"+framesFirst+"-"+framesLast;
	run("Select None");
	run("Duplicate...", "title="+imgAnalyzed+" duplicate range="+framesFirst+"-"+framesLast+" use");
	duplicated = true;
}

// Loop the analyzis through all ROIs in the ROI manager
for (a=0; a<roiManager("count");a++) {
	selectWindow(imgAnalyzed);
	roiManager("select", a);
	roiName = Roi.getName;
	//This directed the analysis for the current ROI to the funciton "findPeaks"
	maxFFTPeriod = 60/minHR; //convert minimum HR to detect to seconds/beat
	r = findPeaks(a,pkFitMaxTolerance,maxFFTPeriod,FFTPeakToleranceFraction,returnArray, doPlot, plotROIs, smoothFactor );
}

//SAVE RESULTS, EXTRACTED TRACE WITH PEAKS AND FFT ANALYSIS TO SEPARATE CSV FILES

if (duplicated==true){
	close(imgAnalyzed);
}

tElapsed = (getTime - t0)/1000;
print("\\Update0: run time " +tElapsed + " s");
// END OF PRIMARY CODE BLOCK ===========================================================================================================




// FUNCTIONS ==========================================================================================================================
function findPeaks(roi,pkFitMaxTolerance,maxFFTPeriod,FFTPeakToleranceFraction,returnArray, doPlot, plotROIs, smoothFactor ) {	

	doPlot = doPlot;
	pkFitMaxTolerance = pkFitMaxTolerance;
	maxFFTPeriod = maxFFTPeriod;
	FFTPeakToleranceFraction = FFTPeakToleranceFraction;

	run("Plot Z-axis Profile");
	plot = getTitle();
	Plot.getValues(x, Ys);

	//Part1: Fit individual peaks
	//===============================================================================
	//===============================================================================
	len = Ys.length;

	// check the frame rate encoded by timeStamps in the metaData. This will be evident by FIJI/ImageJ 
	// labelling the X-axis on the plots as something other than integer values.
	//If this is true, then assume that the chosen frame rate can be overridden.

	t = Array.copy(x);
	oIntervals = newArray(1);
	for (k=0; k<x.length;k++) {
		t[k] = x[k]/frameRate;
		if (k>0){
			oIntervals = Array.concat(oIntervals,x[k]-x[k-1]);	
		}
	}
	oIntervals = Array.slice(oIntervals,1,oIntervals.length);
	Array.getStatistics(oIntervals,min,max, meanInterval,sdInterval);
	if ( round(meanInterval)!=1) {
		print("\\Update4: detected mean frame interval "+meanInterval+" s. Overriding chosen frame rate of "+ frameRate);
		frameRate = 1/meanInterval;
		t =x;
	}
	
	selectWindow(plot);
	close;
	
	Array.getStatistics(Ys, YsMin, YsMax, YsMean, YsStdDev); 
	Array.getStatistics(t, tMin, tMax, tMean, tSD); 
	rawRange = YsMax-YsMin;
	rawTolerance = YsMin*pkFitMaxTolerance;

	// find Maxima in trace
	rawMaxLocs = Array.findMaxima(Ys,rawTolerance); 	
	xRawMaxima = newArray(rawMaxLocs.length);
	yRawMaxima = newArray(rawMaxLocs.length);
	for (jj= 0; jj < rawMaxLocs.length; jj++){
		xRawMaxima[jj]= t[rawMaxLocs[jj]];
		yRawMaxima[jj] = Ys[rawMaxLocs[jj]];
  	}
  	
  	//find Minima in trace
  	rawMinLocs = Array.findMinima(Ys,rawTolerance);
	xRawMinima = newArray(rawMinLocs.length);
	yRawMinima = newArray(rawMinLocs.length);
	for (jj= 0; jj < rawMinLocs.length; jj++){
		xRawMinima[jj]= t[rawMinLocs[jj]];
		yRawMinima[jj] = Ys[rawMinLocs[jj]];
  	}

	//sort maxima by X
	rankPosArr = Array.rankPositions(xRawMaxima);
	ranks = Array.rankPositions(rankPosArr);
	xTemp = Array.copy(ranks);
	yTemp = Array.copy(ranks);
	for (k=0;k<xRawMaxima.length;k++) {
		xTemp[ranks[k]] = xRawMaxima[k];
		yTemp[ranks[k]] = yRawMaxima[k];
	}
	xRawMaxima = xTemp;
	yRawMaxima = yTemp;

	
	//sort minima by X
	rankPosArr = Array.rankPositions(xRawMinima);
	ranks = Array.rankPositions(rankPosArr);
	
	xTemp = Array.copy(ranks);
	yTemp = Array.copy(ranks);
	for (k=0;k<xRawMinima.length;k++) {
		xTemp[ranks[k]] = xRawMinima[k];
		yTemp[ranks[k]] = yRawMinima[k];
	}
	xRawMinima = xTemp;
	yRawMinima = yTemp;

	if (doPlot == true) {
		Plot.create("Raw trace: roi "+roi, "time (s)", "AFU", t, Ys);
		Plot.setColor("red","red");
	  	Plot.setLineWidth(5);
		if (rawMaxLocs.length != 0) Plot.add("circles",xRawMaxima,yRawMaxima);
		Plot.setLineWidth(1);
		if (rawMaxLocs.length != 0) Plot.add("line",xRawMaxima,yRawMaxima);
		Plot.setColor("blue","blue");
	  	Plot.setLineWidth(5);
		if (rawMinLocs.length != 0) Plot.add("circles",xRawMinima,yRawMinima);
		Plot.setLineWidth(1);
		if (rawMinLocs.length != 0) Plot.add("line",xRawMinima,yRawMinima);
		Plot.setColor("black");
		Plot.setLineWidth(1);
		Plot.add("line",t,Ys);
		Plot.setFontSize(16);
		Plot.show();
	}

	nExp =0;
	nGauss = 0;
	tPeak0 = -1;
	tPeak1 = -1;
	nPeaks = xRawMinima.length-1;

	if (nPeaks>minPksAnalyzed) {	

			pkHeightMean = -1;
		pkHeightSD = -1;


		firstMaximaUsed = 0;
		if (xRawMaxima[0] < xRawMinima[0]) firstMaximaUsed = 1;
	
		print("\\Update5: nPeaks = " + nPeaks);
		if (xRawMaxima[xRawMaxima.length-1] < xRawMinima[xRawMinima.length-1]) nPeaks = xRawMinima.length-2;
/*
		//dataHolders for each peak
		preferredModel = newArray(nPeaks);
		Array.fill(preferredModel,-1);
		fitParam1 = Array.copy(preferredModel);
		fitParam2 = Array.copy(preferredModel);
		fitParam3 = Array.copy(preferredModel);
		fitParam4 = Array.copy(preferredModel);
		fitParam5 = Array.copy(preferredModel);
		fitParam6 = Array.copy(preferredModel);
		Array.fill(fitParam6,1);
		fitParam7 = Array.copy(preferredModel);
		expr2 = Array.copy(preferredModel);
		gaussr2 = Array.copy(preferredModel);

*/		
		pkHeight =  newArray(nPeaks);
		
		tPk2Pk =  newArray(nPeaks);

		for (p = 0; p<nPeaks; p++) {
	
			//find time indices of peaks
			for (tt= 0; tt< t.length; tt++) {
				if(t[tt]<=xRawMinima[p]) tIndex1 = tt;
				if(t[tt]<=xRawMinima[p+1]) tIndex2 = tt;
				if(t[tt]<=xRawMaxima[p+firstMaximaUsed]) tIndexPk = tt;
			}
			if (p==0) tPeak0 = xRawMaxima[p+firstMaximaUsed];
			if (p==1) tPeak1 = xRawMaxima[p+firstMaximaUsed];
			
			if (p>0) {
				tPk2Pk[p] = xRawMaxima[p+firstMaximaUsed] - tPkPrevious;				
			}
			tPkPrevious = xRawMaxima[p+firstMaximaUsed];
				
			pkHeight[p] =  Ys[tIndexPk] - (Ys[tIndex1]+Ys[tIndex2])/2;
	
		}
		
		if (tPk2Pk.length>1) { 
			tPk2Pk = Array.trim(tPk2Pk,xRawMaxima.length-1);
			Array.getStatistics(tPk2Pk,min, max, tPk2PkMean, tPk2PkSD);
			print("\\Update3: tPk2PkMean " +tPk2PkMean);			
		} else {
			tPk2PkMean = -1;
			tPk2PkSD = -1;
		}

		// create some kind of correction by slicing out periods greater than 1.5x mean perdiod, and recalculating
		// mean period of those peaks vs the peaks that are less than 1.5 x mean
		tPk2PkNorm = newArray(1);
		tPk2PkHi = newArray(1);
		
		if (tPk2Pk.length>5) { 
			for (q=0;q<tPk2Pk.length;q++) {
				if (tPk2Pk[q] <= 1.5*tPk2PkMean) tPk2PkNorm= Array.concat(tPk2PkNorm,tPk2Pk[q]);
				if (tPk2Pk[q] > 1.5*tPk2PkMean) tPk2PkHi= Array.concat(tPk2PkHi,tPk2Pk[q]);
			}
			tPk2PkNorm	= Array.slice(tPk2PkNorm,1,tPk2PkNorm.length);
			tPk2PkHi	= Array.slice(tPk2PkHi,1,tPk2PkHi.length);
			Array.getStatistics(tPk2PkNorm,mintPk2PkNorm, maxtPk2PkNorm, meantPk2PkNorm, sdtPk2PkNorm);
			Array.getStatistics(tPk2PkHi,mintPk2PkHi, maxtPk2PkHi, meantPk2PkHi, sdtPk2PkHi);
		} else {
			meantPk2PkNorm = -1;
			meantPk2PkHi = -1;
		}

		Array.getStatistics(pkHeight,min, max, pkHeightMean, pkHeightSD);
		//Array.getStatistics(preferredModel,min, max, preferredModelMean, preferredModelSD);
		Array.getStatistics(yRawMinima,min,max,baselineMean,baselineSD);

			
	Ysqrd = newArray(lengthOf(Ys));  
	for (i=0; i<lengthOf(Ys); i++) {
		Ysqrd[i] = (Ys[i]-YsMin)*(Ys[i]-YsMin);	
	}
	Array.getStatistics(Ysqrd, YsqrdMin, YsqrdMax, YsqrdMean, YsqrdStdDev); 
	rawRMS = sqrt(YsqrdMean);

	//Part2: FFT analysis of temporal plot 
	//===============================================================================
	//===============================================================================

	Array.getStatistics(Ys,Ysmin, Ysmax, Ysmean,YsSD);
	normYs = Array.copy(Ys);

	if (normPS == true) {
		for (j=0;j<normYs.length;j++) {
			normYs[j] = (normYs[j] - Ysmin)/(Ysmax - Ysmin);
			}
		YsForFFY = normYs;	
	} else {
		YsForFFY = Ys;	
	}

	//frequ=len/frameRate;       //cycles per array length
	windowType="None"; //None, Hamming, Hann or Flattop
	y = Array.fourier(YsForFFY);
 	f = newArray(lengthOf(y)); // f will be a fraction of the length of the trace in frames, 
 	fInHz = newArray(lengthOf(y));
 	sPerFrame = frameRate;
	for (i=0; i<lengthOf(y); i++) {
		a = lengthOf(y)-i;
		//f[i] = 1- (a/ lengthOf(y));
		f[i] = (sPerFrame)*( 1- (a/ lengthOf(y)))/2; //seems to capture accurate frequnce for given frame rate
		//f[i] = f[i] * lengthOf(y)/frameRate ;  // fraction of x length, length in seconds = frames / fps
	}
	for (ff =0; ff<f.length; ff++) {
		if (f[ff] > 1/maxFFTPeriod) {
			peaksOffset = ff;
			peaksOffsetY = y[ff];
			ff  =	f.length;
		}
		
	}
	 
	y = Array.resample(y, floor(f.length/ smoothFactor));
	y = Array.resample(y, f.length);

	// find minima in fft trace with lowest freq as the lowest freq for valid HR / override maxFFTperiod
	tolerance = FFTPeakToleranceFraction*y[0]/2.8;
	minLocs = Array.findMinima(y,tolerance);
	Array.getStatistics(minLocs, minFFTMinima, maxFFTMinima, meanFFTMinimam, sdFFTMinima);
	if (f[minFFTMinima] < peaksOffset) { 
		peaksOffset = minFFTMinima;
		peaksOffsetY = y[minFFTMinima];
	}
	xMinima = newArray(minLocs.length);
	yMinima = newArray(minLocs.length);
	for(m=0;m<minLocs.length;m++){
		xMinima[m] = f[minLocs[m]];
		yMinima[m] = y[minLocs[m]];
	}


	// find minima with lowest frequency as an estimate / check point for maxFFTPeriod or peakOffset

		
	//find maxima in FFT power spectra
	ySliced = Array.slice(y,peaksOffset,y.length);
	fSliced = Array.slice(f,peaksOffset,f.length);
	
	// peak-to-peak amplitude of sin wave is ~2.8*RMS
	//FFTPeakToleranceFraction = 0.05;   //moved to start for inclusion as function
	tolerance = FFTPeakToleranceFraction*y[0]/2.8;
	maxLocs = Array.findMaxima(ySliced,tolerance);
	

	if (maxLocs.length == 0) {
		freqMain = -1;
		perMain = -1;
		amplMain = -1;
		per2 = -1;
		ampl2 = -1;
		
	} else {

	//find the maxima in the Power Spectrum trace

		
		xMaxima = newArray(maxLocs.length);
		yMaxima = newArray(maxLocs.length);
		nMaxima = 0;
		
		for (jj= 0; jj < maxLocs.length; jj++){
				tempYmax = y[maxLocs[jj]+peaksOffset];
			
			if ( tempYmax*2.8 <= rawRange) {
				xMaxima[nMaxima]= f[maxLocs[jj]+peaksOffset];
				yMaxima[nMaxima] = tempYmax;
				nMaxima++;
			}
	  	}

		freqMain = xMaxima[0];
		perMain = 1/freqMain;
		amplMain = yMaxima[0];
		extraTxt = "";
		per2 = -1;
		ampl2 = -1;
			
		if (nMaxima>1) {
			per2 =  1/xMaxima[1];
			ampl2 = 2.8*yMaxima[1];
			extraTxt = " per2 " + 1/xMaxima[1] +  " Amp2 " + 2.8*yMaxima[1] + ""; 
		}
	}
		
	  if (doPlot == true) {
		Plot.create("Fourier amplitudes: "+windowType+" roi "+roi, "frequency (Hz)", "log(amplitude)", f, y);
		Plot.setLogScaleY(true);
		Array.getStatistics(y, min, max, mean, stdDev);
		yMin = pow(10,floor(log(min)/log(10)));
		yMax = pow(10,1+floor(log(max)/log(10)));
		Array.getStatistics(f, min, max, mean, stdDev);
		xMax = max;
		xMin = min;
        Plot.setLimits(xMin,xMax,yMin,yMax);
	  	Plot.setColor("red","red");
	  	Plot.setLineWidth(3);
		if (maxLocs.length != 0) Plot.add("circles",xMaxima,yMaxima);
	  	Plot.setColor("blue","blue");
	  	Plot.setLineWidth(3);
		if (minLocs.length != 0) Plot.add("circles",xMinima,yMinima);

		Plot.setColor("black");
		Plot.setLineWidth(1);
		Plot.setFontSize(16);
		if (maxLocs.length != 0) Plot.addText("Primary period "+ (1/freqMain) + " Amp " + amplMain*2.8 + " RMS " + rawRMS + " Range " + rawRange + extraTxt, 0.1, 1);
		if (maxLocs.length == 0) Plot.addText("no maxima found", 0.1, 1);
		Plot.show();
		}

	BPMpeaks =  d2s(60 / tPk2PkMean,1); 
	minutesAnalyzed = (1/60)*x.length/frameRate;
	//BPMfft =  d2s(frameRate * 60 / perMain ,1); 
	BPMfft =  d2s(60 * freqMain ,1); 
	label = imgAnalyzed+":"+roiName;
	nFrames = framesLast-framesFirst;
	BPMfft2 = d2s(60/per2,1);
	BPMcorrected = d2s(60/meantPk2PkNorm,2);
	
    returnArray = newArray(label, nFrames,""+framesFirst+"-"+framesLast, BPMpeaks, BPMfft, BPMfft2, nPeaks,tPk2PkMean, tPk2PkSD, BPMcorrected,meantPk2PkNorm, meantPk2PkHi,tPk2PkHi.length );
    returnArrayLabels = newArray("Label","Frames", "range", "BPM_peaks","BPM_fft","BPM_fft2","nPeaks","tPk2PkMean", "tPk2PkSD", "BPM_corrected","meantPk2PkCorrected", "meantPk2PkHi","n long periods");
    row = nResults;
    
	for (r=0; r<returnArray.length; r++) {
		setResult(returnArrayLabels[r], row, returnArray[r]);
	}
	print("\\Update1: Heart rate is " +BPMpeaks + " BPM (peak-to-peak time)");
	print("\\Update2: Heart rate is " +BPMfft + " BPM (fft)");
	//Array.print(xMaxima);
	
	
	return returnArray;

	
}

