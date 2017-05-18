var AallROIs;var Exptype;
macro "Sensitized emission correction Macro" {
	macroname="Sensitized emission correction macro";
	scwidth=screenWidth;
	scheight=screenHeight;
	plotwidth=round(450*1.5);
	plotheight=round(200*1.5);
	Acolor=newArray("black","red","blue","green","orange","magenta","pink","yellow","darkGray","gray","lightGray");
	Atrch_method=newArray("None","built-in 'Subtract background'","use outer rim of each cell");
	Asetmeasure=newArray("mean","modal","median","integrated");
	Ameasure=newArray("Mean","Mode","Median","RawIntDen");
	Atoanalch=newArray(0);
	Achnames=newArray(0);
	if(nImages==0)exit("No open images. Please open images or stacks to analyze and start macro again.");
	if(roiManager("count")==0)exit("No ROIs in ROI Manager. Please segment cells first and start macro again.");
	frames=0;
	for(i=0;i<nImages;i++){
		selectImage(i+1);
		Atoanalch=Array.concat(Atoanalch,getTitle());
		framesi=nSlices;
		if(framesi>frames)frames=framesi;
	};
	Dialog.create(macroname);
	Dialog.addMessage("Details to calculate sensitized emission channel (SE): \n  SE  = (S-"+fromCharCode(946,183)+"D-("+fromCharCode(947)+"-"+fromCharCode(945,946)+")"+fromCharCode(183)+"A)/(1-"+fromCharCode(946,948)+")\nD = donor excitation, donor emission channel\nS = donor excitation, acceptor emission channel\nA = acceptor excitation, acceptor detection channel\nCorrection factors from donor only images: \n  "+fromCharCode(946)+" = S/D\nCorrection factors from acceptor only images: \n  "+fromCharCode(945)+" = D/A; "+fromCharCode(947)+" = S/A; "+fromCharCode(948)+" D/S\n(setting "+fromCharCode(945)+" and "+fromCharCode(948)+" to 0 results in the simpler SE calculation with bleed-through and cross-excitation correction only)\n("+fromCharCode(945)+" and "+fromCharCode(948)+" coefficients correct for acceptor signal that may be detectable in the donor channel)");
	Dialog.addNumber(""+fromCharCode(945)+":",0);
	Dialog.addNumber(""+fromCharCode(946)+" (bleed-through coefficient):",0.245);
	Dialog.addNumber(""+fromCharCode(947)+" (cross-excitation coefficient):",0.072);
	Dialog.addNumber(""+fromCharCode(948)+":",0);
	Dialog.addChoice("Donor channel (D)",Atoanalch,Atoanalch[0]);
	Dialog.addChoice("Acceptor channel (A)",Atoanalch,Atoanalch[2]);
	Dialog.addChoice("Transfer/FRET channel (S)",Atoanalch,Atoanalch[1]);
	Dialog.addMessage(" ");
	Dialog.addCheckbox("Time-series data or Z-stack.",true);
	Dialog.show();
	List.set("SEalpha",Dialog.getNumber());
	List.set("SEbeta",Dialog.getNumber());
	List.set("SEgamma",Dialog.getNumber());
	List.set("SEdelta",Dialog.getNumber());
	List.set("SEdonor",Dialog.getChoice());
	List.set("SEacceptor",Dialog.getChoice());
	List.set("SEtransfer",Dialog.getChoice());
	diff=abs(parseFloat(List.get("SEdelta"))-(parseFloat(List.get("SEalpha")))/(parseFloat(List.get("SEgamma"))));
	if(!isNaN(diff)){
		if(parseFloat(List.get("SEdelta"))>=(parseFloat(List.get("SEalpha")))/(parseFloat(List.get("SEgamma")))){
			diff_p=diff/parseFloat(List.get("SEdelta"))*100;
		}else{
			diff_p=diff/(parseFloat(List.get("SEalpha")))/(parseFloat(List.get("SEgamma")))*100;
		};
		if(diff_p>1)waitForUser("Warning: equation "+fromCharCode(948)+"="+fromCharCode(945)+"/"+fromCharCode(947)+" is not true! Difference is "+diff_p+"%. Please check your coefficients again.");
	};
	Exptype=Dialog.getCheckbox();
	Achnames=Array.concat(Achnames,List.get("SEdonor"));
	Achnames=Array.concat(Achnames,List.get("SEacceptor"));
	Achnames=Array.concat(Achnames,List.get("SEtransfer"));
	rationame="SE corr "+List.get("SEtransfer")+"-"+List.get("SEdonor")+"-ratio";
	SEtransfer="SE corrected "+List.get("SEtransfer");
	Achnames=Array.concat(Achnames,rationame);	
	List.set("SEratiochannel",""+rationame);
	List.set("SEtransfer_corr",SEtransfer);
	Achnames=Array.concat(Achnames,SEtransfer);
	
	SE_channel=calculate_SE(List.get("SEtransfer_corr"),List.get("SEdonor"),List.get("SEacceptor"),List.get("SEtransfer"),parseFloat(List.get("SEalpha")),parseFloat(List.get("SEbeta")),parseFloat(List.get("SEgamma")),parseFloat(List.get("SEdelta")));
	imageCalculator("Divide create 32-bit stack", SE_channel,List.get("SEdonor"));
	run("Rename...", "title=["+rationame+"]");
	run("Colors...", "foreground=white background=black selection=white");
	
	run("Set Measurements...", "  mean redirect=None decimal=3");
	amountrois=roiManager("count");
	if(Exptype==0){
		frames=nSlices;
		run("Clear Results");
		amountrois=roiManager("count");
		icAresults=newArray(Achnames.length*amountrois);
		icASliceinfo=newArray(amountrois);
		icASD=newArray(Achnames.length*amountrois);
		icAresultsno=newArray(amountrois);
		if(amountrois>0){
			counter=0;
			for(channelno=0;channelno<Achnames.length;channelno++){
				ex=Achnames[channelno];	
				if(isOpen(ex)){
					setmeasure="mean";
					run("Set Measurements...", " mean standard stack redirect=["+ex+"] decimal=3");
					run("Select None");
					run("Clear Results");
					roiManager("Deselect");
					roiManager("Show All");
					roiManager("Measure");
					selectWindow("Results");
					for (ROIno =0; ROIno < nResults; ROIno++){//amountrois
						icAresults[ROIno+channelno*amountrois]=getResult("Mean", ROIno);
						icASD[ROIno+channelno*amountrois]=getResult("StdDev", ROIno);
						if(channelno==0){
							counter++;
							slice=getResult("Slice", ROIno)-1;
							icASliceinfo[ROIno+channelno*amountrois]=getResult("Slice", ROIno);
							icAresultsno[ROIno]=counter;
						};
					};
				};
			};
		};
		run("Set Measurements...", "  mean standard redirect=None decimal=3");
		run("Clear Results");
		cASliceinfo=icASliceinfo;
		cAresultsno=icAresultsno;
		mean=newArray(Achnames.length);
		SDmean=newArray(Achnames.length);
		SDEmean=newArray(Achnames.length);
		median=newArray(Achnames.length);
		SDmedian=newArray(Achnames.length);
		minimum=newArray(Achnames.length);
		maximum=newArray(Achnames.length);
		for(channelno=0;channelno<Achnames.length;channelno++){
			ex=Achnames[channelno];
			if(cAresultsno.length>0){
				cAresultsnoname=add_strings_to_array(cAresultsno,"ROI ","");
				PlotRArray(cAresultsnoname,cAresultsno,extract_array2(icAresults,channelno,amountrois),extract_array2(icASD,channelno,amountrois),"Barplot overview of single ROI measurements of "+ex,"ROI","Mean","Sensitized emission correction macro");
				//if(amountrois>1)PlotArray(cAresultsno,extract_array2(icAresults,channelno,amountrois),extract_array2(icASD,channelno,amountrois),"Lineplot overview of single ROI measurements of "+ex,"ROI",""+List.get("Measure"),List.get("origtitle"));
				print_in_results("ROI No",cAresultsno);
				print_in_results("ROI from slice",cASliceinfo);
				x=extract_array2(icAresults,channelno,amountrois);
				print_in_results("Mean of "+ex,extract_array2(icAresults,channelno,amountrois));
				print_in_results("StdDev of "+ex,extract_array2(icASD,channelno,amountrois));
			};
		};
		for(channelno=0;channelno<Achnames.length;channelno++){
			ex=Achnames[channelno];
			if(cAresultsno.length>0){
				mean[channelno]=calMean(extract_array2(icAresults,channelno,amountrois));
				SDmean[channelno]=calSD(extract_array2(icAresults,channelno,amountrois));
				SDEmean[channelno]=calSDE(extract_array2(icAresults,channelno,amountrois));
			};
			if(cAresultsno.length>3){
				median[channelno]=calMedian(extract_array2(icAresults,channelno,amountrois));
				SDmedian[channelno]=calQuartilsdiff(extract_array2(icAresults,channelno,amountrois));
				minimum[channelno]=calMin(extract_array2(icAresults,channelno,amountrois));
				maximum[channelno]=calMax(extract_array2(icAresults,channelno,amountrois));
			};
			if(cAresultsno.length<=3){
				median[channelno]=0;SDmedian[channelno]=0;minimum[channelno]=0;maximum[channelno]=0;		
			};
			print("Results of "+ex);
			print("Total mean = "+mean[channelno]+" "+fromCharCode(177)+" "+SDmean[channelno]);
			print("Total median = "+median[channelno]+" "+fromCharCode(177)+" "+SDmedian[channelno]);
			print("Total minimum = "+minimum[channelno]);
			print("Total maximum = "+maximum[channelno]);
			print("Total StdErr = "+SDEmean[channelno]);
			print(" ");
		};
		if(cAresultsno.length>0){
			Atitle=newArray(Achnames.length);
			xValues=newArray(Achnames.length);
			for(channelno=0;channelno<Achnames.length;channelno++){
				ex=Achnames[channelno];
				chno=channelno+1;
				xValues[channelno]=chno;
				Atitle[channelno]=ex;
			};
			PlotRArray(Atitle,xValues,mean,SDmean,"Overview Mean of all ROIs","Channel","Mean intensity",List.get("origtitle"));
		};
	};
	if(Exptype==1){	
		AallROIs=newArray(frames*amountrois*Achnames.length);
		for(ch=0;ch<Achnames.length;ch++){
			Atime=getFrameArray(Achnames[ch]);
			multimeasureresultsplot(Achnames[ch],"Plot - "+Achnames[ch]+" Mean intensity of all ROIs versus time",""+Achnames[ch]+" Mean intensity",Atime,ch);	
		};
		run("Clear Results");
		for(ch=0;ch<Achnames.length;ch++){
			for(ROIno=0;ROIno<amountrois;ROIno++){
				c=ROIno+1;
				ROI="ROI No "+c;
				print_in_results(""+Achnames[ch]+" "+ROI,extract_array(AallROIs,ROIno,ch,amountrois,frames));
			};
		};
	};
	beep();
	waitForUser("Analysis complete!");
};
function add_strings_to_array(array,before,after){
	array2=newArray(array.length);
	for(i=0;i<array.length;i++){
		array2[i]=""+before+array[i]+after;		
	};
	return array2;
};
function PlotRArray(Atitle,xValues,yValues,ASD,plottitle,xaxis,yaxis,origtitles){
	if(xValues.length>1&&xValues.length<200){	
		Alimits=removeNaN(xValues);
		Array.getStatistics(Alimits, xMin, xMax, mean, stdDev);
		Alimits=Array.concat(yValues,ASD);
		Alimits=removeNaN(Alimits);
		Alimits=Array.concat(Alimits,0);
		Array.getStatistics(Alimits, yMin, yMax, mean, stdDev);
		Alimits=removeNaN(ASD);
		Array.getStatistics(Alimits, min, Errorbars, mean, stdDev);
		xMaxorig=xMax;
		xMinorig=xMin;
		xspace=abs(xMin-xMax)*0.05;
		if(Errorbars>0)yspace=abs(Errorbars*1.05);
		if(Errorbars==0)yspace=abs(yMin-yMax)*0.05;
		xMin = xMin-xspace-0.4;xMax=xMax+xspace+0.4;yMin=yMin-yspace;yMax=yMax+yspace;
		plotname=plottitle;
		if(plotheight/Atitle.length<14){
			plotheight=Atitle.length*14.2;
		};
		stringlength1=0;
		for(f=0;f<Atitle.length;f++){
			if(lengthOf(Atitle[f])>stringlength1)stringlength1=lengthOf(Atitle[f]);		
		};
		for(l=0;l<2;l++){	
			heightofchar=14/plotheight;
			widthofchar=7/plotwidth;
			begincharheight=1.2*heightofchar;
			stringlength=stringlength1*widthofchar;
			textlength=stringlength+10*widthofchar;
			textwidth=1-textlength;
			timeswider=textlength+1+widthofchar;
			if(l==0)xMax=xMax*timeswider;
			linewidth=0.05*xMaxorig;
			if(l==0){
				plotwidthnew=timeswider*plotwidth;
				plotwidth=plotwidthnew;
			};
		};
		Plot.create(plotname, xaxis, yaxis);
		Plot.setFrameSize(plotwidth, plotheight);
		Plot.setLimits(xMin, xMax, yMin, yMax);
		Plot.setColor("black");
		Plot.drawLine(xMinorig-0.45, 0, xMaxorig+0.45, 0);
		for(i=0;i<xValues.length;i++){
			if(i==0)barwidth=abs((xValues[i+1]-xValues[i]))*0.75;
			if(i>0&&i<xValues.length-2){
				barwidth1=abs(xValues[i+1]-xValues[i]);
				barwidth2=abs(xValues[i]-xValues[i-1]);
				if(barwidth1>=barwidth2)barwidth=barwidth2*0.75;
				if(barwidth2>barwidth1)barwidth=barwidth1*0.75;	
			};
			if(i==xValues.length-1)barwidth=abs((xValues[i]-xValues[i-1]))*0.75;
			drawbar(xValues[i],yValues[i],ASD[i],barwidth);	
		};
		Plot.setColor("black");
		setJustification("right");
		for(i=1;i<=Atitle.length;i++){
			ii=i-1;
			Plot.addText("No: "+i+" - "+Atitle[ii], 1-0.5*widthofchar, begincharheight+(ii*heightofchar));
		};
		setJustification("center");
		Plot.addText(plotname, 0.5, 0);
		Plot.setColor("black");
		Plot.setLineWidth(1);
		setJustification("right");
		Plot.addText(origtitles, 1, 1+2.5*heightofchar);
		setJustification("left");
		Plot.setColor("gray");
		Plot.show();
		//run("Profile Plot Options...", "width="+plotwidth+" height="+plotheight+" minimum=0 maximum=0 interpolate draw");
		function drawbar(x,y,yerr,barwidth){
			Plot.setColor("lightgray");
			Plot.setLineWidth(1);
			//barwidth=0.75;
			rep=500;
			for(i=0;i<rep;i++){
				add=(barwidth/rep)*i;
				Plot.drawLine(x-barwidth/2+add, 0, x-barwidth/2+add, y);	
			};
			Plot.setColor("darkgray");
			Plot.drawLine(x-barwidth/2, 0, x-barwidth/2, y);
			Plot.drawLine(x+barwidth/2, 0, x+barwidth/2, y);
			Plot.drawLine(x-barwidth/2, 0, x+barwidth/2, 0);
			Plot.drawLine(x-barwidth/2, y, x+barwidth/2, y);
			Plot.setColor("black");
			Plot.setLineWidth(1);
			Plot.drawLine(x, y-yerr/2, x, y+yerr/2);
			Plot.drawLine(x-barwidth/3, y-yerr/2, x+barwidth/3, y-yerr/2);
			Plot.drawLine(x-barwidth/3, y+yerr/2, x+barwidth/3, y+yerr/2);
		};
	};
};
function extract_array2(array,channelno,rows){//[x+y*xmax]
	if(array.length<channelno*rows)exit("Mistake occured. Array is not long enough for multiple dimensions!");
	sarray=newArray(rows);
	for(i=0;i<rows;i++){
		sarray[i]=array[i+channelno*rows];	
	};
	return sarray;
};
function extract_array(array,ROIno,channelno,amountrois,frames){//[x+y*xmax+z*xmax*ymax]
	if(array.length<channelno*amountrois*frames)exit("Mistake occured. Array is not long enough for multiple dimensions!");
	sarray=newArray(frames);
	for(slice=0;slice<frames;slice++){
		sarray[slice]=array[slice+ROIno*frames+channelno*amountrois*frames];			
	};
	return sarray;
};
function getfromResults(columname,rows){
	Avalues=newArray(rows);
	Array.fill(Avalues,NaN);
	if(rows<=nResults){
		for (y =0; y < rows; y++){
			Avalues[y] = getResult(columname, y);
		};
	};
	return Avalues;
};
function print_in_results(columnname,array){
	if(calMean(array)!=0){
		rows=array.length;
		for(i=0;i<rows;i++){
			setResult (columnname,i,array[i]);
		};
		updateResults();
	};		
};
function multimeasureresultsplot(channel,plotname,yaxis,Atime,channelno){//AallROIs
	if(isOpen(channel)){
		run("Clear Results");
		selectWindow(channel);
		wait(100);
		frames=nSlices;
		run("Set Measurements...", "mean redirect=["+channel+"] decimal=3");
		roiManager("Deselect");
		roiManager("Multi Measure");
		amountrois=roiManager("count");
		ROItraces=newArray(frames*amountrois);
		for(c = 1; c <= amountrois; c++){
			ROIno=c-1;
			value="Mean"+c;
			AROI=getfromResults(value,frames);
			if(channelno<=Achnames.length){
				for(slice=0;slice<frames;slice++){
					AallROIs[slice+ROIno*frames+channelno*amountrois*frames]=AROI[slice];	
				};	
			};
			for(slice=0;slice<frames;slice++){
				ROItraces[slice+ROIno*frames]=AROI[slice];	
			};
		};
		ROInames=newArray(amountrois);
		AanalROIs=newArray(amountrois);
		for(i=0;i<amountrois;i++){
			ROI=i+1;
			ROInames[i]="ROI "+ROI;	
			AanalROIs[i]=i;
		};
		if(Atime.length>1)PlotmultipleArrays(Atime,ROItraces,ROInames,AanalROIs,plotname,"Frame",yaxis);			
		run("Clear Results");
	};
};
function removeNaN(aA){
	c=0;
	while(c<aA.length){
		if(isNaN(aA[c])){
			bA=Array.slice(aA,0,c);
			cA=Array.slice(aA,c+1,aA.length);
			aA=Array.concat(bA,cA);			
		}else c++;
	};
	return aA;
};
function calMax(array){
	array=removeNaN(array);
	Array.getStatistics(array,min,max,mean,stdDev);
	return max;	
};
function PlotmultipleArrays(xValues,yValues,Awindownames,Atoanal,plottitle,xaxis,yaxis){
	Alimits=removeNaN(xValues);
	Array.getStatistics(Alimits, xMin, xMax, mean, stdDev);
	Alimits=removeNaN(yValues);
	Array.getStatistics(Alimits, yMin, yMax, mean, stdDev);
	xMaxorig=xMax;
	xMinorig=xMin;
	xspace=abs(xMin-xMax)*0.05;
	yspace=abs(yMin-yMax)*0.05;
	xMin = xMin-xspace;xMax=xMax+xspace;yMin=yMin-yspace;yMax=yMax+yspace;
	if(plotheight/Atoanal.length<14){
		plotheight=Atoanal.length*14.2;
	};
	stringlengthes=newArray(Atoanal.length);
	for(l=0;l<2;l++){	
		heightofchar=14/plotheight;//*********************************
		widthofchar=7/plotwidth;//*********************************
		begincharheight=1.2*heightofchar;//*********************************
		for(i=0;i<Atoanal.length;i++){
			c=Atoanal[i];
			stringlengthes[i]=lengthOf(Awindownames[c])*widthofchar;		
		};
		stringlength=calMax(stringlengthes);
		textlength=stringlength+4*widthofchar;
		textwidth=1-textlength;//*********************************
		timeswider=textlength+1+widthofchar;//*********************************
		if(l==0)xMax=xMax*timeswider;//*********************************
		linewidth=0.05*xMaxorig;//*********************************
		if(l==0){
			plotwidthnew=timeswider*plotwidth;
			plotwidth=plotwidthnew;
			
		};
	};
	Plot.create(plottitle, xaxis, yaxis);
	Plot.setFrameSize(plotwidth, plotheight);
	Plot.setLineWidth(1);
	Plot.setLimits(xMin, xMax, yMin, yMax);
	col=0;
	for(line=0;line<Atoanal.length;line++){
		channelno=Atoanal[line];
		Aline=extract_array2(yValues,line,frames);
		add_line(xValues,Aline,Acolor[col],Awindownames[channelno],line);
		col++;
		if(col==Acolor.length)col=0;
	};
	//Captions
	setJustification("center");
	Plot.addText(plottitle, 0.5, 0);
	Plot.setLineWidth(1);
	Plot.show();
};
function add_line(xvalues,Aline,color,name,line){
	if(calMean(Aline)!=0){	
		xvalues=Array.trim(xvalues, Aline.length);
		Plot.setColor(color);
		Plot.add("line",xvalues,Aline);
		setJustification("right");
		Plot.addText(""+fromCharCode(9472,9472),1-stringlength-widthofchar,begincharheight+(line*heightofchar));
		Plot.setColor("black");
		Plot.addText(name, 1-0.5*widthofchar, begincharheight+(line*heightofchar));	
	};
};	
function getFrameArray(channel){
	selectWindow(channel);
	wait(100);
	frames=nSlices;
	array=newArray(nSlices);
	for(i=0;i<nSlices;i++){
		array[i]=i+1;
	};
	return array;
};
function get_theor_maxvalue(){
	Bitd=bitDepth();
	maxvalue=pow(2, Bitd)-1;
	return maxvalue;	
};
function getmaximumpixel(channel){
	selectWindow(channel);
	wait(100);
	if(nSlices==1)getStatistics(area, mean, min, max, std, histogram);
	if(nSlices>1)Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	return max;
};
function resize(){
	getLocationAndSize(x, y, width, height);
	wratio=height/width;
	hratio=width/height;
	wfactor=(scwidth-220)/scwidth;
	hfactor=(scheight-220)/scheight;
	wspace=screenWidth*wfactor;
	hspace=screenHeight*hfactor;
	if(parseFloat(List.get("chnumber"))>2){
		setLocation(x,y,wspace/2,hspace/2);
	};
	if(parseFloat(List.get("chnumber"))<=2){
		setLocation(x,y,wspace/2,hspace);
	};
};
function calMean(arrayf){ //Caclulates the Mean value of arrayf
	arrayf=removeNaN(arrayf);
	Array.getStatistics(arrayf,min,max,mean,stdDev);
	return mean;
};
function removeNaN(aA){
	c=0;
	aA=Array.concat(aA);
	while(c<aA.length){
		if(isNaN(aA[c])){
			bA=Array.slice(aA,0,c);
			cA=Array.slice(aA,c+1,aA.length);
			aA=Array.concat(bA,cA);			
		}else c++;
	};
	return aA;
};
function arrange_and_wait(threshold,channel,channel2,message,updateRM,tool){
	selectWindow(channel);
	wait(100);
	getLocationAndSize(x, y, width, height);
	resetMinAndMax();
	run("Enhance Contrast", "saturated=0.35");
	xorig=x;
	yorig=y;
	worig=width;
	horig=height;
	if(isOpen(channel2)){	
		selectWindow(channel2);
		wait(100);
		getLocationAndSize(x, y, width, height);
		xorig2=x;
		yorig2=y;
		worig2=width;
		horig2=height;
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");
	};
	wfactor=(screenWidth-220)/screenWidth;
	hfactor=(screenHeight-220)/screenHeight;
	wspace=screenWidth*wfactor;
	hspace=screenHeight*hfactor;
	if(isOpen("Synchronize Windows")){
		selectWindow("Synchronize Windows");
		setLocation(wspace/2,hspace);		
	};
	if(isOpen(channel2)){
		selectWindow(channel2);
		setLocation(wspace/2,0);
		roiManager("Show All with labels");
	};
	selectWindow(channel);
	wait(100);
	setLocation(0,0,wspace/2,hspace);	
	getLocationAndSize(x, y, width, height);
	selectWindow(channel);
	wait(100);
	roiManager("Show All with labels");
	if(threshold==1){
		run("Threshold...");
		setAutoThreshold("Triangle dark stack");
		getThreshold(lower, upper);
		setThreshold(lower, upper);
		setThreshold(lower, upper);
	};
	set_Tool(tool);
	waitForUser(message);
	if(isOpen(channel2)){	
		selectWindow(channel2);
		wait(100);
		setLocation(xorig2,yorig2,worig2,horig2);	
		if(updateRM)roiManager("Show All with labels");
	};
	selectWindow(channel);
	wait(100);
	setLocation(xorig,yorig,worig,horig);
	if(updateRM)roiManager("Show All with labels");
};
function set_Tool(tool){
	initial=IJ.getToolName();
	found=0;
	for(i=0;i<=22;i++){
		setTool(i);
		x=IJ.getToolName();
		if(x==tool){
			i=22;
			found=1;
		};
	};
	if(found==0){
		setTool(initial);
	};
};
function set_Tool(tool){
	initial=IJ.getToolName();
	found=0;
	for(i=0;i<=22;i++){
		setTool(i);
		x=IJ.getToolName();
		if(x==tool){
			i=22;
			found=1;
		};
	};
	if(found==0){
		setTool(initial);
	};
};
function calMedian(array){ //Calculates the Median of Amedianrun
	counter=0;
	for(i=0;i<array.length;i++){
		if(!isNaN(array[i]))counter++;
	};
	if(counter>3){
		medianrun=newArray(counter);
		counter=0;
		for(i=0;i<array.length;i++){
			if(!isNaN(array[i])){
				medianrun[counter]=array[i];
				counter++;	
			};
		};
		Array.sort(medianrun);
		l=lengthOf(medianrun)+1;
		Median=l/2;
		Medianf=floor(Median);
		Mediand=Median-Medianf;
		if(Mediand!=0){
			median=(medianrun[Medianf-1]+medianrun[Medianf])/2; 
			};
		if(Mediand==0){
			median=medianrun[Medianf-1];
		};
		return median;
	};
	if(counter<=3){
		median=0;
		return median;	
	};
};
function calQuartilsdiff(array){ //Calculates the difference between the first and third Quartile
	counter=0;
	for(i=0;i<array.length;i++){
		if(!isNaN(array[i]))counter++;
	};
	if(counter>3){
		Amedianrun=newArray(counter);
		counter=0;
		for(i=0;i<array.length;i++){
			if(!isNaN(array[i])){
				Amedianrun[counter]=array[i];
				counter++;	
			};
		};
		Array.sort(Amedianrun);
		l=lengthOf(Amedianrun)+1;
		Q1=l/4;
		Q1f=floor(Q1);
		Q1d=Q1-Q1f;
		Q3=l*3/4;
		Q3f=floor(Q3);
		Q3d=Q3-Q3f;
		if(Q1d==0){
			qdiff=abs((Amedianrun[Q1f-1]-Amedianrun[Q3f-1])/2);	
		};
		if(Q1d!=0){
			Q1v=abs(Amedianrun[Q1f-1]+((Amedianrun[Q1f]-Amedianrun[Q1f-1])*Q1d));			
			Q3v=abs(Amedianrun[Q3f-1]+((Amedianrun[Q3f]-Amedianrun[Q3f-1])*Q3d));	
			qdiff=abs((Q3v-Q1v)/2);
		};
		return qdiff;
	};
	if(counter<=3){
		qdiff=0;
		return qdiff;	
	};
};
function calMean(arrayf){ //Caclulates the Mean value of arrayf
	arrayf=removeNaN(arrayf);
	Array.getStatistics(arrayf,min,max,mean,stdDev);
	return mean;
};
function calSD(arrayf){ //Calculates the Standard Deviation of arrayf
	arrayf=removeNaN(arrayf);
	Array.getStatistics(arrayf,min,max,mean,stdDev);
	if(arrayf.length<=2)stdDev=0;
	return stdDev;
};
function calSDE(array){ //Calculates the Standard Deviation of arrayf
	arrayf=removeNaN(array);
	Array.getStatistics(arrayf,min,max,mean,stdDev);
	if(arrayf.length<=2)return 0;
	sderr=stdDev/sqrt(arrayf.length);
	return sderr;
};
function calMin(array){
	array=removeNaN(array);
	Array.getStatistics(array,min,max,mean,stdDev);
	return min;
};
function calMax(array){
	array=removeNaN(array);
	Array.getStatistics(array,min,max,mean,stdDev);
	return max;	
};
function calculate_SE(transfer_seCorr,donor,acceptor,transfer,alpha,beta,gamma,delta){//Rbleedthrough,Rcrossexitation
	batchmode=is("Batch Mode");
	if(!batchmode){
		setBatchMode(true);
	};
	alpha=parseFloat(alpha);
	beta=parseFloat(beta);
	gamma=parseFloat(gamma);
	delta=parseFloat(delta);
	donor_corr="Corrected "+donor;
	acceptor_corr="Corrected "+acceptor;
	selectWindow(donor);
	wait(100);
	run("Duplicate...", "title=["+donor_corr+"] duplicate range=1-["+nSlices+"]");
	selectWindow(acceptor);
	wait(100);
	run("Duplicate...", "title=["+acceptor_corr+"] duplicate range=1-["+nSlices+"]");
	selectWindow(transfer);
	wait(100);
	run("Duplicate...", "title=["+transfer_seCorr+"] duplicate range=1-["+nSlices+"]");
	selectWindow(donor_corr);
	wait(100);
	run("Multiply...", "value=["+beta+"] stack");
	selectWindow(acceptor_corr);
	wait(100);
	acceptor_factor=gamma-alpha*beta;
	run("Multiply...", "value=["+acceptor_factor+"] stack");
	imageCalculator("Subtract stack", transfer_seCorr,donor_corr);
	imageCalculator("Subtract stack", transfer_seCorr,acceptor_corr);
	selectWindow(transfer_seCorr);
	wait(100);
	corr=1-beta*delta;
	run("Divide...", "value=["+corr+"] stack");
	if(isOpen(donor_corr)){
		selectWindow(donor_corr);
		wait(100);
		close();
	};
	if(isOpen(acceptor_corr)){
		selectWindow(acceptor_corr);
		close();
	};
	setBatchMode("exit and display");
	return transfer_seCorr;
};