/*
 * This script is part of the publication:  Martinek J. et al., "ARP2/3 complex associates with peroxisomes to participate in pexophagy in plants"
 * and is available under the CC-BY-4.0 License
 * 
 * This script loops through all the sample-containg subfolders and: 
 * 	1. finds image series named ch00 (corresponding to the green channel in out particular case, change accordingly if autofluorescence is detected in a different channel) obtained from the deconvolution step. 
 * 	2. Finds the corresponding .czi file in the parent folder containing the sample name in its filename and stores its parameters (scale and timeframe) 
 * 	3. Detects autofluorescence by automatic thresholdin, substracts thresholded regions from TrackMate label images and merges them into a multichannel image.
 * 	4. Runs ComDet plugin to detect colocalization events with specified parameters for original images and rotated control (in case it's required). 
 *
 * It requires the following file structure and minimum content:
 * 	> .../decon/.../Parent_folder
 * 		*sample_folder.czi
 * 		>Sample_folder
 * 			>ch00
 * 			>ch01
 * 			LblImg_ch00.tif
 * 			LblImg_ch01.tif
 * 			
 * 	Author: Judith García-González
 * 	Year: 2023
 * 	License: BSD-3
 * 	
 * 	Copyright (c) 2023, Judith García-González, Department of Experimental Plant Biology (Faculty of Science, Charles University in Prague)
 *	All rights reserved.

 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 			
 */


//Function that iterates through sample subfolders and process each of them
function processSequence (inputFolder) {
	//Get sample name
	File.openSequence(inputFolder);
	rename("ch00");
	//run(“Image Sequence…”, “open=&inputFolder sort”);
	current_img=getTitle();
	imgseqpath=getDirectory("image");
	parent_path=File.getParent(imgseqpath);
	print (imgseqpath);
	print (parent_path);
	coloc_path=parent_path+"\\coloc\\";
	File.makeDirectory(coloc_path);
	imgseqpath_split = split(imgseqpath, "\\");
	foldernum=imgseqpath_split.length;
	samplename= imgseqpath_split[imgseqpath_split.length-2];
	getPixelSize(unit, pixelWidth, pixelHeight);
	scale = 1/pixelWidth;
	print(scale);
	print(samplename);

	//Get folder with .czi files and find corresponding file to get time frame data
	path_split = split(imgseqpath, "(decon)");
	czi_folder = path_split[0];
	list = getFileList(czi_folder);
	f=0;
	for (f=0; f<list.length; f++) {
		if (endsWith(list[f], samplename+".czi")){
			print(czi_folder+list[f]);
			run("Bio-Formats Importer", "open=["+czi_folder+list[f]+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			frame_int=Stack.getFrameInterval();
			Stack.getUnits(X, Y, Z, Time, Value);
			close(list[f]);
		}
	}


	//save scale and time frame data into file
	param_f=File.open(parent_path+"\\coloc\\"+samplename+"_param.txt");
	print(param_f, "scale"+"\t"+"px/"+X+"\t"+scale+"\n"+"frameint"+"\t"+Time+"\t"+frame_int);
	File.close(param_f);


	//Determine autofluorescence, filter particles larger than 5um2
	selectWindow(current_img);
	run("Auto Threshold", "method=RenyiEntropy stack");
	/////////////run("Invert", "stack");
	run("Analyze Particles...", "size=5-Infinity show=Masks stack");
	run("Invert LUT");
	run("Invert", "stack");
	rename("Mask of MASK_ch00");
	
	//Substract autofluorescence mask from label images for channel "ch00" and "ch01". Save autofluorescence mask stack
	open(parent_path+"\\trackmate\\LblImg_ch00.tif");
	imageCalculator("AND create stack", "Mask of MASK_ch00","LblImg_ch00.tif");
	selectWindow("Result of Mask of MASK_ch00");
	rename("LblImg_ch00_AF");
	close("LblImg_ch00.tif");

	open(parent_path+"\\trackmate\\LblImg_ch01.tif");
	imageCalculator("AND create stack", "Mask of MASK_ch00","LblImg_ch01.tif");
	selectWindow("Result of Mask of MASK_ch00");
	rename("LblImg_ch01_AF");
	close("LblImg_ch01.tif");
	selectWindow("Mask of MASK_ch00");
	rename("ch00_AF");
	saveAs("tiff",parent_path+"\\trackmate\\ch00_AF.tif");


	//Create and save new hyperstack where red channel is ch01 and green channel is ch00 
	run("Merge Channels...", "c1=LblImg_ch01_AF c2=LblImg_ch00_AF create keep");
	//saveAs("tiff", parent_path+"\\coloc\\"+samplename+"_coloc_label.tif");
	saveAs("tiff", coloc_path+samplename+"_coloc_label.tif");
	
	
	run("Set Measurements...", "area mean min centroid center shape integrated stack display redirect=None decimal=3");

	//run colocalization analysis with ComDet using the merged two-channel hyperstack
	run("Detect Particles", "calculate max=15 plot rois=Ovals add=[All detections] summary=Reset ch1a=12 ch1s=1 ch2a=8 ch2s=1");
	selectWindow("Results");
	saveAs("results", coloc_path+samplename+"_coloc_results.csv"); 
	close("Results");
	selectWindow(samplename+"_coloc_label.tif");
	run("To ROI Manager");
	roiManager("Measure");
	saveAs("results", coloc_path+samplename+"_coloc_overlay.csv");
	roiManager("Deselect");
	roiManager("Save",coloc_path+samplename+"_coloc_overlay.zip");
	roiManager("Deselect");
	roiManager("Delete");
	close("Results");
	selectWindow("Summary");
	saveAs("results", coloc_path+samplename+"_coloc_summary.csv"); 
	close(samplename+"_coloc_summary.csv"); 
	
	//colocalization analysis control, rotated green channel (90 degrees right), red channel (90 degrees left)
	selectWindow("LblImg_ch00_AF");
	run("Rotate 90 Degrees Right");
	selectWindow("LblImg_ch01_AF");
	run("Rotate 90 Degrees Left");
	run("Merge Channels...", "c1=LblImg_ch01_AF c2=LblImg_ch00_AF create keep");
	saveAs("tiff", coloc_path+samplename+"_coloc_label_ctrl.tif");
	run("Detect Particles", "calculate max=15 plot rois=Ovals add=[All detections] summary=Reset ch1a=12 ch1s=1 ch2a=8 ch2s=1");
	selectWindow("Results");
	saveAs("results", coloc_path+samplename+"_coloc_results_ctrl.csv"); 
	close("Results");
	selectWindow(samplename+"_coloc_label_ctrl.tif");
	run("To ROI Manager");
	roiManager("Measure");
	saveAs("results", coloc_path+samplename+"_coloc_overlay_ctrl.csv");
	roiManager("Deselect");
	roiManager("Save",coloc_path+samplename+"_coloc_overlay_ctrl.zip");
	roiManager("Deselect");
	roiManager("Delete");
	
	close("Results");
	selectWindow("Summary");
	saveAs("results", coloc_path+samplename+"_coloc_summary_ctrl.csv");
	close(samplename+"_coloc_summary_ctrl.csv");
}

//User input of parent folder with sample subfolders
dir = getDirectory("Choose a Directory ");
//print(dir);
subdir_list = getFileList(dir);
s=0;
l=0;

//Iterate over sample subfolders and process ch00 image sequence
for (s = 0; s < subdir_list.length; s++) {
	//print (subdir_list[s]);
	if (endsWith(subdir_list[s], "/") && subdir_list[s] != "ch00/"){
		//print("entering");
		replace(subdir_list[s], "/", "\\");
		sub_lvl_list = getFileList(dir+subdir_list[s]);

		for (l = 0; l < sub_lvl_list.length; l++) {
			if (endsWith(sub_lvl_list[l], "/") && sub_lvl_list[l] == "ch00/"){
			inputFolder=dir+subdir_list[s]+"ch00\\";
			inputFolder=replace(inputFolder, "/", "\\");
			print(inputFolder);
			processSequence (inputFolder);
			//print("found it! sublevel!");
			close("*");
			
			}
		}
	}		
}




