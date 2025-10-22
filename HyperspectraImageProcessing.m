(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["HyperspectralImageProcessing`"];


(* ::Input::Initialization:: *)
readInHeadwallFileImages::"usage"="This function takes two arguments - the first is the path to the header (extension .hdr) of the hyperspectral image of the leaf, the second is the path to the header (.hdr) of the image used for dark calibration (an image recorded with the camera objective blocked by opaque material). The files with data saved in a binary form corresponding to the header, should be located in the same directory. If the second argument is not specified, the dark calibration is not performed.";
getImageData::"usage"="";
getImageForWavelengthWithoutInterpolation::"usage"="";
getRandomPixelBlock::"usage"="Arguments: (1) region ((usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) within which we are seekig a random pixel block,(2) length of the sought block in pixels. Returns a rectangle, representing a random block of pixels within the specified region.";


(*keys (method names)*)
getSpectralBandCount::usage="Method of the hyperspectral image object. Arguments: none. Returns: the number of spectral bands."; 
getLineCount::usage="Method of the hyperspectral image object. Arguments: none. Returns: the number of lines in each image.";
getWavelengthsInNanometers::usage="Method of the hyperspectral image object. Arguments: none. Returns: the list of wavelengths of the hyperspectral image in nanometers.";
getImageDimensions::"usage"="Method of the hyperspectral image object. Returns the list with two elements: number of ";
calculateChloroplastMovementIndex::"usage"="Method of the hyperspectral image object. Arguments: the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area for which chloroplast movement index is sought. Returns: avarage chloroplast movement indexs calculated for area within the specified region.";
displayChloroplastMovementIndexHistogram::"usage"="Method of the hyperspectral image object. ";
displayChloroplastMovementIndex::"usage"="Method of the hyperspectral image object. Arguments: none. Returns: an image showing chloroplast movement index calculated for each pixel of the hyperspectral image. The image contains a legend showing the range color-coded values of chloroplast movement index.";
displayChloroplastMovementIndexImageCropped::"usage"="Method of the hyperspectral image object. Arguments: (1) the first, (2) the second row of the image to keep while cropping, (3) the first, (4) the second column of the image to keep while cropping. Arguments: Returns: a cropped image showing chloroplast movement index calculated for each pixel of the hyperspectral image. The image contains a legend showing the range color-coded values of chloroplast movement index.";
displayChloroplastMovementIndexImageCroppedWithoutLegend::"usage"="Method of the hyperspectral image object. Method of the hyperspectral image object. Arguments: (1) the first, (2) the second row of the image to keep while cropping, (3) the first, (4) the second column of the image to keep while cropping. Arguments: Returns: a cropped image showing chloroplast movement index calculated for each pixel of the hyperspectral image. No legend is displayed.";
displayChloroplastMovementIndexImageCroppedWithoutLegendAndBackground::"usage"="Method of the hyperspectral image object. ";
getVISReflectanceImage::"usage"="Method of the hyperspectral image object. Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area occupied by the 100% reflectance standard. Returns: a matrix of reflectance values in the visible (400 - 700 nm) range calculated for each pixel of the hyperspectral image.";
displayVISReflectanceImage::"usage"="Method of the hyperspectral image object. Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area occupied by the 100% reflectance standard. Returns: an image showing reflectance in the visible (400 - 700 nm) range calculated for each pixel of the hyperspectral image.";
displayVISReflectanceImageCropped::"usage"="Method of the hyperspectral image object. Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area occupied by the 100% reflectance standard, (2) the first, (3) the second row of the image to keep while cropping, (4) the first, (5) the second column of the image to keep while cropping. Returns: a cropped image showing reflectance in the visible (400 - 700 nm) range calculated for each pixel of the hyperspectral image.";
displayReflectanceImageCropped::"usage"="Method of the hyperspectral image object. Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area occupied by the 100% reflectance standard, (2) the lowest wavelength in nm to consider during calculation of reflectance, (3) the greatest wavelength, (4) the first, (5) the second row of the image to keep while cropping, (6) the first, (7) the second column of the image to keep while cropping. Returns: a cropped image showing reflectance in the arbitrary range calculated for each pixel of the hyperspectral image.";
displayNDVImageCropped::"usage"="Method of the hyperspectral image object.  Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area occupied by the 100% reflectance standard, (2) the first, (3) the second row of the image to keep while cropping, (4) the first, (5) the second column of the image to keep while cropping. Returns: a cropped image the showing normalized difference vegetation indedx range calculated for each pixel of the hyperspectral image.";
getImageDataByIndex::"usage"="Method of the hyperspectral image object. Arguments: (1) the index of the image (along the spectral direction). Returns: a metrix of pixel values for a single wavlenegth value, corresponding to the specified index.";
displayImageByIndex::"usage"="Method of the hyperspectral image object. Arguments: (1) the index of the image (along the spectral direction). Returns: an image for a single wavelength value, corresponding to the specified index.";
displayImageByIndexCustomNormalized::"usage"="Method of the hyperspectral image object. ";
getImageDataByWavelength::"usage"="Method of the hyperspectral image object. Arguments: (1) wavelength in nanometers. Returns a matrix of image values corresponding to the specified wavelength.";
displayImageDataByWavelength::"usage"="Method of the hyperspectral image object. Arguments: (1) wavelength in nanometers. Returns an image corresponding to the specified wavelength";
displayImageDataByWavelengthConstrastAdjusted::"usage"="Method of the hyperspectral image object. ";
extractRawCountsInRegionForEachWavelength::"usage"="Method of the hyperspectral image object. ";
extractMeanSpectrumWithinWholeImage::"usage"="Method of the hyperspectral image object. ";
extractSpectraWithinWholeImage::"usage"="Method of the hyperspectral image object. ";
extracRawSpectrumWithoutDarkCalibration::"usage"="Method of the hyperspectral image object. ";
extracRawSpectrum::"usage"="Method of the hyperspectral image object. Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area for which spectrum is required. Returns a list of mean image values for each wavelength.";
extractReflectanceSpectrumWithoutDarkCalibration::"usage"="Method of the hyperspectral image object. ";
extractReflectanceSpectrum::"usage"="Method of the hyperspectral image object. Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area for which reflectance is required, (2) the region occupied by the 100% reflectance standard. Returns a list of mean image values for each wavelength, normalized by the values calculated for the reference region.";
extractReflectanceWithinRegionUsingWhiteReferenceCurve::"usage"="Method of the hyperspectral image object. Arguments: (1) the region (usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area for which reflectance spectrum is required, (2) the curve showing absolute spectrum of light reflected from 100% reflectance standard. Returns a list of mean image values for each wavelength, normalized by the values calculated for the reference region.";
extractDifferenceInReflectanceBetweenRegionsUsingReferenceWithoutDarkCalibration::"usage"="Method of the hyperspectral image object. ";
extractDifferenceInReflectenceSpectrum::"usage"="Method of the hyperspectral image object. Arguments: (1) the first region (region A; usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area to use in calculation of reflectance, (2) the second region (region B), (3) he region occupied by the 100% reflectance standard. Returns: a curve showing difference in reflectance between region B and A (B minus A).";
extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurve::"usage"="Method of the hyperspectral image object. Arguments: (1) the first region (region A; usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) of the image area to use in calculation of reflectance, (2) the second region (region B), (3) the curve showing absolute spectrum of light reflected from 100% reflectance standard. Returns: a curve showing difference in reflectance between region B and A (B minus A).";
extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurveWithCalibration::"usage"="Method of the hyperspectral image object. ";
extractDifferenceInMeanSpectrumBetweenRegionsWithRespectToTwoReferences::"usage"="Method of the hyperspectral image object. ";
extractDifferenceInReflectanceBetweenRegionsUsingTwoWhiteReferenceRegions::"usage"="Method of the hyperspectral image object. ";
classifyChloroplastResponse::"usage"=" Arguments: (1) the region occupied by dark-adapted leaf part (region A; usually a Disk or Polygon, can be other object for which the Mathematica built-in command RegionQ returns True) within the image area, (2) the region occupied by the irradiated leaf part (region B), (3) he region occupied by the 100% reflectance standard. Returns: 1 if the irradiated curve shows ";


(* ::Input::Initialization:: *)
Begin["`Private`"];


(* ::Input::Initialization:: *)
(*keys (method names)*)
getSpectralBandCount="getSpectralBandCount";
getLineCount="getLineCount";
getWavelengthsInNanometers="getWavelengthsInNanometers";
getImageDataByIndex="getImageDataByIndex";
getImageDimensions="getImageDimensions";
calculateChloroplastMovementIndex="calculateChloroplastMovementIndex";
displayChloroplastMovementIndexHistogram="displayChloroplastMovementIndexHistogram";
displayChloroplastMovementIndex="displayChloroplastMovementIndex";
displayChloroplastMovementIndexImageCropped="displayChloroplastMovementIndexImageCropped";displayChloroplastMovementIndexImageCroppedWithoutLegend="displayChloroplastMovementIndexImageCroppedWithoutLegend";displayChloroplastMovementIndexImageCroppedWithoutLegendAndBackground="displayChloroplastMovementIndexImageCroppedWithoutLegendAndBackground";
getVISReflectanceImage="getReflectanceImage";
displayVISReflectanceImage="displayVISReflectanceImage";
displayVISReflectanceImageCropped="displayVISReflectanceImageCropped";displayReflectanceImageCropped="displayReflectanceImageCropped";displayNDVImageCropped="displayNDVImageCropped";displayImageByIndex="displayImageByIndex";displayImageByIndexCustomNormalized="displayImageByIndexCustomNormalized";getImageDataByWavelength="getImageDataByWavelength";displayImageDataByWavelength="displayImageDataByWavelength";displayImageDataByWavelengthConstrastAdjusted="displayImageDataByWavelengthConstrastAdjusted";extractRawCountsInRegionForEachWavelength="extractRawCountsInRegionForEachWavelength";extractMeanSpectrumWithinWholeImage="extractMeanSpectrumWithinWholeImage";extractSpectraWithinWholeImage="extractSpectraWithinWholeImage";extracRawSpectrumWithoutDarkCalibration="extracRawSpectrumWithoutDarkCalibration";extracRawSpectrum="extracRawSpectrumWithDarkCalibration";extractReflectanceSpectrumWithoutDarkCalibration="extractReflectanceSpectrumWithoutDarkCalibration";extractReflectanceSpectrum="extractReflectanceSpectrum";extractReflectanceWithinRegionUsingWhiteReferenceCurve="extractReflectanceWithinRegionUsingWhiteReferenceCurve";extractDifferenceInReflectanceBetweenRegionsUsingReferenceWithoutDarkCalibration="extractDifferenceInReflectanceBetweenRegionsUsingReferenceWithoutDarkCalibration";extractDifferenceInReflectenceSpectrum="extractDifferenceInReflectenceSpectrum";extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurve="extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurve";extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurveWithCalibration="extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurveWithCalibration";extractDifferenceInMeanSpectrumBetweenRegionsWithRespectToTwoReferences="extractDifferenceInMeanSpectrumBetweenRegionsWithRespectToTwoReferences";extractDifferenceInReflectanceBetweenRegionsUsingTwoWhiteReferenceRegions="extractDifferenceInReflectanceBetweenRegionsUsingTwoWhiteReferenceRegions";
classifyChloroplastResponse="classifyChloroplastResponse";


(* ::Input::Initialization:: *)
(*It is generally necessary to load any package also on the master kernel*)Quiet[Needs["Combinatorica`"]];


(* ::Input::Initialization:: *)
Quiet[ParallelNeeds["Combinatorica`"]];


(* ::Input::Initialization:: *)
DistributeDefinitions[Combinatorica`RandomPermutation];


convolutionalClassifier=Import[NotebookDirectory[]<>"preTrainedClassifierReflectanceSpectra.wmlf"];


(* ::Input::Initialization:: *)
getGitelsonAbsorptionCofficient[spectrum_]:=With[{},With[{nirReflectance=getMeanReflectanceWithinRange[spectrum,770,800]},Table[{spectrum[[i,1]],nirReflectance/spectrum[[i,2]]-1},{i,1,Length[spectrum]}]]]


(* ::Input::Initialization:: *)
getRandomPixelBlock[region_,blockEdgeInPixels_]:=Module[{currentIteration=0,maxIterations=100,block=With[{randomPoint=If[OddQ[blockEdgeInPixels],Round[RandomPoint[region]],Round[2*RandomPoint[region]]/2]},Rectangle[randomPoint-0.5{blockEdgeInPixels,blockEdgeInPixels},randomPoint+0.5{blockEdgeInPixels,blockEdgeInPixels}]]},While[!RegionWithin[region,block]&&currentIteration<=maxIterations,currentIteration=currentIteration+1;block=With[{randomPoint=Round[RandomPoint[region]]},Rectangle[randomPoint-0.5{blockEdgeInPixels,blockEdgeInPixels},randomPoint+0.5{blockEdgeInPixels,blockEdgeInPixels}]]];block]


(* ::Input::Initialization:: *)
triCube[x_]:=(Abs[Chop[If[ x<= 1,(1 - x^3)^3,0]]]);


(* ::Input::Initialization:: *)
countXsWithinHorizonRight[data_,center_,horizon_]:=Module[{count=0,x0=data[[center,1]]},Do[If[Abs[data[[i,1]]-x0]<=horizon,count++,Break[]],{i,center+1,Length[data]}];count]


(* ::Input::Initialization:: *)
countXsWithinHorizonLeft[data_,center_,horizon_]:=Module[{count=0,x0=data[[center,1]]},Do[If[Abs[data[[i,1]]-x0]<=horizon,count++,Break[]],{i,center-1,1,-1}];count]


(* ::Input::Initialization:: *)
countXsWithinHorizonLeftAndRightCompiled =Compile[{{data,_Real,2},{center,_Integer},{horizon,_Real}},With[{left=If[center>1,Module[{count=0,x0=data[[center,1]]},Do[If[Abs[data[[i,1]]-x0]<=horizon,count++,Break[]],{i,center-1,1,-1}];count],0],right=Module[{count=0,x0=data[[center,1]]},Do[If[Abs[data[[i,1]]-x0]<=horizon,count++,Break[]],{i,center+1,Length[data]}];count]},{left,right}],CompilationTarget->"C"];


(* ::Input::Initialization:: *)
LoessFit[data_,span_]:=Module[{startIndex,endIndex,weights,centerX,y},(Table[(startIndex = Max[1,j - span];endIndex = Min[Length[data], j+ span];centerX = data[[j,1]];weights = Map[(triCube[Abs[(centerX - #[[1]])/Max[Abs[centerX-data[[startIndex,1]]],Abs[centerX-data[[endIndex,1]]]]]])&,data[[startIndex;;endIndex]]];With[{modelFit=LinearModelFit[data[[startIndex;;endIndex]],{1,y,y^2},y,Weights->weights]},{centerX,modelFit[centerX]}]),{j,1,Length[data]}])]


(* ::Input::Initialization:: *)
cropGraphics[g_,x_,y_,w_,h_]:=Graphics[Inset[g,{x,y},{0,0}],PlotRange->{{0,1},{0,1}},ImageSize->{w,h},AspectRatio->Full]


(* ::Input::Initialization:: *)
getMeanReflectanceWithinRange[data_,wavelengthMin_,wavelengthMax_]:=Mean[Reap[Do[If[data[[i,1]]>=wavelengthMin&&data[[i,1]]<=wavelengthMax,Sow[data[[i,2]]]],{i,1,Length[data]}]][[2,1]]]


(* ::Input::Initialization:: *)
panelLabelsStyle=Directive[FontFamily->"Arial",9,SingleLetterItalics->False];


(* ::Input::Initialization:: *)
axisXStyleSmall=Directive[Black,FontFamily->"Arial",7,AbsoluteThickness[0.6]];


(* ::Input::Initialization:: *)
axisYStyleSmall=Directive[Black,FontFamily->"Arial",7,AbsoluteThickness[0.6]];


(* ::Input::Initialization:: *)
plotTextStyleSmall=Directive[FontFamily->"Arial",7];


(* ::Input::Initialization:: *)
axisXLabelLetteringStyleSmall=Directive[Black,FontFamily->"Arial",8];


(* ::Input::Initialization:: *)
axisYLabelLetteringStyleSmall=Directive[Black,FontFamily->"Arial",8];


(* ::Input::Initialization:: *)
plotLineThicknessSmall=AbsoluteThickness[0.6];


(* ::Input::Initialization:: *)
lineStyleALegend=Directive[plotLineThicknessSmall,Dashing[{0.01/1.4,0.013/1.4}],Black];


(* ::Input::Initialization:: *)
panelLetterStyleSmall=Directive[FontFamily->"Arial",8];


(* ::Input::Initialization:: *)
lineStyleA=Directive[plotLineThicknessSmall,Dashing[{0.01,0.013}],Black];
lineStyleB=Directive[plotLineThicknessSmall,Dashing[{.05,.02}],Black];
lineStyleC=Directive[plotLineThicknessSmall,Black,Dashing[{}]];
lineStyleD=Directive[plotLineThicknessSmall,Dashing[{0.01,0.04,0.05,0.04}],Black];
lineStyleE =With[{l1=0.025,l2=0.0075,s=0.025/2}, Directive[plotLineThicknessSmall,Dashing[{l1,s,l1,s,l2,s}],Black]];


(* ::Input::Initialization:: *)
diamond[color_]:=Graphics[{EdgeForm[{Black,AbsoluteThickness[0.5]}],FaceForm[{color,Opacity[1]}],Polygon[{{0,1},{1,0},{0,-1},{-1,0}}]}];


(* ::Input::Initialization:: *)
disk[color_]:=Graphics[{EdgeForm[{Black,AbsoluteThickness[0.5]}],FaceForm[color],Disk[]}];


(* ::Input::Initialization:: *)
rectangle[color_]:=Graphics[{EdgeForm[{Black,AbsoluteThickness[0.5]}],FaceForm[color],Rectangle[]}];


(* ::Input::Initialization:: *)
triangle[color_]:=Graphics[{EdgeForm[{Black,AbsoluteThickness[0.5]}],FaceForm[color],Polygon[{{-1,0},{1,0},{0,-2}}]}];


(* ::Input::Initialization:: *)
triangleUpwards[color_]:=Graphics[{EdgeForm[{Black,AbsoluteThickness[0.5]}],FaceForm[color],Polygon[{{-1,0},{1,0},{0,2}}]}]


(* ::Input::Initialization:: *)
star[color_]:=Graphics[{EdgeForm[{Black,AbsoluteThickness[0.5]}],FaceForm[color],Polygon[Table[If[OddQ[i],1,.4]*{Cos[2\[Pi]*i/10+54*Pi/180],Sin[2\[Pi]*i/10+54*Pi/180]},{i,1,10}]]}];


(* ::Input::Initialization:: *)
extractSampleCount[headerReadIn_]:=With[{keyValueString=Flatten[Select[headerReadIn,StringQ[#[[1]]]&&StringStartsQ[#[[1]],"samples"]&,1]][[1]]},ToExpression[StringSplit[keyValueString," = "][[2]]]]


(* ::Input::Initialization:: *)
extractLineCount[headerReadIn_]:=With[{keyValueString=Flatten[Select[headerReadIn,StringQ[#[[1]]]&&StringStartsQ[#[[1]],"lines"]&,1]][[1]]},ToExpression[StringSplit[keyValueString," = "][[2]]]]


(* ::Input::Initialization:: *)
extractBandCount[headerReadIn_]:=With[{keyValueString=Flatten[Select[headerReadIn,StringQ[#[[1]]]&&StringStartsQ[#[[1]],"bands"]&,1]][[1]]},ToExpression[StringSplit[keyValueString," = "][[2]]]]


(* ::Input::Initialization:: *)
getWavelengths[headerReadIn_]:=Flatten[With[{wavelengthStartIndex=1+Flatten[Position[headerReadIn,{"wavelength = {"},1]][[1]],bandCount=extractBandCount[headerReadIn]},Take[headerReadIn,{wavelengthStartIndex,wavelengthStartIndex+bandCount-1}]]]


(* ::Input::Initialization:: *)
(*Replaing Table with ParallelTable only slows down calculation*)readInHeadwallFileImages[hdrFile_]:=With[{headerReadIn=Import[hdrFile,"TSV"],headerFileNameSplit=FileNameSplit[hdrFile]},With[{dataFileName=FileNameJoin[Join[Take[headerFileNameSplit,{1,-2}],{FileBaseName[headerFileNameSplit[[-1]]]}]]},(With[{sampleCount=extractSampleCount[headerReadIn],lineCount=extractLineCount[headerReadIn],bandCount=extractBandCount[headerReadIn],wavelengthsInNanometers=getWavelengths[headerReadIn]},If[FileExistsQ[dataFileName],With[{rawData=Import[dataFileName,"UnsignedInteger16"]},createHyperspectralmageObject[Table[getImageData[rawData,bandIndex,bandCount,sampleCount],{bandIndex,0,bandCount-1}],wavelengthsInNanometers,bandCount,lineCount,sampleCount]],readInAndCreateImageDataFromCSV[dataFileName<>".csv",wavelengthsInNanometers,bandCount,lineCount,sampleCount]]])]]


(* ::Input::Initialization:: *)
(*Replaing Table with ParallelTable only slows down calculation*)extractSpectraFromHDRFile[hdrFile_]:=With[{headerReadIn=Import[hdrFile,"TSV"],headerFileNameSplit=FileNameSplit[hdrFile]},With[{dataFileName=FileNameJoin[Join[Take[headerFileNameSplit,{1,-2}],{FileBaseName[headerFileNameSplit[[-1]]]}]]},With[{sampleCount=extractSampleCount[headerReadIn],lineCount=extractLineCount[headerReadIn],bandCount=extractBandCount[headerReadIn],wavelengthsInNanometers=getWavelengths[headerReadIn]},With[{rawData=Import[dataFileName,"UnsignedInteger16"]},With[{imagesData=Table[getImageData[rawData,bandIndex,bandCount,sampleCount],{bandIndex,0,bandCount-1}]},extractSpectrum[imagesData,wavelengthsInNanometers,sampleCount,lineCount]]]]]]


(* ::Input::Initialization:: *)
(*Replaing Table with ParallelTable only slows down calculation*)readInHeadwallFileImages[hdrFile_,calibrationHDRFile_]:=With[{headerReadIn=Import[hdrFile,"TSV"],headerFileNameSplit=FileNameSplit[hdrFile]},With[{dataFileName=FileNameJoin[Join[Take[headerFileNameSplit,{1,-2}],{FileBaseName[headerFileNameSplit[[-1]]]}]]},(With[{sampleCount=extractSampleCount[headerReadIn],lineCount=extractLineCount[headerReadIn],bandCount=extractBandCount[headerReadIn],wavelengthsInNanometers=getWavelengths[headerReadIn]},If[FileExistsQ[dataFileName],With[{rawData=Import[dataFileName,"UnsignedInteger16"],calibrationSpectra=(*we shoudl flatten spectra across lines if there are more than 1*)extractSpectraFromHDRFile[calibrationHDRFile]},createHyperspectralmageObject[Table[getImageData[rawData,bandIndex,bandCount,sampleCount],{bandIndex,0,bandCount-1}],wavelengthsInNanometers,bandCount,lineCount,sampleCount,calibrationSpectra]],readInAndCreateImageDataFromCSV[dataFileName<>".csv",calibrationHDRFile,wavelengthsInNanometers,bandCount,lineCount,sampleCount]]])]]


(* ::Input::Initialization:: *)
readInAndCreateImageDataFromCSV[filePath_String,wavelengthsInNanometers_List,bandCount_?NumericQ,lineCount_?NumericQ,sampleCount_?NumericQ]:=With[{rawData=Import[filePath,"CSV"]},With[{imageDataAll=createImageDataFromCSV[rawData,bandCount,lineCount,sampleCount]},createHyperspectralmageObject[imageDataAll,wavelengthsInNanometers,bandCount,lineCount,sampleCount]]]


(* ::Input::Initialization:: *)
readInAndCreateImageDataFromCSV[filePath_String,calibrationHDRFile_String,wavelengthsInNanometers_List,bandCount_?NumericQ,lineCount_?NumericQ,sampleCount_?NumericQ]:=With[{rawData=Import[filePath,"CSV"]},With[{imageDataAll=createImageDataFromCSV[rawData,bandCount,lineCount,sampleCount],calibrationSpectra=(*we shoudl flatten spectra across lines if there are more than 1*)extractSpectraFromHDRFile[calibrationHDRFile]},createHyperspectralmageObject[imageDataAll,wavelengthsInNanometers,bandCount,lineCount,sampleCount,calibrationSpectra]]]


(* ::Input::Initialization:: *)
createImageDataFromCSV[readIn_List,bandCount_?NumericQ,lineCount_?NumericQ,sampleCount_?NumericQ]:=Module[{imageData=ConstantArray[0,{bandCount,lineCount,sampleCount}]},Do[With[{line=readIn[[i,1]],sample=readIn[[i,2]]},Do[imageData[[j-2,line+1,sample+1]]=readIn[[i,j]],{j,3,bandCount+2}]],{i,2,Length[readIn]}];imageData]


(* ::Input::Initialization:: *)
createHyperspectralmageObject[imagesData_List,wavelengthsInNanometers_List,bandCount_Integer,lineCount_Integer,sampleCount_Integer,darkCalibrationData_:Null]:=With[{},Module[{self},self=Association[
"images"->imagesData,
"wavelengthsInNanometers"->wavelengthsInNanometers,
"bandCount"->bandCount,
"sampleCount"->sampleCount,
"lineCount"->lineCount,
"darkCalibrationData"->darkCalibrationData,
getLineCount->(lineCount&),
getSpectralBandCount->(bandCount&),
getWavelengthsInNanometers->(wavelengthsInNanometers&),
getImageDimensions->(Dimensions[imagesData[[1]]]&),
getImageDataByIndex->(imagesData[[#]]&),calculateChloroplastMovementIndex->(With[{wavelength635IndexCeiling=Ceiling[getIndexOfWavelength[wavelengthsInNanometers,635]],wavelength635IndexFloor=Floor[getIndexOfWavelength[wavelengthsInNanometers,635]],wavelength555IndexCeiling=Ceiling[getIndexOfWavelength[wavelengthsInNanometers,555]],wavelength555IndexFloor=Floor[getIndexOfWavelength[wavelengthsInNanometers,555]],referenceSpectrum=self[extracRawSpectrum][#2]},(*extractRawCountsInRegionForSelectedWavelengthsPrivate[imagesData,wavelengthsInNanometers,{wavelength635IndexCeiling,wavelength635IndexFloor,wavelength555IndexCeiling,wavelength555IndexFloor},#,sampleCount,lineCount,darkCalibrationData]*)With[{countsForSelectedWavelengths=extractRawCountsInRegionForSelectedWavelengthsPrivate[imagesData,{wavelength635IndexCeiling,wavelength635IndexFloor,wavelength555IndexCeiling,wavelength555IndexFloor},#1,sampleCount,lineCount,darkCalibrationData]},Mean[Map[(N[(#[[1]]/referenceSpectrum[[wavelength635IndexCeiling,2]]+#[[2]]/referenceSpectrum[[wavelength635IndexFloor,2]]-#[[3]]/referenceSpectrum[[wavelength555IndexCeiling,2]]-#[[4]]/referenceSpectrum[[wavelength555IndexFloor,2]])/(#[[1]]/referenceSpectrum[[wavelength635IndexCeiling,2]]+#[[2]]/referenceSpectrum[[wavelength635IndexFloor,2]]+#[[3]]/referenceSpectrum[[wavelength555IndexCeiling,2]]+#[[4]]/referenceSpectrum[[wavelength555IndexFloor,2]])]&),countsForSelectedWavelengths]](*Mean[Map[(#[[1]]+#[[2]]-#[[3]]-#[[4]])/(#[[1]]+#[[2]]+#[[3]]+#[[4]]),countsForSelectedWavelengths]]*)]]&),
displayChloroplastMovementIndexHistogram->(Histogram[Flatten[createMatrixOfChloroplastMovementIndex[imagesData,#1,#2,#3,#4,sampleCount,lineCount]]]&),
(*https://mathematica.stackexchange.com/questions/89422/how-can-i-change-the-thickness-of-tick-marks-in-barlegend*)displayChloroplastMovementIndex->(ArrayPlot[With[{},createMatrixOfChloroplastMovementIndex[imagesData,sampleCount,lineCount]],PlotRange->{-0.4,0.1},ColorFunction->(ColorData["TemperatureMap"][#]&),PlotLegends->BarLegend[Automatic,Method->{Frame->True,Ticks->{{-0.4,"0.4"},{-0.2,"-0.2"},{0.0,"0.0"},{0.1,"0.1"}},TicksStyle->Directive[AbsoluteThickness[0.5]]},LegendFunction->(#/. Directive[AbsoluteThickness[_]]->Directive[AbsoluteThickness[0.5],Opacity[1]]&),LegendLayout->"Row",LegendMarkerSize->200,LegendLabel->"ND\[ThinSpace]\[ThinSpace](635 nm, 555 nm)",LabelStyle->Directive[Black,FontFamily->"Arial",FontSize->14]],Frame->False]&),displayChloroplastMovementIndexImageCropped->(ArrayPlot[With[{},createMatrixOfChloroplastMovementIndex[imagesData,#1,#2,#3,#4,sampleCount,lineCount]],PlotRange->{-0.4,0.1},ColorFunction->(ColorData["TemperatureMap"][#]&),PlotLegends->BarLegend[Automatic,Method->{Frame->True,Ticks->{{-0.4,"0.4"},{-0.2,"-0.2"},{0.0,"0.0"},{0.1,"0.1"}},TicksStyle->Directive[AbsoluteThickness[0.5]]},LegendFunction->(#/. Directive[AbsoluteThickness[_]]->Directive[AbsoluteThickness[0.5],Opacity[1]]&),LegendLayout->"Row",LegendMarkerSize->200,LegendLabel->"ND\[ThinSpace]\[ThinSpace](635 nm, 555 nm)",LabelStyle->Directive[Black,FontFamily->"Arial",FontSize->14]],Frame->False]&),displayChloroplastMovementIndexImageCroppedWithoutLegend->(ArrayPlot[With[{},createMatrixOfChloroplastMovementIndex[imagesData,#1,#2,#3,#4,sampleCount,lineCount]],PlotRange->{-0.4,0.1},ColorFunction->(ColorData["TemperatureMap"][#]&),Frame->False]&),displayChloroplastMovementIndexImageCroppedWithoutLegendAndBackground->(ArrayPlot[With[{},createMatrixOfChloroplastMovementIndexAndHideBackground[imagesData,#1,#2,#3,#4,sampleCount,lineCount,#5,#6]],PlotRange->{-0.4,0.1},ColorFunction->(ColorData["TemperatureMap"][#]&),Frame->False]&),

getVISReflectanceImage->(With[{meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#,sampleCount,lineCount,darkCalibrationData]]},createMatrixOfVISReflectance[imagesData,wavelengthsInNanometers,meansReference,sampleCount,lineCount,darkCalibrationData]]&),
displayVISReflectanceImage->(Image[With[{meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]]},createMatrixOfVISReflectance[imagesData,wavelengthsInNanometers,meansReference,sampleCount,lineCount,darkCalibrationData]],ColorSpace->"Grayscale"]&),
displayVISReflectanceImageCropped->(ArrayPlot[With[{meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]]},createMatrixOfVISReflectance[imagesData,#2,#3,#4,#5,wavelengthsInNanometers,meansReference,sampleCount,lineCount,darkCalibrationData]],PlotRange->{0,0.3},ColorFunction->(GrayLevel[#]&),PlotLegends->BarLegend[Automatic,LegendLayout->"Row",LegendMarkerSize->150,LegendLabel->"VIS Reflectance",LabelStyle->Directive[Black,FontSize->15]],Frame->False]&),
displayReflectanceImageCropped->(ArrayPlot[With[{meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]]},createMatrixOfReflectanceArbitraryWavelengthRange[imagesData,#6,#7,#2,#3,#4,#5,wavelengthsInNanometers,meansReference,sampleCount,lineCount,darkCalibrationData]],PlotRange->{0,0.3},ColorFunction->(GrayLevel[#]&),PlotLegends->BarLegend[Automatic,Method->{Frame->False,TicksStyle->Directive[AbsoluteThickness[2]]},LegendLayout->"Row",LegendMarkerSize->150,LegendLabel->"VIS Reflectance",LabelStyle->Directive[Black,FontSize->15]],Frame->False]&),displayNDVImageCropped->(ArrayPlot[With[{meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]]},createMatrixOfNDVIndex[imagesData,#2,#3,#4,#5,wavelengthsInNanometers,meansReference,sampleCount,lineCount,darkCalibrationData]],ColorFunction->(GrayLevel[#]&),PlotRange->{0,1},PlotLegends->BarLegend[Automatic,LegendMarkerSize->150,LegendLabel->"NDVI",LabelStyle->Directive[Black,FontSize->15]],Frame->False]&),
displayImageByIndex->(Image[imagesData[[#]]/4095.,ColorSpace->"Grayscale"]&),
displayImageByIndexCustomNormalized->(Image[imagesData[[#1]]/#2,ColorSpace->"Grayscale"]&),
getImageDataByWavelength->(imagesData[[getIndexOfWavelength[wavelengthsInNanometers,#]]]&),
displayImageDataByWavelength->((Image[getImageForWavelengthWithoutInterpolation[imagesData,wavelengthsInNanometers,bandCount,sampleCount,#1]/4095.,ColorSpace->"Grayscale"])&),displayImageDataByWavelengthConstrastAdjusted->((Image[getImageForWavelengthWithoutInterpolation[imagesData,wavelengthsInNanometers,bandCount,sampleCount,#1]/#2,ColorSpace->"Grayscale"])&),extractRawCountsInRegionForEachWavelength->(extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#,sampleCount,lineCount]&),extractMeanSpectrumWithinWholeImage->(calculateMeansOfSpectra[extractSpectrum[imagesData,wavelengthsInNanometers,sampleCount,lineCount]]&),extractSpectraWithinWholeImage->(extractSpectrum[imagesData,wavelengthsInNanometers,sampleCount,lineCount]&),extracRawSpectrumWithoutDarkCalibration->(calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#,sampleCount,lineCount]]&),extracRawSpectrum->(Quiet[calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#,sampleCount,lineCount,darkCalibrationData]]]&),extractReflectanceSpectrumWithoutDarkCalibration->(With[{meansRegionOfInterest=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount]],meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount]]},Table[{meansRegionOfInterest[[i,1]],meansRegionOfInterest[[i,2]]/meansReference[[i,2]]},{i,2,Length[meansRegionOfInterest]}]]&),extractReflectanceSpectrum->(With[{meansRegionOfInterest=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]],meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount,darkCalibrationData]]},Table[{meansRegionOfInterest[[i,1]],meansRegionOfInterest[[i,2]]/meansReference[[i,2]]},{i,2,Length[meansRegionOfInterest]}]]&),extractReflectanceWithinRegionUsingWhiteReferenceCurve->(With[{meansRegionOfInterest=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]],referenceFunction=Interpolation[#2]},Table[{meansRegionOfInterest[[i,1]],meansRegionOfInterest[[i,2]]/referenceFunction[meansRegionOfInterest[[i,1]]]},{i,2,Length[meansRegionOfInterest]}]]&),
extractDifferenceInReflectanceBetweenRegionsUsingReferenceWithoutDarkCalibration->(With[{meansRegionOfInterestA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount]],meansRegionOfInterestB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount]],meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#3,sampleCount,lineCount]]},Table[{meansRegionOfInterestA[[i,1]],(meansRegionOfInterestA[[i,2]]-meansRegionOfInterestB[[i,2]])/meansReference[[i,2]]},{i,2,Length[meansRegionOfInterestA]}]]&),extractDifferenceInReflectenceSpectrum->(With[{meansRegionOfInterestA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]],meansRegionOfInterestB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount,darkCalibrationData]],meansReference=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#3,sampleCount,lineCount,darkCalibrationData]]},Table[{meansRegionOfInterestA[[i,1]],(meansRegionOfInterestA[[i,2]]-meansRegionOfInterestB[[i,2]])/meansReference[[i,2]]},{i,2,Length[meansRegionOfInterestA]}]]&),extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurve->((With[{meansRegionOfInterestA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount]],meansRegionOfInterestB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount]],referenceFunction=Interpolation[#3]},Table[{meansRegionOfInterestA[[i,1]],(meansRegionOfInterestA[[i,2]]-meansRegionOfInterestB[[i,2]])/referenceFunction[meansRegionOfInterestA[[i,1]]]},{i,2,Length[meansRegionOfInterestA]}]])&),
extractDifferenceInReflectanceSpectrumUsingWhiteReferenceCurveWithCalibration->((With[{meansRegionOfInterestA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]],meansRegionOfInterestB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount,darkCalibrationData]],referenceFunction=Interpolation[#3]},Table[{meansRegionOfInterestA[[i,1]],(meansRegionOfInterestA[[i,2]]-meansRegionOfInterestB[[i,2]])/referenceFunction[meansRegionOfInterestA[[i,1]]]},{i,2,Length[meansRegionOfInterestA]}]])&),
extractDifferenceInMeanSpectrumBetweenRegionsWithRespectToTwoReferences->((With[{meansRegionOfInterestA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount]],meansRegionOfInterestB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount]],meansReferenceA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#3,sampleCount,lineCount]],meansReferenceB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#4,sampleCount,lineCount]]},Table[{meansRegionOfInterestA[[i,1]],meansRegionOfInterestA[[i,2]]/meansReferenceA[[i,2]]-meansRegionOfInterestB[[i,2]]/meansReferenceB[[i,2]]},{i,2,Length[meansRegionOfInterestA]}]])&),
extractDifferenceInReflectanceBetweenRegionsUsingTwoWhiteReferenceRegions->((With[{meansRegionOfInterestA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#1,sampleCount,lineCount,darkCalibrationData]],meansRegionOfInterestB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#2,sampleCount,lineCount,darkCalibrationData]],meansReferenceA=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#3,sampleCount,lineCount,darkCalibrationData]],meansReferenceB=calculateMeansOfSpectra[extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthsInNanometers,#4,sampleCount,lineCount,darkCalibrationData]]},Table[{meansRegionOfInterestA[[i,1]],meansRegionOfInterestA[[i,2]]/meansReferenceA[[i,2]]-meansRegionOfInterestB[[i,2]]/meansReferenceB[[i,2]]},{i,2,Length[meansRegionOfInterestA]}]])&),
classifyChloroplastResponse->(With[{differenceSpectrum=self[extractDifferenceInReflectenceSpectrum][#1,#2,#3],labels={"Accumulation","Avoidance","Dark"}},With[{classification=convolutionalClassifier[Take[Map[#[[2]]&,differenceSpectrum],{1,158}]]},labels[[classification]]]]&)]]]


(* ::Input::Initialization:: *)
extractSpectrumWithinRegionOld1[imagesData_List,wavelengthInNanoMeters_List,region_?RegionQ,imageSizeX_Integer,imageSizeY_Integer]:=With[{boundingBox=RegionBounds[region]},Module[{spectralData=Table[{wavelengthInNanoMeters[[i]],{}},{i,1,Length[wavelengthInNanoMeters]}]},Do[(*If[IntegerQ[y/20]&&IntegerQ[x/20],Print["x "<>ToString[x]<>", y "<>ToString[y]]];*)If[RegionMember[region,{x,y}],Do[spectralData[[bandIndex,2]]=Append[spectralData[[bandIndex,2]],imagesData[[bandIndex,y,x]]],{bandIndex,1,Length[imagesData]}]],{y,Max[1,Floor[boundingBox[[2,1]]]],Min[imageSizeY,Ceiling[boundingBox[[2,2]]]]},{x,Max[1,Floor[boundingBox[[1,1]]]],Min[imageSizeX,Ceiling[boundingBox[[1,2]]]]}];spectralData]]


(* ::Input::Initialization:: *)
extractSpectrumWithinRegionOld2[imagesData_List,wavelengthInNanoMeters_List,region_?RegionQ,imageSizeX_Integer,imageSizeY_Integer]:=With[{boundingBox=RegionBounds[region]},Reap[Do[(*If[IntegerQ[y/20]&&IntegerQ[x/20],Print["x "<>ToString[x]<>", y "<>ToString[y]]];*)If[RegionMember[region,{x,y}],Do[Sow[imagesData[[bandIndex,y,x]],bandIndex],{bandIndex,1,Length[imagesData]}]],{y,Max[1,Floor[boundingBox[[2,1]]]],Min[imageSizeY,Ceiling[boundingBox[[2,2]]]]},{x,Max[1,Floor[boundingBox[[1,1]]]],Min[imageSizeX,Ceiling[boundingBox[[1,2]]]]}],Range[1,Length[imagesData]],(wavelengthInNanoMeters[[#]]&)]]


(* ::Input::Initialization:: *)
extractSpectrumWithinRegionOld3[imagesData_List,wavelengthInNanoMeters_List,region_?RegionQ,imageSizeX_Integer,imageSizeY_Integer]:=With[{boundingBox=RegionBounds[region]},With[{minX=Max[1,Floor[boundingBox[[1,1]]]],maxX=Min[imageSizeX,Ceiling[boundingBox[[1,2]]]],minY=Max[1,Floor[boundingBox[[2,1]]]],maxY=Min[imageSizeY,Ceiling[boundingBox[[2,2]]]]},With[{regionInsidness=ParallelTable[RegionMember[region,{x,y}],{y,minY,maxY},{x,minX,minY}]},Reap[Do[(*If[IntegerQ[y/20]&&IntegerQ[x/20],Print["x "<>ToString[x]<>", y "<>ToString[y]]];*)If[regionInsidness[[y-minY+1,x-minX+1]],Do[Sow[imagesData[[bandIndex,y,x]],bandIndex],{bandIndex,1,Length[imagesData]}]],{y,minY,maxY},{x,minX,maxX}],Range[1,Length[imagesData]],(wavelengthInNanoMeters[[#]]&)]]]]


(* ::Input::Initialization:: *)
extractRawCountsInRegionForEachWavelengthPrivate[imagesData_List,wavelengthInNanoMeters_List,region_?RegionQ,imageSizeX_Integer,imageSizeY_Integer]:=With[{boundingBox=RegionBounds[region],regionMemberFunction=RegionMember[region]},With[{minX=Max[1,Floor[boundingBox[[1,1]]]],maxX=Min[imageSizeX,Ceiling[boundingBox[[1,2]]]],minY=Max[1,imageSizeY-Ceiling[boundingBox[[2,2]]]],maxY=Min[imageSizeY,imageSizeY-Floor[boundingBox[[2,1]]]]},With[{regionInsidness=ParallelTable[regionMemberFunction[{x,imageSizeY-y}],{y,minY,maxY},{x,minX,maxX}]},Flatten[Reap[Do[(*If[IntegerQ[y/20]&&IntegerQ[x/20],Print["x "<>ToString[x]<>", y "<>ToString[y]]];*)If[regionInsidness[[y-minY+1,x-minX+1]],Do[Sow[imagesData[[bandIndex,y,x]],bandIndex],{bandIndex,1,Length[imagesData]}]],{y,minY,maxY},{x,minX,maxX}],Range[1,Length[imagesData]],({wavelengthInNanoMeters[[#1]],#2}&)][[2]],1]]]]


(* ::Input::Initialization:: *)
extractRawCountsInRegionForEachWavelengthPrivate[imagesData_List,wavelengthInNanoMeters_List,region_?RegionQ,imageSizeX_Integer,imageSizeY_Integer,calibrationSpectra_]:=If[calibrationSpectra===Null,extractRawCountsInRegionForEachWavelengthPrivate[imagesData,wavelengthInNanoMeters,region,imageSizeX,imageSizeY],With[{boundingBox=RegionBounds[region],regionMemberFunction=RegionMember[region]},With[{minX=Max[1,Floor[boundingBox[[1,1]]]],maxX=Min[imageSizeX,Ceiling[boundingBox[[1,2]]]],minY=Max[1,imageSizeY-Ceiling[boundingBox[[2,2]]]],maxY=Min[imageSizeY,imageSizeY-Floor[boundingBox[[2,1]]]]},With[{regionInsidness=ParallelTable[regionMemberFunction[{x,imageSizeY-y}],{y,minY,maxY},{x,minX,maxX}]},Flatten[Reap[Do[(*If[IntegerQ[y/20]&&IntegerQ[x/20],Print["x "<>ToString[x]<>", y "<>ToString[y]]];*)If[regionInsidness[[y-minY+1,x-minX+1]],Do[Sow[imagesData[[bandIndex,y,x]]-calibrationSpectra[[bandIndex,2,x]],bandIndex],{bandIndex,1,Length[imagesData]}]],{y,minY,maxY},{x,minX,maxX}],Range[1,Length[imagesData]],({wavelengthInNanoMeters[[#1]],#2}&)][[2]],1]]]]]


(* ::Input::Initialization:: *)
extractRawCountsInRegionForSelectedWavelengthsPrivate[imagesData_List,bandIndices_List,region_?RegionQ,imageSizeX_Integer,imageSizeY_Integer,calibrationSpectra_]:=With[{boundingBox=RegionBounds[region],regionMemberFunction=RegionMember[region]},With[{minX=Max[1,Floor[boundingBox[[1,1]]]],maxX=Min[imageSizeX,Ceiling[boundingBox[[1,2]]]],minY=Max[1,imageSizeY-Ceiling[boundingBox[[2,2]]]],maxY=Min[imageSizeY,imageSizeY-Floor[boundingBox[[2,1]]]]},With[{regionInsidness=ParallelTable[regionMemberFunction[{x,imageSizeY-y}],{y,minY,maxY},{x,minX,maxX}]},Flatten[Reap[Do[(*If[IntegerQ[y/20]&&IntegerQ[x/20],Print["x "<>ToString[x]<>", y "<>ToString[y]]];*)If[regionInsidness[[y-minY+1,x-minX+1]],Sow[Table[imagesData[[bandIndex,y,x]]-calibrationSpectra[[bandIndex,2,x]],{bandIndex,bandIndices}]]],{y,minY,maxY},{x,minX,maxX}]][[2]],1]]]]


(* ::Input::Initialization:: *)
extractSpectrum[imagesData_List,wavelengthInNanoMeters_List,imageSizeX_Integer,imageSizeY_Integer]:=With[{minX=1,maxX=imageSizeX,minY=1,maxY=imageSizeY},Flatten[Reap[Do[(*If[IntegerQ[y/20]&&IntegerQ[x/20],Print["x "<>ToString[x]<>", y "<>ToString[y]]];*)Do[Sow[imagesData[[bandIndex,y,x]],bandIndex],{bandIndex,1,Length[imagesData]}],{y,minY,maxY},{x,minX,maxX}],Range[1,Length[imagesData]],({wavelengthInNanoMeters[[#1]],#2}&)][[2]],1]]


(* ::Input::Initialization:: *)
createMatrixOfSpectra[imagesData_List,wavelengthInNanoMeters_List,imageSizeX_Integer,imageSizeY_Integer]:=With[{minX=1,maxX=imageSizeX,minY=1,maxY=imageSizeY},Table[Table[{wavelengthInNanoMeters[[bandIndex]],imagesData[[bandIndex,y,x]]},{bandIndex,1,Length[imagesData]}],{y,minY,maxY},{x,minX,maxX}]]


(* ::Input::Initialization:: *)
createMatrixOfSpectraWithRespectToReferenceWithCalibration[imagesData_List,wavelengthInNanoMeters_List,referenceSpectrum_,imageSizeX_Integer,imageSizeY_Integer,calibrationSpectra_]:=With[{minX=1,maxX=imageSizeX,minY=1,maxY=imageSizeY},Table[Table[{wavelengthInNanoMeters[[bandIndex]],(imagesData[[bandIndex,y,x]]-calibrationSpectra[[bandIndex,2,x]])/referenceSpectrum[[bandIndex,2]]},{bandIndex,1,Length[imagesData]}],{y,minY,maxY},{x,minX,maxX}]]


(* ::Input::Initialization:: *)
createMatrixOfChloroplastMovementIndex[imagesData_List,imageSizeX_Integer,imageSizeY_Integer]:=With[{minX=1,maxX=imageSizeX,minY=1,maxY=imageSizeY},createMatrixOfChloroplastMovementIndex[imagesData,minX,maxX,minY,maxY,imageSizeX,imageSizeY]]


(* ::Input::Initialization:: *)
createMatrixOfChloroplastMovementIndex[imagesData_List,minX_,maxX_,minY_,maxY_,imageSizeX_Integer,imageSizeY_Integer]:=With[{wavelength635IndexCeiling=Ceiling[getIndexOfWavelength[635]],wavelength635IndexFloor=Floor[getIndexOfWavelength[635]],wavelength555IndexCeiling=Ceiling[getIndexOfWavelength[555]],wavelength555IndexFloor=Floor[getIndexOfWavelength[555]]},Table[(imagesData[[wavelength635IndexCeiling,y,x]]+imagesData[[wavelength635IndexFloor,y,x]]-imagesData[[wavelength555IndexCeiling,y,x]]-imagesData[[wavelength555IndexFloor,y,x]])/(imagesData[[wavelength635IndexCeiling,y,x]]+imagesData[[wavelength635IndexFloor,y,x]]+imagesData[[wavelength555IndexCeiling,y,x]]+imagesData[[wavelength555IndexFloor,y,x]]),{y,minY,maxY},{x,minX,maxX}]]


(* ::Input::Initialization:: *)
createMatrixOfChloroplastMovementIndexAndHideBackground[imagesData_List,minX_,maxX_,minY_,maxY_,imageSizeX_Integer,imageSizeY_Integer,minGreenReflectance_,backgroundCode_]:=With[{wavelength635IndexCeiling=Ceiling[getIndexOfWavelength[635]],wavelength635IndexFloor=Floor[getIndexOfWavelength[635]],wavelength555IndexCeiling=Ceiling[getIndexOfWavelength[555]],wavelength555IndexFloor=Floor[getIndexOfWavelength[555]]},Table[If[imagesData[[wavelength555IndexCeiling,y,x]]>=minGreenReflectance,(imagesData[[wavelength635IndexCeiling,y,x]]+imagesData[[wavelength635IndexFloor,y,x]]-imagesData[[wavelength555IndexCeiling,y,x]]-imagesData[[wavelength555IndexFloor,y,x]])/(imagesData[[wavelength635IndexCeiling,y,x]]+imagesData[[wavelength635IndexFloor,y,x]]+imagesData[[wavelength555IndexCeiling,y,x]]+imagesData[[wavelength555IndexFloor,y,x]]),backgroundCode],{y,minY,maxY},{x,minX,maxX}]]


(* ::Input::Initialization:: *)
createMatrixOfVISReflectance[imagesData_List,wavelengthInNanoMeters_List,referenceSpectrum_,imageSizeX_Integer,imageSizeY_Integer,calibrationSpectra_]:=With[{minX=1,maxX=imageSizeX,minY=1,maxY=imageSizeY},createMatrixOfVISReflectance[imagesData,minX,maxX,minY,maxY,wavelengthInNanoMeters,referenceSpectrum,imageSizeX,imageSizeY,calibrationSpectra]]


(* ::Input::Initialization:: *)
createMatrixOfVISReflectance[imagesData_List,minX_,maxX_,minY_,maxY_,wavelengthInNanoMeters_List,referenceSpectrum_,imageSizeX_Integer,imageSizeY_Integer,calibrationSpectra_]:=With[{minWavelengthIndex=Max[1,Round[getIndexOfWavelength[400]]],maxWavelengthIndex=Min[Length[wavelengthInNanoMeters],Round[getIndexOfWavelength[700]]]},Table[Max[0,Mean[Table[(imagesData[[bandIndex,y,x]]-calibrationSpectra[[bandIndex,2,x]])/referenceSpectrum[[bandIndex,2]],{bandIndex,minWavelengthIndex,maxWavelengthIndex}]]],{y,minY,maxY},{x,minX,maxX}]]


(* ::Input::Initialization:: *)
createMatrixOfReflectanceArbitraryWavelengthRange[imagesData_List,minWavelength_,maxWavelength_,minX_,maxX_,minY_,maxY_,wavelengthInNanoMeters_List,referenceSpectrum_,imageSizeX_Integer,imageSizeY_Integer,calibrationSpectra_]:=With[{minWavelengthIndex=Max[1,Round[getIndexOfWavelength[minWavelength]]],maxWavelengthIndex=Min[Length[wavelengthInNanoMeters],Round[getIndexOfWavelength[maxWavelength]]]},Table[Max[0,Mean[Table[(imagesData[[bandIndex,y,x]]-calibrationSpectra[[bandIndex,2,x]])/referenceSpectrum[[bandIndex,2]],{bandIndex,minWavelengthIndex,maxWavelengthIndex}]]],{y,minY,maxY},{x,minX,maxX}]]


(* ::Input::Initialization:: *)
createMatrixOfNDVIndex[imagesData_List,minX_,maxX_,minY_,maxY_,wavelengthInNanoMeters_List,referenceSpectrum_,imageSizeX_Integer,imageSizeY_Integer,calibrationSpectra_]:=With[{redIndex=Max[1,Round[getIndexOfWavelength[680]]],NRIIndex=Min[Length[wavelengthInNanoMeters],Round[getIndexOfWavelength[800]]]},Table[With[{reflectanceRed=(imagesData[[redIndex,y,x]]-calibrationSpectra[[redIndex,2,x]])/referenceSpectrum[[redIndex,2]],reflectanceNIR=(imagesData[[NRIIndex,y,x]]-calibrationSpectra[[NRIIndex,2,x]])/referenceSpectrum[[NRIIndex,2]]},(reflectanceNIR-reflectanceRed)/(reflectanceNIR+reflectanceRed)],{y,minY,maxY},{x,minX,maxX}]]


(* ::Input::Initialization:: *)
calculateMeansOfSpectra[spectra_List]:=Table[{spectra[[i,1]],N[Mean[spectra[[i,2]]]]},{i,1,Length[spectra]}]


(* ::Input::Initialization:: *)
wavelengths={398.523999314,400.757190915,402.990382516,405.2235741170001,407.456765718,409.689957319,411.92314892,414.156340521,416.389532122,418.622723723,420.855915324,423.089106925,425.3222985260001,427.555490127,429.788681728,432.021873329,434.25506493,436.488256531,438.721448132,440.954639733,443.187831334,445.4210229350001,447.654214536,449.887406137,452.120597738,454.353789339,456.58698094,458.820172541,461.053364142,463.286555743,465.5197473440001,467.752938945,469.986130546,472.219322147,474.4525137480001,476.685705349,478.91889695,481.152088551,483.385280152,485.6184717530001,487.851663354,490.084854955,492.318046556,494.5512381570001,496.784429758,499.017621359,501.25081296,503.484004561,505.7171961620001,507.950387763,510.183579364,512.4167709650001,514.6499625660001,516.883154167,519.116345768,521.349537369,523.5827289700001,525.8159205710001,528.049112172,530.282303773,532.515495374,534.7486869750001,536.9818785760001,539.215070177,541.448261778,543.681453379,545.91464498,548.1478365810001,550.381028182,552.614219783,554.847411384,557.080602985,559.3137945860001,561.546986187,563.7801777880001,566.013369389,568.24656099,570.4797525910001,572.7129441920001,574.9461357930001,577.179327394,579.412518995,581.6457105960001,583.8789021970001,586.112093798,588.345285399,590.578477,592.8116686010001,595.0448602020001,597.278051803,599.511243404,601.744435005,603.9776266060001,606.2108182070001,608.444009808,610.677201409,612.91039301,615.1435846110001,617.3767762120001,619.609967813,621.843159414,624.076351015,626.309542616,628.5427342170001,630.775925818,633.0091174190001,635.24230902,637.475500621,639.7086922220001,641.9418838230001,644.1750754240001,646.408267025,648.641458626,650.8746502270001,653.1078418280001,655.341033429,657.57422503,659.807416631,662.0406082320001,664.2737998330001,666.5069914340002,668.740183035,670.973374636,673.2065662370001,675.4397578380001,677.6729494390001,679.90614104,682.139332641,684.3725242420001,686.6057158430001,688.8389074440001,691.072099045,693.305290646,695.538482247,697.7716738480001,700.0048654490001,702.23805705,704.471248651,706.704440252,708.9376318530001,711.1708234540001,713.4040150549999,715.637206656,717.870398257,720.1035898580001,722.3367814590001,724.5699730599999,726.803164661,729.036356262,731.2695478630001,733.5027394640001,735.7359310650002,737.969122666,740.202314267,742.4355058680001,744.6686974690001,746.9018890700002,749.135080671,751.368272272,753.6014638730001,755.8346554740001,758.0678470750001,760.301038676,762.534230277,764.7674218780001,767.0006134790001,769.2338050800001,771.466996681,773.700188282,775.933379883,778.1665714840001,780.3997630850001,782.6329546859999,784.866146287,787.099337888,789.3325294890001,791.5657210900001,793.7989126909999,796.032104292,798.265295893,800.4984874940001,802.7316790950001,804.9648706960002,807.198062297,809.431253898,811.6644454990001,813.8976371000001,816.1308287010002,818.364020302,820.597211903,822.8304035040001,825.0635951050001,827.2967867060001,829.529978307,831.763169908,833.9963615090001,836.2295531100001,838.4627447110001,840.695936312,842.929127913,845.162319514,847.3955111150001,849.6287027160001,851.861894317,854.095085918,856.328277519,858.5614691200001,860.7946607210001,863.0278523219999,865.261043923,867.494235524,869.7274271250001,871.9606187260001,874.1938103270002,876.427001928,878.660193529,880.8933851300001,883.1265767310001,885.3597683320002,887.592959933,889.826151534,892.0593431350001,894.2925347360001,896.5257263370002,898.758917938,900.992109539,903.2253011400001,905.4584927410001,907.6916843420001,909.924875943,912.158067544,914.391259145,916.6244507460001,918.8576423470001,921.090833948,923.324025549,925.55721715,927.7904087510001,930.0236003520001,932.2567919529999,934.489983554,936.723175155,938.9563667560001,941.1895583570001,943.4227499580002,945.655941559,947.88913316,950.1223247610001,952.3555163620001,954.5887079630002,956.821899564,959.055091165,961.2882827660001,963.5214743670001,965.7546659680002,967.987857569,970.22104917,972.4542407710001,974.6874323720001,976.9206239730001,979.153815574,981.387007175,983.6201987760001,985.8533903770001,988.0865819780001,990.319773579,992.55296518,994.786156781,997.0193483820001,999.2525399830001,1001.485731584,1003.718923185};


(* ::Input::Initialization:: *)
getIndexOfWavelength[wavelengthInNanoMeters_?Positive]:=BinarySearch[wavelengths,wavelengthInNanoMeters]


(* ::Input::Initialization:: *)
getIndexOfWavelength[wavelengths_List,wavelengthInNanoMeters_?Positive]:=BinarySearch[wavelengths,wavelengthInNanoMeters]


(* ::Input::Initialization:: *)
getImageData[rawData_,(*from 0 to totalBandCount - 1*)bandIndex_,totalBandCount_,sampleCount_]:=With[{dataLength=Length[rawData]},With[{lineCount=Length[rawData]/(totalBandCount*sampleCount)},Table[Take[rawData,{1+lineIndex*sampleCount*totalBandCount+bandIndex*sampleCount,Min[1+lineIndex*sampleCount*totalBandCount+(bandIndex+1)*sampleCount,dataLength]}],{lineIndex,0,lineCount-1}]]]


(* ::Input::Initialization:: *)
extractImageForWavelengthWithInterpolation[allImages_,wavelengthInNanoMeters_]:=With[{index1=Floor[getIndexOfWavelength[wavelengthInNanoMeters]],index2=Ceiling[getIndexOfWavelength[wavelengthInNanoMeters]]},With[{image1=allImages[[index1]],image2=allImages[[index2]]},With[{d1=wavelengthInNanoMeters-wavelengths[[index1]],d2=wavelengths[[index2]]-wavelengthInNanoMeters},Table[(1-d1/(d1+d2))image1[[i]]+image1[[2]](d1/(d1+d2)),{i,1,Length[image1]}]]]]


(* ::Input::Initialization:: *)
getImageForWavelengthWithInterpolation[rawData_,totalBandCount_,sampleCount_,wavelengthInNanoMeters_]:=With[{index1=Floor[getIndexOfWavelength[wavelengthInNanoMeters]],index2=Ceiling[getIndexOfWavelength[wavelengthInNanoMeters]]},With[{image1=getImageData[rawData,index1-1,totalBandCount,sampleCount],image2=getImageData[rawData,index2-1,totalBandCount,sampleCount]},With[{d1=wavelengthInNanoMeters-wavelengths[[index1]],d2=wavelengths[[index2]]-wavelengthInNanoMeters},Table[(1-d1/(d1+d2))image1[[i]]+image2[[i]](d1/(d1+d2)),{i,1,Length[image1]}]]]]


(* ::Input::Initialization:: *)
getImageForWavelengthWithInterpolation[images_List,wavelengths_List,totalBandCount_?NumericQ,sampleCount_?NumericQ,wavelengthInNanoMeters_?NumericQ]:=With[{index1=Floor[getIndexOfWavelength[wavelengths,wavelengthInNanoMeters]],index2=Ceiling[getIndexOfWavelength[wavelengths,wavelengthInNanoMeters]]},With[{image1=images[[index1]],image2=images[[index2]]},With[{d1=wavelengthInNanoMeters-wavelengths[[index1]],d2=wavelengths[[index2]]-wavelengthInNanoMeters},Table[(1-d1/(d1+d2))image1[[i]]+image2[[i]](d1/(d1+d2)),{i,1,Length[image1]}]]]]


(* ::Input::Initialization:: *)
getImageForWavelengthWithoutInterpolation[rawData_,totalBandCount_?NumericQ,sampleCount_?NumericQ,wavelengthInNanoMeters_?NumericQ]:=With[{index1=Floor[getIndexOfWavelength[wavelengthInNanoMeters]],index2=Ceiling[getIndexOfWavelength[wavelengthInNanoMeters]]},With[{d1=wavelengthInNanoMeters-wavelengths[[index1]],d2=wavelengths[[index2]]-wavelengthInNanoMeters},If[d1<d2,getImageData[rawData,index1-1,totalBandCount,sampleCount],getImageData[rawData,index2-1,totalBandCount,sampleCount]]]]


(* ::Input::Initialization:: *)
getImageForWavelengthWithoutInterpolation[images_List,wavelengths_List,totalBandCount_?NumericQ,sampleCount_?NumericQ,wavelengthInNanoMeters_?NumericQ]:=With[{index1=Floor[getIndexOfWavelength[wavelengths,wavelengthInNanoMeters]],index2=Ceiling[getIndexOfWavelength[wavelengths,wavelengthInNanoMeters]]},With[{d1=wavelengthInNanoMeters-wavelengths[[index1]],d2=wavelengths[[index2]]-wavelengthInNanoMeters},If[d1<d2,images[[index1]],images[[index2]]]]]


(* ::Input::Initialization:: *)
Options[createReferenceSpectrum]={PlotColor->GrayLevel[0],YLabel->"Signal intensity (a.u.)"};


(* ::Input::Initialization:: *)
createReferenceSpectrum[curve_List,minY_,maxY_,tickIntervalY_,opts:OptionsPattern[{createReferenceSpectrum}]]:=With[{minIndex=2,plotColor=OptionValue[PlotColor]},With[{},Show[ListPlot[curve,Joined->True,PlotStyle->{Thick,plotColor},PlotRange->{minY,maxY}],PlotRange->{minY,maxY},Ticks->{LinTicks[0,1000,100,1,MajorTickLength->{0.013,0},MinorTickLength->{0.006,0}],LinTicks[minY,maxY,tickIntervalY,1,MajorTickLength->{0.011,0},MinorTickLength->{0.006,0}]},Epilog-> {Inset[Style[OptionValue[YLabel],Directive[GrayLevel[0],FontFamily->"Arial",12]],Scaled@{-0.15,0.5},Automatic,Automatic,{0,1}],Inset[Style["Wavelength (nm)",Directive[GrayLevel[0],FontFamily->"Arial",12],SingleLetterItalics->True],Scaled@{0.5,-0.16},Automatic,Automatic]},TicksStyle->Directive[GrayLevel[0],FontFamily->"Arial",12],PlotRangeClipping->False,PlotRangePadding->0,ImagePadding->{{50,12},{45,10}},AxesStyle->Directive[AbsoluteThickness[1]],ImageSize->{350,Automatic}]]]


(* ::Input::Initialization:: *)
End[];


(* ::Input::Initialization:: *)
EndPackage[]
