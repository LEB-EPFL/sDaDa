%-------------------------------------------------------------------------%
%-            THIS IS THE MAIN FUNCTION: 'SIManalysis'                   -%
%-------------------------------------------------------------------------%
%
% This function analyses SIM dataset for doing edge detection..
% format:[StackOut, maskStackOut, objectsOut] = SIManalysis(FileName, PathName, PathNameResults, FileNameZ)
% INPUT:
% - FileName: file name of image stack
% - PathName: path for the image stack
% - PathNameResults: path to save output figures to
% - FileNameZ: optional, filename for second channel (divisome protein)
% This function returns in OUTPUT:
% - imageStack: grayscaled, normalized imageStack (for debugging/testing)
% - cellInfo: a cell cointaining all the info af all the bacteria at each
%             frame
% - num_frames: number of frames
% - divededVarEachFrm: n-by-3 matrix, 
%
% Authors: Aster Vanhecke and Anna Archetti
%--------------------------------------------------------------------------

function  [imageStack,  num_frames, cellInfo, divededVarEachFrm] = SIManalysis(varargin)
narginchk(3,4);
FileName=varargin{1,1};
PathName=varargin{1,2};
PathNameResults=varargin{1,3};
if size(varargin,2) == 3
    DualC=false;
elseif size(varargin,2) == 4
    FileNameZ=varargin{1,4};
    DualC=true;
else
    error('Wrong number of input parameters, there should be 3 or 4') 
end
% Load the image stack and convert from rgb to gray scale
[imageStack, imgInfo, num_frames] = imageLoad(FileName, PathName);
width = imgInfo.Width;
height = imgInfo.Height;
% figure, imshow(imageStack(:,:,1)) % debug: plot first frame loaded image
if DualC
    [imageStackZ, ~, num_framesZ] = imageLoad(FileNameZ, PathName);
    if num_frames ~= num_framesZ
       error('The two channels do not hold the same amount of frames!')
    end
end
% Segmentation of the image
if ~DualC
    [cellInfo, divededVarEachFrm] = imageSegmentation(width,height,num_frames,imageStack,PathNameResults); 
else
    [cellInfo, divededVarEachFrm] = imageSegmentation(width,height,num_frames,imageStack,PathNameResults,imageStackZ);
end
end

%-------------------------------------------------------------------------%
%-                      THESE ARE THE OTHER METHODS                      -%
%-------------------------------------------------------------------------%
 
% Functions for loading and preparing imageStacks 
function [imageStack, imgInfo, num_frames] = imageLoad(FileName, PathName)
        % This function allows to load  a stack of .tif images in matlab as
        % a 3D matlab matrix.
        % If the format is RGB, it will convert to grayscale, if not, it
        % will just load the image stack.

        filesPathAndName = strcat(PathName, FileName);
        
        % Grab file info.
        imgInfo = imfinfo(filesPathAndName);
        num_frames = numel(imgInfo);
        width = imgInfo.Width;
        hight = imgInfo.Height;
        
        % Display info
        disp('Number of frames in the first file:')
        disp(num_frames)
        disp('Image Width and Hight:')
        disp(width)
        disp(hight)
        disp(['Bit depth: ' num2str(imgInfo(1).BitDepth)])
        disp(['Color type: ' imgInfo(1).ColorType])
        
        % preallocation
        imageStack = zeros(imgInfo(1).Height, imgInfo(1).Width, num_frames);
        
        tic
        switch imgInfo(1).ColorType
            case 'truecolor'
               for frIdx = 1:num_frames % if RGB, convert to grayscale during loading
                RGB=imread(filesPathAndName,frIdx);
                imageStack(:,:,frIdx)= rgb2gray(RGB);
               end 
               max=2^(imgInfo(1).BitDepth/3);
            case 'grayscale'
                for frIdx=1:num_frames
                imageStack(:,:,frIdx)=imread(filesPathAndName,frIdx);
                end
                max=2^(imgInfo(1).BitDepth);
            otherwise
                error('Unsupported filetype')
        end
        imageStack=NormalizeISt(imageStack, imgInfo, num_frames,max);
        toc
end
%__________________________________________________________________________

% Functions for normalize the images 
function imageStackNorm = NormalizeISt(imageStack, imgInfo, num_frames,max)
    % Format function imageStackNorm = NormalizeISt(imageStack,num_frames) 
    % Function to normalize imagestacks:
    % For each frame it divides the value of each pixel by the max possible 
    % value for that bitdepth, after substracting the background (mean2).
    % Background substraction is necessary for meaningful FWHM measurements.
    
    %preallocation
    imageStackNorm = zeros(imgInfo(1).Height, imgInfo(1).Width, num_frames);
    
    for frIdx2 = 1:num_frames
        imageStackNorm(:,:,frIdx2) = (imageStack(:,:,frIdx2) - mean2(imageStack(:,:,frIdx2)))./max;
    end
    imageStackNorm(imageStackNorm<0)=0;
    
end

%__________________________________________________________________________

% Segmentation of the images
function [segmentedAllCellInfo, divededVarEachFrm] = imageSegmentation(varargin) 
% subfunction to loop the segmentation/measurement process over all frames
narginchk(5,6);
width=varargin{1,1};
height=varargin{1,2};
num_frames=varargin{1,3};
originalImage=varargin{1,4};
PathNameResults=varargin{1,5};
if size(varargin,2) == 6
    originalImageZ=varargin{1,6};
    DualC=true;
else
    DualC=false;
end
segmentedAllCellInfo{num_frames,1}=[];%preallocation

format long g;
format compact;

% ScreenSize is a four-elementvector: [left, bottom, width, height]:
scrsz = get(0,'ScreenSize');

% Initializations of variables
%     allSegmentedCellInfo = cell(num_frames, 1);
maskStack = zeros(height,width,num_frames);
num_frames
allSegmentedCellInfo = [];

divededVarPriv = [];   
for frIdx = 1:num_frames

    % First segmentation to eliminate the blob that are not cells
    [maskedImageCell, coloredGoodCell] = generalSeg(originalImage, frIdx, maskStack);

    % First segmentation of the first full field of view
    % To choose which cell to follow      
    if frIdx == 1
        getMeshFig=figure('Position',[scrsz(1,3)/2 50 scrsz(1,3)/2 scrsz(1,4)/2]);
        allSegmentedCellInfo = semiAutoSegForSIM(maskedImageCell, PathNameResults,getMeshFig);
    end
    
    % Segmentation check cell by cell
    try
        segmentedAllCellInfoForInput=segmentedAllCellInfo{frIdx-1, 1};
    catch
        segmentedAllCellInfoForInput=[];
    end
    if DualC
        if frIdx == 1
        FtsZFig=figure;
        end
        [segmentedCellInfo, divededVar] = segCheckCellByCell(allSegmentedCellInfo, maskedImageCell, originalImage, frIdx, divededVarPriv, PathNameResults, getMeshFig, segmentedAllCellInfoForInput, originalImageZ(:,:,frIdx),FtsZFig);
    else
        [segmentedCellInfo, divededVar] = segCheckCellByCell(allSegmentedCellInfo, maskedImageCell, originalImage, frIdx, divededVarPriv, PathNameResults,getMeshFig, segmentedAllCellInfoForInput);
    end
    
    segmentedAllCellInfo{frIdx, 1} = segmentedCellInfo;
    
    %% Single cell drift/growth correction
    % Code to move roibox etc. when cell drifts --> use box and contour of
    % last frame, instead of first frame.
    if frIdx > 1
    for cellIdx = 1:size(allSegmentedCellInfo,2)
        if divededVar(cellIdx,3) == 0 % Only try if the cell is not divided.
        try
        % box: to coordinates, add difference between last two boxes(drift)
        allSegmentedCellInfo{1,cellIdx}.box(1) = allSegmentedCellInfo{1,cellIdx}.box(1)+segmentedAllCellInfo{frIdx,1}{1,cellIdx}.box(1)-segmentedAllCellInfo{frIdx-1,1}{1,cellIdx}.box(1); 
        allSegmentedCellInfo{1,cellIdx}.box(2) = allSegmentedCellInfo{1,cellIdx}.box(2)+segmentedAllCellInfo{frIdx,1}{1,cellIdx}.box(2)-segmentedAllCellInfo{frIdx-1,1}{1,cellIdx}.box(2);
        % size becomes the size of last box
        allSegmentedCellInfo{1,cellIdx}.box(3) = segmentedAllCellInfo{frIdx,1}{1,cellIdx}.box(3);
        allSegmentedCellInfo{1,cellIdx}.box(4) = segmentedAllCellInfo{frIdx,1}{1,cellIdx}.box(4);
        % contour: replace old contour (x,y) by new contour, but add
        % relative position of box. (original box (e.g. 528)- relative
        % box (e.g. 41.5).
        allSegmentedCellInfo{1,cellIdx}.contour=[]; % first clear the contour to prevent matrix dimension mismatch
        allSegmentedCellInfo{1,cellIdx}.contour(:,2) = segmentedAllCellInfo{frIdx,1}{1,cellIdx}.contour(:,2)+allSegmentedCellInfo{1,cellIdx}.box(2)-segmentedAllCellInfo{frIdx,1}{1,cellIdx}.box(2); % y
        allSegmentedCellInfo{1,cellIdx}.contour(:,1) = segmentedAllCellInfo{frIdx,1}{1,cellIdx}.contour(:,1)+allSegmentedCellInfo{1,cellIdx}.box(1)-segmentedAllCellInfo{frIdx,1}{1,cellIdx}.box(1); % x

        allSegmentedCellInfo{1,cellIdx}.centroid = segmentedCellInfo{1,cellIdx}.centroid; % centroid is relative to box
        catch
            fprintf('Drift correction failed for cell %d!\n', cellIdx);
        end
        end
    end
    end
    %%
    divededVarPriv = divededVar;
    if frIdx == 1
         divededVarEachFrm = divededVar;
    else 
        divededVarEachFrm(end+1 : end + size(divededVar, 1), :) = divededVar;
    end

end % end frames
    

 end
%__________________________________________________________________________


function [maskedImageCell, coloredGoodCell] = generalSeg(originalImage, frIdx, maskStack)

    format long g;
    format compact;
    
    % Global image threshold using Otsu's method
        % The output is a normalized intensity value that lies in the range [0, 1].
        bg = graythresh(originalImage(:,:,frIdx));
        
        % Convert image to binary image
        % The output image BW replaces all pixels in the input image with 
        % luminance greater than 'bg' with the value 1 (white) and replaces 
        % all other pixels with the value 0 (black).
        maskStack(:,:,frIdx)=im2bw(originalImage(:,:,frIdx),bg);

        % Do a "hole fill" to get rid of any background pixels or "holes" inside the blobs.
        maskStack(:,:,frIdx) = imfill(maskStack(:,:,frIdx), 'holes');

        % Identify individual blobs by seeing which pixels are connected to each other.
        % Each group of connected pixels will be given a label, a number, to identify it and distinguish it from the other blobs.
        % Do connected components labeling with either bwlabel() or bwconncomp().
        labeledImage = bwlabel(maskStack(:,:,frIdx), 8);     % Label each blob so we can make measurements of it
        % labeledImage is an integer-valued image where all pixels in the blobs have values of 1, or 2, or 3, or ... etc.
        
        % Get all the cells properties.  Can only pass in originalImage in version R2008a and later.
        cellMeasurements = regionprops(labeledImage, originalImage(:,:,frIdx), 'all');
        
        % Select certain cells based using the ismember() function.
        allCellAreas = [cellMeasurements.Area];
        
        % Get a list of the cells that meet our criteria and we need to keep.
        % These will be logical indices - lists of true or false depending on 
        % whether the feature meets the criteria or not.
        
        %allowableIntensityIndexes = (allCellIntensities > 150) & (allCellIntensities < 220);
        allowableAreaIndexes = allCellAreas > 600; % Take the big objects.
        
        % Now let's get actual indexes, rather than logical indexes, of the features that meet the criteria.
        %keeperIndexes = find(allowableIntensityIndexes & allowableAreaIndexes);
        keeperIndexes = find(allowableAreaIndexes);
        
        % Extract only those cells that meet our criteria, and
        % eliminate those cells that don't meet our criteria.
        % Note how we use ismember() to do this.  
        % Result will be an image - the same as labeledImage but with only the cells
        % listed in keeperIndexes in it.
        keeperCellsImage = ismember(labeledImage, keeperIndexes);
        
        % Re-label with only the keeper cells kept.
        labeledCellsImage = bwlabel(keeperCellsImage, 8);     % Label each cell so we can make measurements of it
        % Now we're done.  We have a labeled image of cells that meet our specified criteria.
        
        % Let's assign each blob a different color to visually show the user the distinct blobs.
        coloredGoodCell = label2rgb (labeledCellsImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
        
        % Now use the keeper cells as a mask on the original image.
        % This will let us display the original image in the regions of the keeper cells.
        maskedImageCell = originalImage(:,:,frIdx); % Simply a copy at first.
        maskedImageCell(~keeperCellsImage) = 0;  % Set all non-keeper pixels to zero.
end
