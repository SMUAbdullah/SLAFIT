% Written: 02-March, 2020
% Author: M Saqib Sohail

% update 06-May-2020
%         Linear interpolation option added
% update 3-Jun-2020
%         streamlined inputs
% update 16-Jun-2020
%         option to chose reference sequence by FLAG_UserProvidedRefSeq
%         FLAG_UserProvidedRefSeq = true, user supplies refSeq as an ACGT
%         sequence as the 4rth input parameter to this function
%         FLAG_UserProvidedRefSeq = false, ref seq is the consensus seq of
%         time point 1


% function that implements MPL (binary approx) no epistasis, with
% bootstrap samples.

% 1. Load data
% 2. Calculate frequencies from MSA
% 3. Calculate terms required for MPL estimate
% 4. Set file names to save analysis results
% 5. Calculating MPL estimates
% 6. Save estimates, other data to file
% function analysisStep1_v2(fileNameContainingDirPath, priorConstSC, FLAG_stratonovich, FLAG_MarkAccessibility, FLAG_SaveFile, FLAG_SaveIntCovMtx, FLAG_useFreqEntry, FLAG_troubleShoot, FLAG_linearInt);
function analysis_MPL_10_1(fileNameContainingDirPath,pathdirw0,pathdir_est,priorConst,FLAG_vector,refernceSequence,thisSet,patID,thisProt,noisethresh,varargin)

if(length(priorConst) == 1)
	priorConstSC = priorConst(1);
	priorConstEpi = 0;
elseif(length(priorConst) == 2)
	priorConstSC = priorConst(1);
	priorConstEpi = priorConst(2);
	disp('Warning: Recombination probability not provided. Setting it to 0...press anykey to continue...')
	pause
	recombProb = 0;
elseif(length(priorConst) == 3)
	priorConstSC = priorConst(1);
	priorConstEpi = priorConst(2);
	recombProb = priorConst(3);
else
	disp('Wraning: priorConst should either contain 1 or 2 reg values only. using first two values as \gamma_{sc} and \gamma_{epi}.')
	priorConstSC = priorConst(1);
	priorConstEpi = priorConst(2);
end

if(length(FLAG_vector) == 7)
	FLAG_stratonovich = FLAG_vector(1);
	FLAG_MarkAccessibility = FLAG_vector(2);
	FLAG_UserProvidedRefSeq = FLAG_vector(3);
	FLAG_SaveIntCovMtx = FLAG_vector(4);
	FLAG_useFreqEntry = FLAG_vector(5);
	FLAG_troubleShoot = FLAG_vector(6);
	FLAG_linearInt = FLAG_vector(7);
	FLAG_Epi = false;
elseif(length(FLAG_vector) == 8)
	FLAG_stratonovich = FLAG_vector(1);
	FLAG_MarkAccessibility = FLAG_vector(2);
	FLAG_UserProvidedRefSeq = FLAG_vector(3);
	FLAG_SaveIntCovMtx = FLAG_vector(4);
	FLAG_useFreqEntry = FLAG_vector(5);
	FLAG_troubleShoot = FLAG_vector(6);
	FLAG_linearInt = FLAG_vector(7);
	FLAG_Epi = FLAG_vector(8);
end
userRefSequence = refernceSequence; % comment after testing
FLAG_SaveFile = true;
% --------------------------  initializations -----------------------------
if(ispc)
	chosenSlash = '\';
elseif(isunix)
	chosenSlash = '/';
else
	disp('Error: system is not unix and not PC...')
	pause
end
%----------------------------------------------------------------------
% 1.1 load dir names
[~, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
runAnalysisCode = true;
if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash], 'dir') == 0)
	disp('Can not find Analysis_Misc folder. Run preprocessing Step 1.')
	disp(' Exiting without running Analysis code.')
	runAnalysisCode = false;
end

if(runAnalysisCode)
	disp('Running analysis Step 1 - ')
	disp('-----------------------------------------------------------')
	fileNamesThisPat_inAnalysisDir = getFolderContent(dirNameAnalysis, 'files');
		
	%----------------------------------------------------------------------
	% 1.2 load the .txt file name that contains list of HeaderUpdated fasta files
	if(isempty(fileNamesThisPat_inAnalysisDir))
		disp('Error: Analysis directory is empty. Run preprocess Step 1 before running Step 2')
		pause
	else
		extensionCell{1} = 'txt';
		fileNamesListDataThisPat = findFileNamesWithGivenExtension(dirNameAnalysis, extensionCell);
		numFilesWithTXTExt = length(fileNamesListDataThisPat);
		% find the filename that has '...HU.txt' This is the
		FLAG_fileNameHAsHUTXT = false(1,numFilesWithTXTExt);
		for e = 1:numFilesWithTXTExt
			thisFile = fileNamesListDataThisPat{e};
			indOfHUTXT = strfind(thisFile, 'HU.txt');
			if(~isempty(indOfHUTXT))
				FLAG_fileNameHAsHUTXT(e) = true;
			end
		end
		if(sum(FLAG_fileNameHAsHUTXT) > 1)
			disp('Error: More than 1 files in the Analysis folder has the string HU.txt in it')
			pause
		elseif(sum(FLAG_fileNameHAsHUTXT) == 0)
			disp('Error: No .txt file has string HU.txt in it. Run preprocessing step 1')
			pause
		else
			fileNameContainingHUFiles = fileNamesListDataThisPat{find(FLAG_fileNameHAsHUTXT)};
		end
	end	
	fileNameHUFilesCell = loadDataFileNames([dirNameAnalysis fileNameContainingHUFiles]);
	numFileNameHUFilesCell = length(fileNameHUFilesCell);	
	%----------------------------------------------------------------------
	% 1.3 get timePointVecFromFileName and bsampleVecFromFileName from the
	%     filenames in fileNameHUFilesCell
	
	timePointVecFromFileName = -1*ones(1, numFileNameHUFilesCell);
	bsampleVecFromFileName = -1*ones(1, numFileNameHUFilesCell);
	for f = 1:numFileNameHUFilesCell
		thisFileNameHUFiles = fileNameHUFilesCell{f};
		indOfDash = strfind(thisFileNameHUFiles, '_');		
		% find time point from file name (time has to be the last _t)
		indOfDashT = strfind(thisFileNameHUFiles, '_t');
		indOfDashT = indOfDashT(end);
		temp1 = (indOfDash - indOfDashT > 0); % (logical)
		allIndOfDashGreatherThanDashT = indOfDash(temp1);
		thisTimePoint = str2double(thisFileNameHUFiles(indOfDashT+2:allIndOfDashGreatherThanDashT(1)-1));
		timePointVecFromFileName(f) = thisTimePoint;
		
		% find bsample from filename
		indOfbSample = strfind(thisFileNameHUFiles, '_bsample');
		if(isempty(indOfbSample)) % filenames do not contain the string bsample---no bootstrapsampling done. 1 data file per time point
			bsampleVecFromFileName(f) = 1;
		else
			temp2 = (indOfDash - indOfbSample) > 0;
			allIndOfDashGreatherThanbSample = indOfDash(temp2);
			str1 = thisFileNameHUFiles(indOfbSample+8:(allIndOfDashGreatherThanbSample(1)-1));
			indOfOF = strfind(str1, 'of');
			bsampleVecFromFileName(f) = str2double(str1(1:indOfOF-1));
		end
		
	end    
	%----------------------------------------------------------------------
	% 2. Calculate frequencies from MSA
	%----------------------------------------------------------------------    
	uniqueTimePoints = unique(timePointVecFromFileName);
	numUniqueTimePoints = length(uniqueTimePoints);
	
	timeStep = uniqueTimePoints(2:end) - uniqueTimePoints(1:end-1);
	
	fprintf('Calculating one and two point probabilities...')
	tp = 1;
	thisTP = uniqueTimePoints(tp);
	indOfFileToLoad = find((timePointVecFromFileName == thisTP));
	thisFileNameHUFiles = fileNameHUFilesCell{indOfFileToLoad};    
	% load the aligned haplotypes
	[Header, seqNT_All] = fastaread([dirNameAnalysis thisFileNameHUFiles]);
	% remove 1st seq as that is ref seq (used for alignment purpose only)
	Header = Header(2:end);
	seqNT_All = seqNT_All(2:end);    
	%------------------------------------------------------------------
	% check if header is consistent with MPL format and extract
	% info from header    
	% header consistency check
	tempHeader = Header{1};
	if(contains(tempHeader, 'freq: ') == 1 )
		%disp('pass.')
	else
		fprintf(['Header format of ' thisFileNameHUFiles])
		disp(' not in proper format. Cannot run MPL analysis code.')
		disp(' ')
		disp('Make sure the headers in FASTA files are in proper format for MPL code. See readme.txt')
	end    
	%------------------------------------------------------------------
	% get information from header
	%[thisSeqTimeVec, thisSeqFreqVec, thisSeqNumReadsVec] = getDayFreqNumReedsFromHeader(Header);
	[thisSeqTimeVec, thisSeqFreqVecTemp, thisSeqNumReadsVec] = getInfoFromHeader(Header);    
	timePointVec = unique(thisSeqTimeVec);
	numTimePoints = length(timePointVec);    
	%------------------------------------------------------------------
	% use numReads to find frequency of each seq in the MSA
	%if(sum(thisSeqNumReadsVec == -1) > 0)
	if(FLAG_useFreqEntry == true)
		% numReeds entry not present, use freq entry
		thisSeqFreqVec = thisSeqFreqVecTemp;
	else
		thisSeqFreqVec = [];
		for i = 1:numTimePoints
			thisTimePoint = timePointVec(i);
			readsVecThisTimePoint = (thisSeqNumReadsVec(thisSeqTimeVec == thisTimePoint));
			totalReadsThisTimePoint = sum(readsVecThisTimePoint);
			seqFreqThisTimePointTemp = readsVecThisTimePoint/totalReadsThisTimePoint;
			temp1 = seqFreqThisTimePointTemp/sum(seqFreqThisTimePointTemp);
			temp2 = round(temp1*10000)/10000; % four digit precision
			diff = sum(temp2) - 1;            
			% remove diff from largest entry of temp2
			[val, ind] = max(temp2);
			temp2(ind) = val - diff;
			if(abs(sum(temp2) - 1) > 1e-4)
				disp('Error: check, freq does not sum to 1')
				disp(num2str(sum(temp2)))
				pause
			end
			seqFreqThisTimePoint = seqFreqThisTimePointTemp;
			thisSeqFreqVec = [thisSeqFreqVec seqFreqThisTimePoint];
		end
	end    
	% make compressed MSA
	msaNT = repmat(' ', size(seqNT_All, 2), length(seqNT_All{1}));
	for k = 1:size(seqNT_All, 2)
		msaNT(k,:) = seqNT_All{k};
	end    
	msaNTInt = nt2int(msaNT);
	numSitesNT = size(msaNT,2);
	numUniqueSeqs = size(msaNT,1);    
	freqCountNT = ntFreqCountCompMSA(msaNTInt, thisSeqFreqVec);
	[~, consensusSeqTP1] = max(freqCountNT);    
	%------------------------------------------------------------------
	% 2.1 make binary MSA    
	%  possible states at each site: [0 1], consensus founder:0, mut:1
	%  each site represented by 1
	%disp(' ')
	%fprintf('Converting NT MSA to binary...')    
	numSitesNTExteded = numSitesNT;%numSitesNT*numNTAtEachSite;
	msaNT_Bin = -1*ones(numUniqueSeqs, numSitesNTExteded);    
	if(FLAG_UserProvidedRefSeq == true)
		referenceSequenceForConversion = nt2int(userRefSequence);
	elseif(FLAG_UserProvidedRefSeq == false)
		referenceSequenceForConversion = consensusSeqTP1;
	end
	for j = 1:numSitesNT
		msaNT_Bin(:,j) = referenceSequenceForConversion(j) ~= msaNTInt(:,j);
	end    
	%------------------------------------------------------------------
	% 2.2 Calculate 1,2 point probabilities, mutational flux    
	% Calculate 1-point freqs
	q = thisSeqFreqVec*msaNT_Bin;
	q(q > 1) = 1;
	q(q < 0) = 0;
	indOfNonZeroQs = find(q);
	indOfNonZeroQsLen = length(indOfNonZeroQs);
	q_sp = sparse(q);
	q_All = zeros(numUniqueTimePoints, numSitesNT);
	q_All(tp,:) = q_All(tp,:) + q;    
	%------------------------------------------------------------------
	% 2.2.1 Calculate 2-point freqs
	% q11_vec contains all non-zero 2 points in vector format
	% q11_sp is the sparse matrix representation of q11    
	weighted_msa = thisSeqFreqVec .* msaNT_Bin(:,indOfNonZeroQs)';  % Weight each row
	q11_mat = weighted_msa * msaNT_Bin(:, indOfNonZeroQs);  % Matrix multiplication    
	% Clip values to [0,1] range
	q11_mat(q11_mat < 0) = 0;
	q11_mat(q11_mat > 1) = 1;	
	% Check for out-of-bounds values (optional, for debugging)
	if any(q11_mat < 0, 'all')
		disp('Calculated q11 is less than 0...');
		disp(q11_mat(q11_mat < 0));
	end    
	if any(q11_mat > 1, 'all')
		disp('Calculated q11 is greater than 1...');
		disp(q11_mat(q11_mat > 1) - 1);
	end    
	% Flatten the matrix row-wise to match original output
	q11_vec = reshape(q11_mat', 1, []);
	firstInd = repmat(indOfNonZeroQs, 1, indOfNonZeroQsLen);
	secondInd = reshape(repmat(indOfNonZeroQs, indOfNonZeroQsLen,1),1, indOfNonZeroQsLen*indOfNonZeroQsLen);
	q11_sp = sparse(firstInd, secondInd, q11_vec, numSitesNTExteded, numSitesNTExteded);        
	%disp('done.')
	q_sp_ext = cell(numUniqueTimePoints,1);
	q11_sp_ext = cell(numUniqueTimePoints,1);    
	q_sp_ext{1} = q_sp';
	q11_sp_ext{1} = q11_sp;
	parfor tp = 2:numUniqueTimePoints
		thisTP = uniqueTimePoints(tp);
		% load MSA, calculate 1,2 pt frequencies from it
		indOfFileToLoad = find((timePointVecFromFileName == thisTP));
		thisFileNameHUFiles = fileNameHUFilesCell{indOfFileToLoad};        
		% load the aligned haplotypes
		[Header, seqNT_All] = fastaread([dirNameAnalysis thisFileNameHUFiles]);
		% remove 1st seq as that is ref seq (used for alignment purpose only)
		Header = Header(2:end);
		seqNT_All = seqNT_All(2:end);        
		%------------------------------------------------------------------
		% check if header is consistent with MPL format and extract
		% info from header        
		% header consistency check
		tempHeader = Header{1};
		if(contains(tempHeader, 'freq: ') == 1 )
			%disp('pass.')
		else
			fprintf(['Header format of ' thisFileNameHUFiles])
			disp(' not in proper format. Cannot run MPL analysis code.')
			disp(' ')
			disp('Make sure the headers in FASTA files are in proper format for MPL code. See readme.txt')
		end        
		%------------------------------------------------------------------
		% get information from header
		%[thisSeqTimeVec, thisSeqFreqVec, thisSeqNumReadsVec] = getDayFreqNumReedsFromHeader(Header);
		[thisSeqTimeVec, thisSeqFreqVecTemp, thisSeqNumReadsVec] = getInfoFromHeader(Header);        
		timePointVec = unique(thisSeqTimeVec);
		numTimePoints = length(timePointVec);        
		%------------------------------------------------------------------
		% use numReads to find frequency of each seq in the MSA
		%if(sum(thisSeqNumReadsVec == -1) > 0)
		if(FLAG_useFreqEntry == true)
			% numReeds entry not present, use freq entry
			thisSeqFreqVec = thisSeqFreqVecTemp;
		else
			thisSeqFreqVec = [];
			for i = 1:numTimePoints
				thisTimePoint = timePointVec(i);
				readsVecThisTimePoint = (thisSeqNumReadsVec(thisSeqTimeVec == thisTimePoint));
				totalReadsThisTimePoint = sum(readsVecThisTimePoint);
				seqFreqThisTimePointTemp = readsVecThisTimePoint/totalReadsThisTimePoint;
				temp1 = seqFreqThisTimePointTemp/sum(seqFreqThisTimePointTemp);
				temp2 = round(temp1*10000)/10000; % four digitl precision
				diff = sum(temp2) - 1;                
				% remove diff from largest entry of temp2
				[val, ind] = max(temp2);
				temp2(ind) = val - diff;
				if(abs(sum(temp2) - 1) > 1e-4)
					disp('Error: check, freq does not sum to 1')
					disp(num2str(sum(temp2)))
					pause
				end
				seqFreqThisTimePoint = seqFreqThisTimePointTemp;
				thisSeqFreqVec = [thisSeqFreqVec seqFreqThisTimePoint];
			end
		end        
		% make compressed MSA
		msaNT = repmat(' ', size(seqNT_All, 2), length(seqNT_All{1}));
		for k = 1:size(seqNT_All, 2)
			msaNT(k,:) = seqNT_All{k};
		end        
		msaNTInt = nt2int(msaNT);
		numSitesNT = size(msaNT,2);
		numUniqueSeqs = size(msaNT,1);        
		%------------------------------------------------------------------
		% 2.1 make binary MSA        
		%  possible states at each site: [0 1], consensus founder:0, mut:1
		%  each site represented by 1    
		referenceSequenceForConversion_internal = referenceSequenceForConversion;
		msaNT_Bin = -1*ones(numUniqueSeqs, numSitesNTExteded);
		for j = 1:numSitesNT
			msaNT_Bin(:,j) = referenceSequenceForConversion_internal(j) ~= msaNTInt(:,j);
		end        
		%------------------------------------------------------------------
		% 2.2 Calculate 1,2 point probabilities, mutational flux        
		% Calculate 1-point freqs
		q = thisSeqFreqVec*msaNT_Bin;
		q(q > 1) = 1;
		q(q < 0) = 0;
		indOfNonZeroQs = find(q);
		indOfNonZeroQsLen = length(indOfNonZeroQs);
		q_sp_ext{tp} = sparse(q');
		q_All(tp,:) = q_All(tp,:) + q;
		%------------------------------------------------------------------
		% 2.2.1 Calculate 2-point freqs
		% q11_vec contains all non-zero 2 points in vector format
		% q11_sp is the sparse matrix representation of q11
		weighted_msa = thisSeqFreqVec .* msaNT_Bin(:,indOfNonZeroQs)';  % Weight each row
		q11_mat = weighted_msa * msaNT_Bin(:, indOfNonZeroQs);  % Matrix multiplication        
		% Clip values to [0,1] range
		q11_mat(q11_mat < 0) = 0;
		q11_mat(q11_mat > 1) = 1;        
		% Check for out-of-bounds values (optional, for debugging)
		if any(q11_mat < 0, 'all')
			disp('Calculated q11 is less than 0...');
			disp(q11_mat(q11_mat < 0));
		end        
		if any(q11_mat > 1, 'all')
			disp('Calculated q11 is greater than 1...');
			disp(q11_mat(q11_mat > 1) - 1);
		end        
		% Flatten the matrix row-wise to match original output
		q11_vec = reshape(q11_mat', 1, []);		
		firstInd = repmat(indOfNonZeroQs, 1, indOfNonZeroQsLen);
		secondInd = reshape(repmat(indOfNonZeroQs, indOfNonZeroQsLen,1),1, indOfNonZeroQsLen*indOfNonZeroQsLen);
		q11_sp_ext{tp} = sparse(firstInd, secondInd, q11_vec, numSitesNTExteded, numSitesNTExteded);
	end   
	q_reads = cell2mat(q_sp_ext);
	vld_ind=~((abs(max(q_reads)-min(q_reads))<=noisethresh)&((max(q_reads)==1)|(min(q_reads)==0)));
	[numtp,loci]=size(q_reads);	
	numTP = numUniqueTimePoints;
	guessOfNonZeroEntries = ceil((numSitesNTExteded^2)*0.07);
	timeSteps = [0 timeStep(1:end-1)]; % Adjust timeStep vector for vectorization	
	% Initialize
	useGPU = gpuDeviceCount > 0;
	q_0 = q_sp_ext{1};  % CPU or GPU, depending on need
	q_T = q_sp_ext{end};
	vEst1 = sparse(numSitesNTExteded, 1);  % Column vector
	vEst2 = sparse(numSitesNTExteded, 1);  % Column vector
	intCovMtx = spalloc(numSitesNTExteded, numSitesNTExteded, guessOfNonZeroEntries);
	% Process first time point (tp=1)
	% last_q_sp = q_sp_ext{1};
	% last_q11_sp = q11_sp_ext{1};
	% Main loop (tp > 1)
	% Preprocess data for GPU (outside loop)
	if useGPU
		q_sp_ext = cellfun(@(x) gpuArray(full(x)), q_sp_ext, 'UniformOutput', false);
		q11_sp_ext = cellfun(@(x) gpuArray(full(x)), q11_sp_ext, 'UniformOutput', false);
		last_q_sp = gpuArray(full(q_sp_ext{1}));
		last_q11_sp = gpuArray(full(q11_sp_ext{1}));
	else
		last_q_sp = q_sp_ext{1};
		last_q11_sp = q11_sp_ext{1};
	end
	for tp = 2:numUniqueTimePoints
		% Transfer data (already dense if GPU)
		this_q_sp = q_sp_ext{tp};
		this_q11_sp = q11_sp_ext{tp};
		timeStep_gpu = max(timeStep(tp-1), eps);		
		% Compute intermediates
		q_sp_mid = (last_q_sp + this_q_sp) / 2;		
		% Covariance (GPU/CPU-agnostic)
		term1 = (last_q11_sp + this_q11_sp) / 2;
		term2 = (last_q_sp * last_q_sp')/3 + (this_q_sp * this_q_sp')/3 + ...
				(last_q_sp * this_q_sp')/6 + (this_q_sp * last_q_sp')/6;						
		thisCovMtx = term1 - term2;				
		% Thresholding
		thisCovMtx(abs(thisCovMtx) < 1e-10) = 0;
		thisCovMtx(abs(thisCovMtx) > 1) = sign(thisCovMtx(abs(thisCovMtx) > 1));		
		% Accumulate results
		if useGPU
			intCovMtx = intCovMtx + gather(thisCovMtx * timeStep_gpu);
			vEst1 = vEst1 + gather(1 - q_sp_mid) * timeStep_gpu;  % Dense
			vEst2 = vEst2 + gather(-q_sp_mid) * timeStep_gpu;
		else
			intCovMtx = intCovMtx + thisCovMtx * timeStep_gpu;
			vEst1 = vEst1 + sparse(1 - q_sp_mid) * timeStep_gpu;
			vEst2 = vEst2 + sparse(-q_sp_mid) * timeStep_gpu;
		end		
		% Update for next iteration
		last_q_sp = this_q_sp;
		last_q11_sp = this_q11_sp;
	end 
	%----------------------------------------------------------------------
	% 4. Set file names to save analysis results
	%----------------------------------------------------------------------
	% check accessibility of estimates and save to file
	if(FLAG_stratonovich == true)
		convention = 'Stratonovich';
	elseif(FLAG_stratonovich == false && FLAG_linearInt == false)
		convention = 'Ito';
	elseif(FLAG_stratonovich == false && FLAG_linearInt == true)
		convention = 'LinearInter';
	end	
	thisFileNameHUFiles = fileNameHUFilesCell{1};
	indOfDash = strfind(thisFileNameHUFiles, '_');	
	fileNameSelEst = 's_LD.txt';
	fileNameSelEstSL = 's_LD_independent.txt';
	fileNameIntCovMtx = 'aggregated_LD.txt';
	fileNameAllTrajs = 'SNV_traj.txt';
	fileNameAllTrajsR = 'SNV_traj_with_time.txt';
	fileNameAccessEst = 'Accessibility.txt';	
	%----------------------------------------------------------------------
	% 4.1 Calculate accessibitily
	% right now, can not handle too large matricies
	if(FLAG_MarkAccessibility && (size(intCovMtx, 1) < 5000))
		disp(' ')
		disp('Calculate accessibility...NOTE: Accessibility code works for additive model with < 5000 poly sites')		
		[~, wellCondSites, illCondSites, cluster] = findIndColsOfHmat(intCovMtx);		
		wellCondSitesLen = length(wellCondSites);
		illCondSitesLen = length(illCondSites);
		numClusters = length(cluster);		
		if(isempty(wellCondSitesLen))
			wellCondSitesLen = 0;
		end		
		if(isempty(illCondSitesLen))
			illCondSitesLen = 0;
		end		
		if(~isempty(numClusters))
			clusterLengths = zeros(1, numClusters);
			for cc = 1:numClusters
				clusterLengths(cc) = length(cluster{cc});
			end
			numColAccDataSave = 2 + numClusters;
		else
			numClusters = 0;
			clusterLengths = 0;
			numColAccDataSave = 2;
		end		
		maxEntriesTemp = max([wellCondSitesLen, illCondSitesLen, clusterLengths]);		
		accDataMtxToWrite = zeros(maxEntriesTemp, numColAccDataSave);
		accDataMtxToWrite(:,1) = [wellCondSites zeros(1, maxEntriesTemp-wellCondSitesLen)]';
		accDataMtxToWrite(:,2) = [illCondSites zeros(1, maxEntriesTemp-illCondSitesLen)]';
		for cc = 1:numClusters
			accDataMtxToWrite(:,2+cc) = [cluster{cc} zeros(1, maxEntriesTemp-clusterLengths(cc))]';
		end		
		if(exist([dirNameAnalysis 'Estimates' chosenSlash], 'dir') == 0)
			mkdir([dirNameAnalysis 'Estimates' chosenSlash])
		end
		dlmwrite([dirNameAnalysis 'Estimates' chosenSlash  fileNameAccessEst], accDataMtxToWrite);
	end
	%----------------------------------------------------------------------
	% 5. Calculating MPL estimates
	%----------------------------------------------------------------------
	%tic;
	% load mutation vectors
	fileNameMutVec = [thisFileNameHUFiles(1:indOfDash(2)-1) '_mutVecs.txt'];	
	mutationVectors = dlmread([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVec]);
	mutVecWT2Mut = mutationVectors(1,:)';
	mutVecMut2WT = mutationVectors(2,:)';	
	% MPL Estimate in sparse matrix format	
	intCovMtx_sub=intCovMtx(vld_ind,vld_ind);	   
	vEstMu = vEst1.*mutVecWT2Mut + vEst2.*mutVecMut2WT;
	regMtx = sparse(1:numSitesNTExteded, 1:numSitesNTExteded, priorConstSC*ones(numSitesNTExteded,1), numSitesNTExteded, numSitesNTExteded);
	regMtx_sub=regMtx(vld_ind,vld_ind);
	numer = (q_T - q_0 - vEstMu);
	numer_sub=numer(vld_ind);
	denom = (intCovMtx_sub + regMtx_sub);
	selEst_sub = denom\numer_sub;
	denomSL = diag(diag(intCovMtx_sub + regMtx_sub));
	selEstSL_sub = denomSL\numer_sub;	
	selEst=zeros(length(intCovMtx),1);
	selEstSL=zeros(length(intCovMtx),1);	
	selEst(vld_ind)=selEst_sub;
	selEstSL(vld_ind)=selEstSL_sub;	
	if(FLAG_troubleShoot == 1)
		qT_sub=q_T(vld_ind);
		q0_sub=q_0(vld_ind);
		selEstNoMu = denom\(qT_sub - q0_sub);
		selEstSLNoMu = denomSL\(qT_sub - q0_sub);
	end
	%-----------------------------------------------------------------
	% 6. Save estimates, other data to file
	%-----------------------------------------------------------------
	if(FLAG_SaveFile)		
		if(exist([pathdir_est chosenSlash], 'dir') == 0)
			mkdir([pathdir_est chosenSlash])
		end		
		if(exist([pathdir_est chosenSlash fileNameSelEst], 'file') == 2)
			delete([pathdir_est chosenSlash fileNameSelEst])
		end
		if(exist([pathdir_est chosenSlash fileNameSelEstSL], 'file') == 2)
			delete([pathdir_est chosenSlash fileNameSelEstSL])
		end
		if(FLAG_SaveIntCovMtx == 1)
			filename = [pathdir_est chosenSlash fileNameIntCovMtx];
			fid = fopen([filename '.bin'], 'wb');
			fwrite(fid, size(full(intCovMtx)), 'int32');
			fwrite(fid, full(intCovMtx), 'double');
			fclose(fid);	
		end
				
		filename = [pathdir_est chosenSlash fileNameSelEst];
		fid = fopen(filename, 'w');
		fprintf(fid, '%.8f\n', full(selEst));
		fclose(fid);
	
		filename = [pathdir_est chosenSlash fileNameSelEstSL];
		fid = fopen(filename, 'w');
		fprintf(fid, '%.8f\n', full(selEstSL));
		fclose(fid);
	
		filename = [pathdir_est chosenSlash fileNameAllTrajs];
		fid = fopen([filename '.bin'], 'wb');
		fwrite(fid, size(q_All), 'int32');
		fwrite(fid, q_All, 'double');
		fclose(fid);
	
		filename = [pathdir_est chosenSlash fileNameAllTrajsR];
		fid = fopen([filename '.bin'], 'wb');
		fwrite(fid, size([uniqueTimePoints' q_All]), 'int32');
		fwrite(fid, [uniqueTimePoints' q_All], 'double');
		fclose(fid);
	end
end