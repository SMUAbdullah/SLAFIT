function MPL_7(main_dir,data_dir,this_set,protein,output_dir,genome_start,genome_end,gamma,thresh)

thisSet=char(this_set);
thisProt=char(protein);
pathdirw0=char(output_dir);
maindir=char(main_dir);
datadir=char(data_dir);
thisGenomicSegStartInd=str2double(genome_start);
thisGenomicSegStopInd=str2double(genome_end);
reg=str2double(gamma);
noisethresh=str2double(thresh);

PreprocessingStep_0_8;

PreProcessingStep_1_9;

AnalysisMPL_10

end
