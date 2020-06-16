% ------------------------------------------------------------------------
%      Analyze non-linear ROIs that were clicked as 'object points'
% ------------------------------------------------------------------------


%% ENTER PARAMETERS AND FILE LOCATION

% file location of object points
points_folder = 'D:\Dropbox (UCL - SWC)\Project_transcriptomics\analysis\PAG_registration\PAG_cells_to_register\processed';

% directory to save results
save_folder = 'D:\Dropbox (UCL - SWC)\Project_transcriptomics\analysis\PAG_registration\PAG_cells_to_register\processed\results';

% directory of reference atlas files
annotation_volume_location = 'D:\PhD\GitHub\allenCCF_philip\annotation_volume_10um_by_index.npy'; % from the allen inst (see readme)
structure_tree_location = 'D:\PhD\GitHub\allenCCF_philip\structure_tree_safe_2017.csv'; % located in github repo
CCF_to_FP_location =  'D:\PhD\GitHub\allenCCF_philip\CCF_to_FP.csv'; % located in github repo
FP_table_location = 'D:\PhD\GitHub\allenCCF_philip\FP_table_Chon_2020.csv'; % located in github repo
chon_images_loc = 'D:\PhD\GitHub\allenCCF_philip\Chon_et_al_Supp_Data_Labels'; % from chon et al (supplementary data 4, https://www.nature.com/articles/s41467-019-13057-w)

% name of the saved object points
object_save_name_suffix = '_PAG_scRNAseq_registered_cells';

% either set to 'all' or a list of indices from the clicked objects in this file, e.g. [2,3]
objects_to_analyze = 'all';

% plane used to view when points were clicked ('coronal' -- most common, 'sagittal', 'transverse')
plane = 'coronal';

% brain figure black or white
black_brain = true;


%% LOAD THE REFERENCE ANNOTATIONS AND PROBE POINTS

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

if ~exist('CCFtoFPtable','var') || ~exist('FPtable','var')
    CCFtoFPtable = loadCCFtoFP(CCF_to_FP_location);
    FPtable = loadFPtable(FP_table_location);
end

% load object points
objectPoints = load(fullfile(points_folder, ['probe_points' object_save_name_suffix]));
% objectPoints.pointList.pointList{1, 3} % This should output the image names of VGAT cells
% objectPoints.pointList.pointList{2, 3} % This should output the image names of vGluT2 cells

% determine which objects to analyze
if strcmp(objects_to_analyze,'all')
    objects = 1:size(objectPoints.pointList.pointList,1);
else
    objects = objects_to_analyze;
end 


%% BRING UP THE RELEVANT DATA FOR EACH PROBE POINTS, FOR FURTHER ANALYSIS

% initialize cell array containing info on each clicked point
if length(objects) > 1
    roi_annotation_CCF = cell(length(objects),1);
    roi_annotation_FP = cell(length(objects),1);
    roi_location_bregma = cell(length(objects),1);
    roi_location_atlas = cell(length(objects),1);
end

% generate needed values
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm

% plot brain grid
% ProbeColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 
% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','purple','orange','red'};

% to have the same VGAT/vGluT2 colors (probe 1 contains VGAT cells, probe 2 contains vGluT2 cells), do the following:
ProbeColors = [1 .502 .502; .059 .6 .698;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0];

fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';

for object_num = objects
    
    selected_object = objects(object_num);
    
    % add the cell type to a new variable (probe 1 contains VGAT cells, probe 2 contains vGluT2 cells)
    if object_num == 1
        curr_cell_type = 'VGAT';
    elseif object_num == 2
        curr_cell_type = 'VGluT2';
    end
        
    % get the object points for the currently analyzed object    
    if strcmp(plane,'coronal')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [3 2 1]);
        curr_slice_id = objectPoints.pointList.pointList{selected_object, 3}';
    elseif strcmp(plane,'sagittal')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [1 2 3]);
        curr_slice_id = objectPoints.pointList.pointList{selected_object, 3}';
    elseif strcmp(plane,'transverse')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [1 3 2]);
        curr_slice_id = objectPoints.pointList.pointList{selected_object, 3}';
    end
    
    % Fix: images were flipped, need to flip back ML coordinates to get the hemispheres correctly
    % CCF coordinates are: 0 to 1320 (AP), 0 to 800 (DV), 0 to 1140 (ML). AP: 0 is more anterior; ML: 0 is more left; DV: 0 is more dorsal
    % A quick fix to swap sides and obtain the cells in the right hemisphere, we can just do ML_coord = 1140 - ML_coord.
    ml_atlas_original = curr_objectPoints(:,3); % save the wrong ML coordinates before fixing them
    curr_objectPoints(:,3) = 1140 - curr_objectPoints(:,3); % flip ML values so hemispheres are correct

    % plot points on the wire frame brain
    figure(fwireframe); hold on
    hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.','linewidth',2, 'color',[ProbeColors(object_num,:) .2],'markers',10);   

    % use the point's position in the atlas to get the AP, DV, and ML coordinates
    ap = -(curr_objectPoints(:,1)-bregma(1))*atlas_resolution;
    dv = (curr_objectPoints(:,2)-bregma(2))*atlas_resolution;
    ml = (curr_objectPoints(:,3)-bregma(3))*atlas_resolution;

    roi_location_curr_bregma = [ap dv ml];
    
    % save the coordinates in atlas space too
    ap_atlas = curr_objectPoints(:,1);
    dv_atlas = curr_objectPoints(:,2);
    ml_atlas = curr_objectPoints(:,3);

    roi_location_curr_atlas = [ap_atlas dv_atlas ml_atlas ml_atlas_original];
    
    % initialize array of region annotations
    roi_annotation_curr = cell(size(curr_objectPoints,1),3);    
    roi_annotation_curr_FP = cell(size(curr_objectPoints,1),3);    
    
    % loop through every point to get ROI locations and region annotations
    for point = 1:size(curr_objectPoints,1)

        % find the annotation, name, and acronym of the current ROI pixel
        ann = av(curr_objectPoints(point,1),curr_objectPoints(point,2),curr_objectPoints(point,3));
        name = st.safe_name{ann};
        acr = st.acronym{ann};

        roi_annotation_curr{point,1} = ann;
        roi_annotation_curr{point,2} = name;
        roi_annotation_curr{point,3} = acr;
        
        % add the cell type
        roi_annotation_curr{point,4} = curr_cell_type;
        
        % add the image name of the current cell
        curr_cell_id = curr_slice_id{point, 1};
        roi_annotation_curr{point,5} = curr_cell_id;
        
        % find the FP annotation, name, and acronym of the current ROI pixel
        [ann_FP, name_FP, acr_FP] = CCF_to_FP(curr_objectPoints(point,1), curr_objectPoints(point,2), curr_objectPoints(point,3), CCFtoFPtable, FPtable, chon_images_loc);

        roi_annotation_curr_FP{point,1} = ann_FP;
        roi_annotation_curr_FP{point,2} = name_FP;
        roi_annotation_curr_FP{point,3} = acr_FP;

    end
    
    % save results in cell array
    if length(objects) > 1
        roi_annotation_CCF{object_num} = roi_annotation_curr;
        roi_annotation_FP{object_num} = roi_annotation_curr_FP;
        roi_location_bregma{object_num} = roi_location_curr_bregma;
        roi_location_atlas{object_num} = roi_location_curr_atlas;
    else
        roi_annotation_CCF = roi_annotation_curr;
        roi_annotation_FP = roi_annotation_curr_FP;
        roi_location_bregma = roi_location_curr_bregma;
        roi_location_atlas = roi_location_curr_atlas;
    end
  
    % display results in a table
    disp(['Clicked points for object ' num2str(selected_object)])
    roi_table = table(roi_annotation_curr(:,5), ...
                      roi_annotation_curr(:,4), ...
                      roi_annotation_curr(:,2), ...
                      roi_annotation_curr(:,3), ...
                      roi_annotation_curr_FP(:,2), ...
                      roi_annotation_curr_FP(:,3), ...
                      roi_location_curr_bregma(:,1), ...
                      roi_location_curr_bregma(:,2), ...
                      roi_location_curr_bregma(:,3), ...
                      roi_location_curr_atlas(:,1), ...
                      roi_location_curr_atlas(:,2), ...
                      roi_location_curr_atlas(:,3), ...
                      roi_location_curr_atlas(:,4), ...
                      roi_annotation_curr(:,1), ...
                      roi_annotation_curr_FP(:,1), ...
                      'VariableNames', {'cell_ID', 'cell_type', ...
                      'CCF_brain_area', 'CCF_acronym', ...
                      'FP_brain_area', 'FP_acronym', ...
                      'AP_location', 'DV_location', 'ML_location', ...
                      'AP_atlas', 'DV_atlas', 'ML_atlas', 'ML_atlas_original', ...
                      'avIndex', 'fpIndex'});
    disp(roi_table)
    
end

% Save combined table with both VGAT and VGluT2 cells:
PAG_cells_table = table([roi_annotation_CCF{1,1}(:,5); roi_annotation_CCF{2,1}(:,5)], ...
                        [roi_annotation_CCF{1,1}(:,4); roi_annotation_CCF{2,1}(:,4)], ...
                        [roi_annotation_CCF{1,1}(:,2); roi_annotation_CCF{2,1}(:,2)], ...
                        [roi_annotation_CCF{1,1}(:,3); roi_annotation_CCF{2,1}(:,3)], ...
                        [roi_annotation_FP{1,1}(:,2); roi_annotation_FP{2,1}(:,2)], ...
                        [roi_annotation_FP{1,1}(:,3); roi_annotation_FP{2,1}(:,3)], ...
                        [roi_location_bregma{1,1}(:,1); roi_location_bregma{2,1}(:,1)], ...
                        [roi_location_bregma{1,1}(:,2); roi_location_bregma{2,1}(:,2)], ...
                        [roi_location_bregma{1,1}(:,3); roi_location_bregma{2,1}(:,3)], ...
                        [roi_location_atlas{1,1}(:,1); roi_location_atlas{2,1}(:,1)], ...
                        [roi_location_atlas{1,1}(:,2); roi_location_atlas{2,1}(:,2)], ...
                        [roi_location_atlas{1,1}(:,3); roi_location_atlas{2,1}(:,3)], ...
                        [roi_location_atlas{1,1}(:,4); roi_location_atlas{2,1}(:,4)], ...
                        [roi_annotation_CCF{1,1}(:,1); roi_annotation_CCF{2,1}(:,1)], ...
                        [roi_annotation_FP{1,1}(:,1); roi_annotation_FP{2,1}(:,1)], ...
                        'VariableNames', {'cell_ID', 'cell_type', ...
                        'CCF_brain_area', 'CCF_acronym', ...
                        'FP_brain_area', 'FP_acronym', ...
                        'AP_location', 'DV_location', 'ML_location', ...
                        'AP_atlas', 'DV_atlas', 'ML_atlas', 'ML_atlas_original', ...
                        'avIndex', 'fpIndex'});
disp(PAG_cells_table)

% Save table and results:
save(fullfile(save_folder, 'PAG_scRNAseq_registered_cells'), 'objectPoints'); % should be the same as ['probe_points' object_save_name_suffix]
save(fullfile(save_folder, 'roi_location_atlas'), 'roi_location_atlas');
save(fullfile(save_folder, 'roi_location_bregma'), 'roi_location_bregma');
save(fullfile(save_folder, 'roi_annotation_CCF'), 'roi_annotation_CCF');
save(fullfile(save_folder, 'roi_annotation_FP'), 'roi_annotation_FP');
writetable(PAG_cells_table, fullfile(save_folder,'PAG_cells_registration_table.csv'),'Delimiter',',','QuoteStrings',true)

% now, use roi_location and roi_annotation for your further analyses
