%this script does timeseries extraction on the fMRI data of subjects
%doing three different tasks. This generates a time series by region of
%interest matrix for each subject and each task. Next it highpass filters that
%timebyROI variable to exclude frequencies lower than 0.008 Hz and saves
%the output as the variable hp_filtered. The script then removes time
%points with framewise displacement values over 0.3 and saves that output
%as FDremoved. Lastly, this script creates a connectivity matrix for each
%subject, saves those variables, and creates one cumulative connectivity
%matrix for each task.

subjects = {'sub-RAD14' 'sub-RAD16' 'sub-RAD17' 'sub-RAD18' 'sub-RAD23' 'sub-RAD24' 'sub-RAD28' 'sub-RAD29' 'sub-RAD30' 'sub-RAD32' 'sub-RAD33' 'sub-RAD35' 'sub-RAD37' 'sub-RAD38' 'sub-RAD39' 'sub-RAD40' 'sub-RAD41' 'sub-RAD43' 'sub-RAD44' 'sub-RAD45' 'sub-RAD47' 'sub-RAD49' 'sub-RAD50' 'sub-RAD51' 'sub-RAD54' 'sub-RAD55' 'sub-RAD56' 'sub-RAD57' 'sub-RAD59' 'sub-RAD60' 'sub-RAD61' 'sub-RAD65' 'sub-RAD66' 'sub-RAD67' 'sub-RAD68' 'sub-RAD69' 'sub-RAD70' 'sub-RAD72' 'sub-RAD73' 'sub-RAD75' 'sub-RAD77' 'sub-RAD78' 'sub-RAD79' 'sub-RAD82' 'sub-RAD85' 'sub-RAD86' 'sub-RAD88' 'sub-RAD90' 'sub-RAD93' 'sub-RAD95' 'sub-RAD97' 'sub-RAD99' 'sub-RAD100' 'sub-RAD101' 'sub-RAD103' 'sub-RAD104' 'sub-RAD105' 'sub-RAD107' 'sub-RAD109' 'sub-RAD111' 'sub-RAD112' 'sub-RAD113' 'sub-RAD115' 'sub-RAD116' 'sub-RAD117' 'sub-RAD119' 'sub-RAD120' 'sub-RAD121' 'sub-RAD125' 'sub-RAD127' 'sub-RAD128' 'sub-RAD129' 'sub-RAD130' 'sub-RAD131' 'sub-RAD132' 'sub-RAD135' 'sub-RAD140' 'sub-RAD141' 'sub-RAD142' 'sub-RAD144' 'sub-RAD145' 'sub-RAD149' 'sub-RAD150' 'sub-RAD157' 'sub-RAD158' 'sub-RAD159' 'sub-RAD160' 'sub-RAD161' 'sub-RAD162' 'sub-RAD164' 'sub-RAD165' 'sub-RAD167' 'sub-RAD170' 'sub-RAD171' 'sub-RAD172' 'sub-RAD176' 'sub-RAD177' 'sub-RAD182' 'sub-RAD186' 'sub-RAD187' 'sub-RAD190' 'sub-RAD191' 'sub-RAD192' 'sub-RAD193' 'sub-RAD194' 'sub-RAD196' 'sub-RAD202' 'sub-RAD203' 'sub-RAD205' 'sub-RAD207' 'sub-RAD208' 'sub-RAD210' 'sub-RAD212' 'sub-RAD217' 'sub-RAD218' 'sub-RAD219' 'sub-RAD220' 'sub-RAD221' 'sub-RAD226' 'sub-RAD229' 'sub-RAD231' 'sub-RAD235' 'sub-RAD236' 'sub-RAD238' 'sub-RAD256' 'sub-RAD261' 'sub-RAD270' 'sub-RAD287' 'sub-RAD298' 'sub-RAD299' 'sub-RAD309' 'sub-RAD311' 'sub-RAD312' 'sub-RAD321' 'sub-RAD322' 'sub-RAD328' 'sub-RAD334' 'sub-RAD335' 'sub-RAD337' 'sub-RAD342' 'sub-RAD346' 'sub-RAD347' 'sub-RAD353' 'sub-RAD357' 'sub-RAD359' 'sub-RAD361' 'sub-RAD362' 'sub-RAD364' 'sub-RAD369' 'sub-RAD377' 'sub-RAD379' 'sub-RAD389' 'sub-RAD390' 'sub-RAD392' 'sub-RAD396' 'sub-RAD401' 'sub-RAD402' 'sub-RAD403' 'sub-RAD404'}
task={ 'con' 'cpt' 'gng'}

%timeseries extraction; output = timebyROI
atlas = niftiread('/Users/leonardotozzi/Desktop/Server_Leo/ChelseaStudy/BN_Atlas_246_1mm_corrected.nii');
num_ROIs = max(atlas(:));
for isub = 1:length(subjects)
    
    for itask = 1:length(task)
        
        try
            
            fMRI_data = niftiread(strcat('/Users/leonardotozzi/Desktop/Server_Leo/RAD_Study_preproc/fmriprep/', subjects{isub}, '/func/', subjects{isub},'_',task{itask},'_ICA_AROMA_cut3/denoised_func_data_nonaggr.nii.gz'));
            
        catch
            
            disp(strcat(subjects{isub},'_',task{itask},{' '}, 'not found' ))
            continue
            
        end
        
        if exist(strcat('/Users/leonardotozzi/Desktop/Server_Leo/ChelseaStudy/Connmats/',subjects{isub},'_', task{itask}, '_connmat.csv'))
            
            disp('Matrix already computed')
            
        else
            
            [x,y,z,time] = size(fMRI_data);
            timebyROI = nan(time, num_ROIs);
            for iROI = 1:num_ROIs
                voxels4ROI = find(atlas==iROI);
                for itime = 1:time
                    curr_fMRI = fMRI_data(:,:,:,itime);
                    itime_mean = mean(curr_fMRI(voxels4ROI));
                    timebyROI(itime,iROI) = itime_mean;
                end
            end
            
            fprintf(strcat('Just finished timeseries extraction of ',subjects{isub},' ', task{itask}, '\n \n'))
            
            hp_filtered = highpass(timebyROI, 0.008,0.5);
            fprintf(strcat('Just finished high pass filtering of ',subjects{isub},' ', task{itask}, '\n \n'))
            
            FD_data = dlmread(strcat('/Users/leonardotozzi/Desktop/Server_Leo/RAD_Study_preproc/fmriprep/',subjects{isub},'/func/',subjects{isub},'_task-',task{itask},'_desc-confounds_regressors.tsv'),'\t',4,0);
            FD2remove = FD_data(:,6) >= 0.3;
            hp_filtered(FD2remove,:) = [];
            FDremoved = hp_filtered;
            fprintf(strcat('Just finished framewise displacement removal of ',subjects{isub},' ', task{itask}, '\n \n'))
            
            conn_matrix = corr(FDremoved,'type','Pearson');
            csvwrite(strcat('/Users/leonardotozzi/Desktop/Server_Leo/ChelseaStudy/Connmats/',subjects{isub},'_', task{itask}, '_connmat.csv'),conn_matrix);
            
        end
    end
end


