clear;
fprintf('batch_copy ver = 1.0.0\n');
rawDir=fullfile('Y:\ZY_multiscale\data\20191112_sbr_lf_wf\sbr_lf_1112');



%%
rawDir2 = strrep(rawDir,'/','\\');
rawDirSp = strsplit(rawDir2,'\\');
outDir = strcat(rawDirSp{end},'_output');
mkdir(outDir);
dirOutput=dir(rawDir);
source={dirOutput.name}';
source = source(3:end);

for sourceIdx = 1:size(source,1)
    try
        ticer = tic;
        dirName = source{sourceIdx};
        if(~exist(strcat(rawDir,'/',dirName,'/MMStack_Default.ome.tif'),'file'))
            continue;
        end
        sz = 0;
        %src = strcat(rawDir,'/',dirName,'/MMStack_Default.ome.tif');
        src = sprintf('%s/%s/MMStack_Default.ome.tif', rawDir, dirName);
        dst = sprintf('%s/%s_1x1_xx.xms_XXXX_Hardware_LaserCount1_xxxxxxxxxxxx.0.tiff', outDir, dirName);
        %dst = strcat(outDir, '/', dirName,'_1x1_xx.xms_XXXX_Hardware_LaserCount1_xxxxxxxxxxxx.0.tiff');
        FileInfo = dir(src);
        sz = sz + FileInfo.bytes;
        copyfile(src, dst);
        for i = 1:999
            src = sprintf('%s/%s/MMStack_Default_%d.ome.tif', rawDir, dirName, i);
            if(~exist(src, 'file'))
                break;
            end
            dst = sprintf('%s/%s_1x1_xx.xms_XXXX_Hardware_LaserCount1_xxxxxxxxxxxx.%d.tiff', outDir, dirName, i);
             FileInfo = dir(src);
            sz = sz + FileInfo.bytes;
            copyfile(src, dst);
        end
        
        fprintf('Copied %s, size : %d MB, cost %f sec\n', dirName, round(sz /1048576), toc(ticer));
    catch
         warning('ньньнь, we met some error!');
    end
end