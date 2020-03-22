function [numOfPickedPar,numOfPickedNoise] = picking_from_scoring_mat(logTestN,patchSzFun,patchSzN,patchSzPickBox,num_of_particles,num_of_noise_images,mgBigSz,mgScale,coordinatsPathParticle,coordinatsPathNoise,microName,thresh)
% Input parameters:
% logTestN              the scoring matrix to pick from.
% patchSzFun            the patch size used to construct the eigenfunctions 
% patchSzFun            the patch size used to center the particles.
% patchSzPickBox        particle box size to use
% num_of_particles      amount of particle to pick. -1 to pick all.
% num_of_noise_images   amount of noise images to pick.
% mgBigSz               size of micrograph
% mgScale               downsampling scalar
% boxPathParticle       path to the particle coordinate files
% boxPathNoise          path to the noise images coordinate files
% microName             name of micrograph  
% thresh                threshold for the picking out of the scoring mat.

% Output parameters:
% coordinate files in coordinatsPathParticle and coordinatsPathNoise.
% numOfPickedPar        number of picked particles.
% numOfPickedPar        number of picked noise images.

    % setting parameters
    idx = 1:size(logTestN,1);
    [colIdx,rowIdx] = meshgrid(idx);
    rDel = floor(patchSzPickBox);
    sz = size(logTestN);
    numOfPickedPar = 0;
    numOfPickedNoise = 0;
    %% particle picking from scoring matrix
    scoringMat = logTestN;
    logMax = max(logTestN(:));
    coordinatsPathParticleBox = [coordinatsPathParticle,'/box'];
    coordinatsPathParticleStar = [coordinatsPathParticle,'/star'];
    if ~exist(coordinatsPathParticleBox, 'dir')
       mkdir(coordinatsPathParticleBox)
    end
    if ~exist(coordinatsPathParticleStar, 'dir')
       mkdir(coordinatsPathParticleStar)
    end
  
    particlesCordinateBox = fopen(fullfile(coordinatsPathParticleBox,[microName,'.box']),'w');
    particlesCordinateStar = fopen(fullfile(coordinatsPathParticleStar,[microName,'.star']),'w');
    % format relion star file
    fprintf(particlesCordinateStar,'%s\n','data_');
    fprintf(particlesCordinateStar,'%s\n','');
    fprintf(particlesCordinateStar,'%s\n','loop_');
    fprintf(particlesCordinateStar,'%s\n','_rlnCoordinateX #1');
    fprintf(particlesCordinateStar,'%s\n','_rlnCoordinateY #2');
    fprintf(particlesCordinateStar,'%s\n','_rlnAutopickFigureOfMerit #3');
    if num_of_particles == -1 % pick all particles
        pMax = thresh+1; % initilizing
        while pMax>thresh
            [pMax,I] = max(scoringMat(:));
            if pMax<=thresh 
                break
            end
            [i_row, i_col] = ind2sub(sz,I);
            i_rowPatch = (i_row-1)+ceil((patchSzFun+patchSzN)/2); i_colPatch = (i_col-1)+ceil((patchSzFun+patchSzN)/2); % change index to middle of the patch.  
            rowIdxB = rowIdx-i_row;
            colIdxB = colIdx-i_col;
            Rsquare = rowIdxB.^2+colIdxB.^2;
            scoringMat(Rsquare<=(rDel^2)) = -inf;
            fprintf(particlesCordinateBox,'%i\t%i\t%i\t%i\n',(1/mgScale)*(i_colPatch - floor(patchSzPickBox/2)),(mgBigSz(1)+1)-(1/mgScale)*(i_rowPatch + floor(patchSzPickBox/2)),(1/mgScale)*patchSzPickBox,(1/mgScale)*patchSzPickBox);
            fprintf(particlesCordinateStar,'%i\t%i\t%i\n',(1/mgScale)*i_colPatch,(mgBigSz(1)+1)-(1/mgScale)*i_rowPatch,pMax/logMax);
            numOfPickedPar = numOfPickedPar + 1;
        end
    else
        iterPick = 1;  
        pMax = thresh+1;
        while and(iterPick <= num_of_particles,pMax>thresh)
            [pMax,I] = max(scoringMat(:));
            if pMax<=thresh 
                break
            end
            [i_row, i_col] = ind2sub(sz,I);
            i_rowPatch = (i_row-1)+ceil((patchSzFun+patchSzN)/2); i_colPatch = (i_col-1)+ceil((patchSzFun+patchSzN)/2); % change index to middle of the patch.  
            rowIdxB = rowIdx-i_row;
            colIdxB = colIdx-i_col;
            Rsquare = rowIdxB.^2+colIdxB.^2;
            scoringMat(Rsquare<=(rDel^2)) = -inf;
            fprintf(particlesCordinateBox,'%i\t%i\t%i\t%i\n',(1/mgScale)*(i_colPatch - floor(patchSzPickBox/2)),(mgBigSz(1)+1)-(1/mgScale)*(i_rowPatch + floor(patchSzPickBox/2)),(1/mgScale)*patchSzPickBox,(1/mgScale)*patchSzPickBox);
            fprintf(particlesCordinateStar,'%i\t%i\n',(1/mgScale)*i_colPatch,(mgBigSz(1)+1)-(1/mgScale)*i_rowPatch);
            iterPick = iterPick + 1;
            numOfPickedPar = numOfPickedPar + 1;
        end
    end
    fclose(particlesCordinateBox); 
    fclose(particlesCordinateStar); 
    %% noise images picking from the scoring matrix
    if num_of_noise_images ~= 0 % pick noise images
        scoringMat = logTestN;
        coordinatsPathNoiseBox = [coordinatsPathNoise,'/box'];
        coordinatsPathNoiseStar = [coordinatsPathNoise,'/star'];
        if ~exist(coordinatsPathNoiseBox, 'dir')
            mkdir(coordinatsPathNoiseBox)
        end
        if ~exist(coordinatsPathNoiseStar, 'dir')
            mkdir(coordinatsPathNoiseStar)
        end
        noiseImagesCordinateBox = fopen(fullfile(coordinatsPathNoiseBox,[microName,'.box']),'w');
        noiseImagesCordinateStar = fopen(fullfile(coordinatsPathNoiseStar,[microName,'.star']),'w');
        % format relion star file
        fprintf(noiseImagesCordinateStar,'%s\n','data_');
        fprintf(noiseImagesCordinateStar,'%s\n','');
        fprintf(noiseImagesCordinateStar,'%s\n','loop_');
        fprintf(noiseImagesCordinateStar,'%s\n','_rlnCoordinateX #1');
        fprintf(noiseImagesCordinateStar,'%s\n','_rlnCoordinateY #2');
        iterPick = 1;  
        pMin = thresh-1;
        while and(iterPick <= num_of_noise_images,pMin<thresh)
            [pMin,I] = min(scoringMat(:));
            if pMin>=thresh 
                break
            end
            [i_row, i_col] = ind2sub(sz,I);
            i_rowPatch = (i_row-1)+ceil((patchSzFun+patchSzN)/2); i_colPatch = (i_col-1)+ceil((patchSzFun+patchSzN)/2); % change index to middle of the patch.  
            rowIdxB = rowIdx-i_row;
            colIdxB = colIdx-i_col;
            Rsquare = rowIdxB.^2+colIdxB.^2;
            scoringMat(Rsquare<=(rDel^2)) = inf;
            fprintf(noiseImagesCordinateBox,'%i\t%i\t%i\t%i\n',(1/mgScale)*(i_colPatch - floor(patchSzPickBox/2)),(mgBigSz(1)+1)-(1/mgScale)*(i_rowPatch + floor(patchSzPickBox/2)),(1/mgScale)*patchSzPickBox,(1/mgScale)*patchSzPickBox);
            fprintf(noiseImagesCordinateStar,'%i\t%i\n',(1/mgScale)*i_colPatch,(mgBigSz(1)+1)-(1/mgScale)*i_rowPatch);
            iterPick = iterPick + 1;
            numOfPickedNoise = numOfPickedNoise + 1;
        end
        fclose(noiseImagesCordinateBox); 
        fclose(noiseImagesCordinateStar); 
    end
    
end