% Copyright (c) 2024 Johannes Sieberer
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% Import cases to be analyzed
 cases = importPatients("DataPatients.csv");
 cases = transpose(cases.Identifier);
%% Iterate trhough all patients
for o = 1:length(cases)
    casename = int2str(cases(o));
    canal = readSurfaceMesh(strcat("Models/","Canal_Surf_",casename,".stl")); %Stl model of the canal
    cortex = readSurfaceMesh(strcat("Models/","Total_Surf_",casename,".stl"));%Stl model of the cortex
    %% Get Points placed on the proximal radius
    points = importPoints(strcat("Measure/",casename,"_R_Prox_RadiusMeasurements.csv"));
    radialHead = [points.P1xMm(points.Name == "RadialHead");points.P1yMm(points.Name == "RadialHead");points.P1zMm(points.Name == "RadialHead")];
    radialShaft = [points.P1xMm(points.Name == "RadialShaft");points.P1yMm(points.Name == "RadialShaft");points.P1zMm(points.Name == "RadialShaft")];
    radialCut = [points.P1xMm(points.Name == "RadialCut");points.P1yMm(points.Name == "RadialCut");points.P1zMm(points.Name == "RadialCut")];
    %% Get Plane Normal
    normalPl = (radialShaft - radialHead)/norm(radialShaft - radialHead);
    %% Get meshes and translate
    meshFull.vertices = cortex.Vertices - transpose(radialHead);
    meshFull.faces = cortex.Faces;
    meshCanal.vertices = canal.Vertices- transpose(radialHead);
    meshCanal.faces = canal.Faces;
    radialShaft =  radialShaft - radialHead;
    radialCut = radialCut - radialHead;
    %% get rot matrix to rotate cut into plane
    zvector = [0 0 1];
    v = cross(normalPl, zvector);
    skew=[0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];
    rotmatrix = eye(3) + skew + skew^2 * 1/(1+dot(normalPl,zvector));
    
    %% Rotate all meshes into plane
    for i = 1:height(meshFull.vertices)
    meshFull.vertices(i,:) = transpose(rotmatrix * transpose(meshFull.vertices(i,:)));
    end
    for i = 1:height(meshCanal.vertices)
    meshCanal.vertices(i,:) = transpose(rotmatrix * transpose(meshCanal.vertices(i,:)));
    end
    radialHead = rotmatrix * radialHead;
    radialCut= rotmatrix * radialCut;
    
    %% Get intersections  Canal
    MeasureDistance = 60; %Define measure interval
    heightp = radialCut(3); 
    heightm = heightp+MeasureDistance;
    interval = 0.1; 
    heights = (heightp:interval:heightm);
    numSlices = length(heights);
    
    ellFull(numSlices) = struct('a', [], 'b', [], 'phi', [], 'X0', [], 'Y0', [], 'X0_in', [], 'Y0_in', [], 'long_axis', [], 'short_axis', [], 'status', [],'error', []);
    ellCanal(numSlices) = struct('a', [], 'b', [], 'phi', [], 'X0', [], 'Y0', [], 'X0_in', [], 'Y0_in', [], 'long_axis', [], 'short_axis', [], 'status', [],'error', []);
    
    eMidPoint = zeros(numSlices,2);
    eFullMidPoint = zeros(numSlices,2);
    eCanal = zeros(numSlices,3); eFull = NaN(numSlices,3);
    errorCanalCirc = zeros(numSlices,1); errorCortCirc = zeros(numSlices,1);
    for (n = 1:numSlices)
        [Surf12] = PlaneIntersection(meshCanal,heights(n));
        [Surf13] = PlaneIntersection(meshFull,heights(n));
        try
            S2=Surf13; S1=Surf12;
            if sum(size(S1)) > 0
                ellCanal(n) = fit_ellipse2(S1(:,1),S1(:,2));
                [~,errorCanalCirc(n)] = circfit(S1(:,1),S1(:,2));
                if(ellCanal(n).a >0)
                    eCanal(n,:) = [ellCanal(n).long_axis,ellCanal(n).short_axis,ellCanal(n).error];
                    eMidPoint(n,:)=[ellCanal(n).X0_in,ellCanal(n).Y0_in];
                end
            else 
                rCanal(n,:) = [nan,nan,nan];
                eCanal(n,:) = [0,0];
            end
            if sum(size(S2)) > 0
                ellFull(n) = fit_ellipse2(S2(:,1),S2(:,2));
                [~,errorCortCirc(n)] = circfit(S2(:,1),S2(:,2));
                if(ellCanal(n).a >0)
                    eFull(n,:) = [ellFull(n).long_axis,ellFull(n).short_axis,ellFull(n).error];
                    eFullMidPoint(n,:)=[ellFull(n).X0_in,ellFull(n).Y0_in];
                end
            else 
                rFull(n,:) = [nan,nan,nan];
                eFull(n,:) = [0,0];
            end
        end
    end
    %% Prepare results for .csv file
    strtrim(cellstr(num2str(transpose(heights),'%.2f')));
    results.realdistance = strtrim(cellstr(num2str(transpose(heights),'%.2f')));
    results.eCanalL = strtrim(cellstr(num2str(eCanal(:,1),'%.2f')));
    results.eCanalS = strtrim(cellstr(num2str(eCanal(:,2),'%.2f')));
    results.eFullL = strtrim(cellstr(num2str(eFull(:,1),'%.2f')));
    results.eFullS = strtrim(cellstr(num2str(eFull(:,2),'%.2f')));
    results.eCanalMidPointX = strtrim(cellstr(num2str( eMidPoint(:,1),'%.2f')));
    results.eCanalMidPointY = strtrim(cellstr(num2str(eMidPoint(:,2),'%.2f')));
    results.eFullMidPointX = strtrim(cellstr(num2str( eFullMidPoint(:,1),'%.2f')));
    results.eFullMidPointY = strtrim(cellstr(num2str(eFullMidPoint(:,2),'%.2f')));
    results.eCanalRMS = strtrim(cellstr(num2str(eCanal(:,3),'%.5f')));
    results.eFullRMS = strtrim(cellstr(num2str(eFull(:,3),'%.5f')));
    results.cCanalRMS = strtrim(cellstr(num2str(errorCanalCirc(:),'%.2f')));
    results.cFullRMS = strtrim(cellstr(num2str(errorCortCirc(:),'%.2f')));
    %% Write .csv file
    T=struct2table(results);
    writetable(T,strcat("Results/",casename,".csv")); 
    strcat("Results/",casename,".csv");
end