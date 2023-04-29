clear all; % Clears old variables.
clc; % Clears command window.
close all; % Closes any open windows.

% Default to docked figures
% set(0,'DefaultFigureWindowStyle','docked')

% LaTeX stuff
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Stop ticks rotation
set(groot, 'defaultAxesXTickLabelRotationMode','manual')
set(groot, 'defaultAxesYTickLabelRotationMode','manual')
set(groot, 'defaultAxesZTickLabelRotationMode','manual')

rng(0,'twister'); % Fixed RNG seed for reproducibility

fontSize = 18; % Font size used in figures
saveFig = 1; % Boolean to turn on export of figures %%% WARNING: The export procedure will DELETE old figures in the folders WITHOUT asking!
prefix = 'EXPORT/'; % Export directory

L1 = 128; L2 = 64; % Size of box (L1,L2)
betaTemp = 0.6; % (inverse) Temperature
cStop = 0.1; % Stopping concentration of solvent
createDisks = 1; % Boolean to create disks before starting full dynamics.
domainMigration = 1; % Boolean to allow dynamics without evaporation and PBC after evaporation finished.

diskSteps = 10^8; % Number of iterations to run the system with partial dynamics to create the disks.
betaConst = 0.6; % Temperature when creating the disks
migrationSteps = 10^8; % Number of migration steps.

% Initial concentrations
c1 = 0.40;     % Acrylate; Blue
c2 = 0.15;     % Rubber; Yellow
cSol = 0.45;
c3 = 0.7;      % Ethylacetate; Red             
c4 = 0.3;      % Benzine; Black

cPrescribed = [c1 c2 cSol*c3 cSol*c4]; % Vector with the prescribed concentrations for debugging

% Interaction matrix components
J11 = 0;
J12 = 6;
J13 = 0.5;
J14 = 1.5;

J22 = -4;
J23 = 1;
J24 = 0.5;

J33 = 0;
J34 = 0.75;

J44 = 0;

JMatrix = [J11 J12 J13 J14; J12 J22 J23 J24; J13 J23 J33 J34; J14 J24 J34 J44]; % Construct symmetric interaction matrix
JMatrixInitial = [0 0 0 0; 0 J22 0 0; 0 0 0 0; 0 0 0 0]; % The interaction matrix used in the disks formation stage

% Error check
if c1 + c2 + cSol ~= 1
    error('Incorrect concentrations.')
end
if c3 + c4 ~= 1
    error('Incorrect solvent concentrations.')
end

% Initialise Lambda
Lambda = initialseLambda(L1,L2,c1,c2,cSol,c3,c4);
cInitial = getConcentrations(Lambda);

% Initial plot
figNum = 1;
figureIteration(figNum) = 0;
figureConcentration(figNum) = 0;
figure(figNum)
plotLambda(Lambda,fontSize)
title(['Initial condition; cSolvent: ' num2str(round(cInitial(3) + cInitial(4),2))])
pause(0.1)
figNum = figNum + 1;

% Produce disks
if createDisks == 1
    fprintf('Creating disks... ')
    for k = 1:diskSteps
        % Pick a radom site
        i = randi(L2);
        j = randi(L1);
    
        % Pick a random bond
        % 1 is horizontal
        % 2 is vertical
        bondType = randi(2);
    
        % Find bonded site
        if bondType == 1
            if j == L1
                jNeighbour = 1; % Right boundary
            else
                jNeighbour = j+1;
            end
            iNeighbour = i;
        else
            if i == 1
                iNeighbour = L2; % Top boundary
            else
                iNeighbour = i-1;
            end
            jNeighbour = j;
        end
    
        % Just move on if the two spins are equivalent.
        if Lambda(i,j) == Lambda(iNeighbour,jNeighbour)
            continue;
        end
    
        % Find neighbours of selected site
        [iAbove, jAbove, iBelow, jBelow, iLeft, jLeft, iRight, jRight] = findNeighbours(i,j,L1,L2);
    
        % Find neighbours of bonded site
        [iNeighbourAbove, jNeighbourAbove, iNeighbourBelow, jNeighbourBelow, iNeighbourLeft, jNeighbourLeft, iNeighbourRight, jNeighbourRight] = findNeighbours(iNeighbour,jNeighbour,L1,L2);
    
        % Compute current and proposed energetic landscapes
        sigmaCurrent = Lambda(i,j);
        sigmaNeighbour = Lambda(iNeighbour,jNeighbour);
        if bondType == 1
            ECurrent = JMatrixInitial(sigmaCurrent,Lambda(iAbove,jAbove)) + JMatrixInitial(sigmaCurrent,Lambda(iBelow,jBelow)) + JMatrixInitial(sigmaCurrent,Lambda(iLeft,jLeft)) + ...
                JMatrixInitial(sigmaNeighbour,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrixInitial(sigmaNeighbour,Lambda(iNeighbourBelow,jNeighbourBelow)) + JMatrixInitial(sigmaNeighbour,Lambda(iNeighbourRight,jNeighbourRight));
            EProposed = JMatrixInitial(sigmaNeighbour,Lambda(iAbove,jAbove)) + JMatrixInitial(sigmaNeighbour,Lambda(iBelow,jBelow)) + JMatrixInitial(sigmaNeighbour,Lambda(iLeft,jLeft)) + ...
                JMatrixInitial(sigmaCurrent,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrixInitial(sigmaCurrent,Lambda(iNeighbourBelow,jNeighbourBelow)) + JMatrixInitial(sigmaCurrent,Lambda(iNeighbourRight,jNeighbourRight));
        else
            ECurrent = JMatrixInitial(sigmaCurrent,Lambda(iBelow,jBelow)) + JMatrixInitial(sigmaCurrent,Lambda(iLeft,jLeft)) + JMatrixInitial(sigmaCurrent,Lambda(iRight,jRight)) + ...
                JMatrixInitial(sigmaNeighbour,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrixInitial(sigmaNeighbour,Lambda(iNeighbourLeft,jNeighbourLeft)) + JMatrixInitial(sigmaNeighbour,Lambda(iNeighbourRight,jNeighbourRight));
            EProposed = JMatrixInitial(sigmaNeighbour,Lambda(iBelow,jBelow)) + JMatrixInitial(sigmaNeighbour,Lambda(iLeft,jLeft)) + JMatrixInitial(sigmaNeighbour,Lambda(iRight,jRight)) + ...
                JMatrixInitial(sigmaCurrent,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrixInitial(sigmaCurrent,Lambda(iNeighbourLeft,jNeighbourLeft)) + JMatrixInitial(sigmaCurrent,Lambda(iNeighbourRight,jNeighbourRight));
        end
    
        % Move spins
        if EProposed < ECurrent
            Lambda(i,j) = sigmaNeighbour;
            Lambda(iNeighbour,jNeighbour) = sigmaCurrent;
        else
            randVar = rand(1);
            if randVar < exp(-betaConst*(EProposed - ECurrent))
                Lambda(i,j) = sigmaNeighbour;
                Lambda(iNeighbour,jNeighbour) = sigmaCurrent;
            end
        end
    end
    
    % Plot resulting disks formation
    figureIteration(figNum) = 0;
    figureConcentration(figNum) = 0;
    figure(figNum)
    plotLambda(Lambda,fontSize)
    title(['Disks formation; cSolvent: ' num2str(round(cInitial(3) + cInitial(4),2))])
    pause(0.1)
    figNum = figNum + 1;

    fprintf('Done!\n')
end

% Start full dynamics
stoppingCriterium = 0; % Initialise stopping variable
currentIteration = 0; % Count "iterations" from after disks formation stage
cPrevious = cInitial(3) + cInitial(4); % Specify initial solvent concentration
fprintf('Running dynamics... \n')
while stoppingCriterium == 0
    % Pick a radom site
    i = randi(L2);
    j = randi(L1);

    % Evporate solvent from top boundry and insert blue sites back
    if i == 1
        if Lambda(i,j) == 3 || Lambda(i,j) == 4
            randVar = rand(1);
            Lambda(i,j) = 1;
            currentIteration = currentIteration + 1;
            continue;
        end
    end

    % Pick a random bond
    % 1 is horizontal
    % 2 is vertical
    bondType = randi(2);

    % Find bonded site
    if bondType == 1
        if j == L1
            jNeighbour = 1; % Right boundary
        else
            jNeighbour = j+1;
        end
        iNeighbour = i;
    else
        if i == 1
            currentIteration = currentIteration + 1;
            continue; % Disallow migration accross top boundary.
        else
            iNeighbour = i-1;
        end
        jNeighbour = j;
    end

%     % Debug plot
%     figure(2)
%     plotSigma(sigma,fontSize)
%     hold on
%     plot(j,i,'*k','MarkerSize',40)
%     plot(jNeighbour,iNeighbour,'*w','MarkerSize',40)
%     hold off

    % Just move on if the two spins are equivalent.
    if Lambda(i,j) == Lambda(iNeighbour,jNeighbour)
        currentIteration = currentIteration + 1;
        continue;
    end

    % Find neighbours of selected site
    [iAbove, jAbove, iBelow, jBelow, iLeft, jLeft, iRight, jRight] = findNeighbours(i,j,L1,L2);

%     % Debug plot
%     figure(3)
%     plotSigma(sigma,fontSize)
%     hold on
%     plot(j,i,'*k','MarkerSize',40)
%     plot(jAbove,iAbove,'.k','MarkerSize',40)
%     plot(jBelow,iBelow,'.k','MarkerSize',40)
%     plot(jLeft,iLeft,'.k','MarkerSize',40)
%     plot(jRight,iRight,'.k','MarkerSize',40)
%     hold off

    % Find neighbours of bonded site
    [iNeighbourAbove, jNeighbourAbove, iNeighbourBelow, jNeighbourBelow, iNeighbourLeft, jNeighbourLeft, iNeighbourRight, jNeighbourRight] = findNeighbours(iNeighbour,jNeighbour,L1,L2);

%     % Debug plot
%     figure(4)
%     plotSigma(sigma,fontSize)
%     hold on
%     plot(jNeighbour,iNeighbour,'*w','MarkerSize',40)
%     plot(jNeighbourAbove,iNeighbourAbove,'.w','MarkerSize',40)
%     plot(jNeighbourBelow,iNeighbourBelow,'.w','MarkerSize',40)
%     plot(jNeighbourLeft,iNeighbourLeft,'.w','MarkerSize',40)
%     plot(jNeighbourRight,iNeighbourRight,'.w','MarkerSize',40)
%     hold off

    % Compute current and proposed energetic landscapes
    sigmaCurrent = Lambda(i,j);
    sigmaNeighbour = Lambda(iNeighbour,jNeighbour);
    if bondType == 1
        ECurrent = JMatrix(sigmaCurrent,Lambda(iAbove,jAbove)) + JMatrix(sigmaCurrent,Lambda(iBelow,jBelow)) + JMatrix(sigmaCurrent,Lambda(iLeft,jLeft)) + ...
            JMatrix(sigmaNeighbour,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourBelow,jNeighbourBelow)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourRight,jNeighbourRight));
        EProposed = JMatrix(sigmaNeighbour,Lambda(iAbove,jAbove)) + JMatrix(sigmaNeighbour,Lambda(iBelow,jBelow)) + JMatrix(sigmaNeighbour,Lambda(iLeft,jLeft)) + ...
            JMatrix(sigmaCurrent,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaCurrent,Lambda(iNeighbourBelow,jNeighbourBelow)) + JMatrix(sigmaCurrent,Lambda(iNeighbourRight,jNeighbourRight));
    else
        ECurrent = JMatrix(sigmaCurrent,Lambda(iBelow,jBelow)) + JMatrix(sigmaCurrent,Lambda(iLeft,jLeft)) + JMatrix(sigmaCurrent,Lambda(iRight,jRight)) + ...
            JMatrix(sigmaNeighbour,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourLeft,jNeighbourLeft)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourRight,jNeighbourRight));
        EProposed = JMatrix(sigmaNeighbour,Lambda(iBelow,jBelow)) + JMatrix(sigmaNeighbour,Lambda(iLeft,jLeft)) + JMatrix(sigmaNeighbour,Lambda(iRight,jRight)) + ...
            JMatrix(sigmaCurrent,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaCurrent,Lambda(iNeighbourLeft,jNeighbourLeft)) + JMatrix(sigmaCurrent,Lambda(iNeighbourRight,jNeighbourRight));
    end

    % Move spins
    if EProposed < ECurrent
        Lambda(i,j) = sigmaNeighbour;
        Lambda(iNeighbour,jNeighbour) = sigmaCurrent;
    else
        randVar = rand(1);
        if randVar < exp(-betaTemp*(EProposed - ECurrent))
            Lambda(i,j) = sigmaNeighbour;
            Lambda(iNeighbour,jNeighbour) = sigmaCurrent;
        end
    end

%     % Debug plot
%     figure(5)
%     plotSigma(sigma,fontSize)
%     hold on
%     plot(j,i,'*k','MarkerSize',40)
%     plot(jNeighbour,iNeighbour,'*w','MarkerSize',40)
%     hold off

    % Give nice feedback to user and plot frames
    if mod(currentIteration,10000) == 0
        cCurrent = getConcentrations(Lambda);

        if round(cCurrent(3) + cCurrent(4),2) <= cPrevious - 0.01
            figureIteration(figNum) = currentIteration;
            figureConcentration(figNum) = round(cCurrent(3) + cCurrent(4),2);

            figure(figNum)
            plotLambda(Lambda,fontSize)
            title(['Iteration: ' num2str(round(currentIteration/10^6,0)) '$\cdot 10^6$' '; cSolvent: ' num2str(round(cCurrent(3) + cCurrent(4),2))])
            pause(0.1)
            figNum = figNum + 1;
            
            cPrevious = round(cCurrent(3) + cCurrent(4),2);

            fprintf('Completion: %.0f%% \n', 100*(cSol-(cCurrent(3) + cCurrent(4)))/(cSol-cStop))
        end
        
        if cCurrent(3) + cCurrent(4) <= cStop
            figureIteration(figNum) = currentIteration;
            figureConcentration(figNum) = round(cCurrent(3) + cCurrent(4),2);

            figure(figNum)
            plotLambda(Lambda,fontSize)
            title(['Iteration: ' num2str(round(currentIteration/10^6,0)) '$\cdot 10^6$' '; cSolvent: ' num2str(round(cCurrent(3) + cCurrent(4),2))])
            pause(0.1)
            figNum = figNum + 1;

            fprintf('Completion: %.0f%% \n', 100*(cSol-(cCurrent(3) + cCurrent(4)))/(cSol-cStop))
            
            break;
        end
    end

    currentIteration = currentIteration + 1; % Increment itaration variable.
end
fprintf('Done!\n')

% Migrate domains
if domainMigration == 1
    fprintf('Migrating domains... \n')
    for k = 1:migrationSteps
        % Pick a radom site
        i = randi(L2);
        j = randi(L1);
    
        % Pick a random bond
        % 1 is horizontal
        % 2 is vertical
        bondType = randi(2);
    
        % Find bonded site
        if bondType == 1
            if j == L1
                jNeighbour = 1; % Right boundary
            else
                jNeighbour = j+1;
            end
            iNeighbour = i;
        else
            if i == 1
                iNeighbour = L2; % Top boundary
            else
                iNeighbour = i-1;
            end
            jNeighbour = j;
        end
    
        % Just move on if the two spins are equivalent.
        if Lambda(i,j) == Lambda(iNeighbour,jNeighbour)
            continue;
        end
    
        % Find neighbours of selected site
        [iAbove, jAbove, iBelow, jBelow, iLeft, jLeft, iRight, jRight] = findNeighbours(i,j,L1,L2);
    
        % Find neighbours of bounded site
        [iNeighbourAbove, jNeighbourAbove, iNeighbourBelow, jNeighbourBelow, iNeighbourLeft, jNeighbourLeft, iNeighbourRight, jNeighbourRight] = findNeighbours(iNeighbour,jNeighbour,L1,L2);
    
        % Compute current and proposed energetic landscapes
        sigmaCurrent = Lambda(i,j);
        sigmaNeighbour = Lambda(iNeighbour,jNeighbour);
        if bondType == 1
            ECurrent = JMatrix(sigmaCurrent,Lambda(iAbove,jAbove)) + JMatrix(sigmaCurrent,Lambda(iBelow,jBelow)) + JMatrix(sigmaCurrent,Lambda(iLeft,jLeft)) + ...
                JMatrix(sigmaNeighbour,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourBelow,jNeighbourBelow)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourRight,jNeighbourRight));
            EProposed = JMatrix(sigmaNeighbour,Lambda(iAbove,jAbove)) + JMatrix(sigmaNeighbour,Lambda(iBelow,jBelow)) + JMatrix(sigmaNeighbour,Lambda(iLeft,jLeft)) + ...
                JMatrix(sigmaCurrent,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaCurrent,Lambda(iNeighbourBelow,jNeighbourBelow)) + JMatrix(sigmaCurrent,Lambda(iNeighbourRight,jNeighbourRight));
        else
            ECurrent = JMatrix(sigmaCurrent,Lambda(iBelow,jBelow)) + JMatrix(sigmaCurrent,Lambda(iLeft,jLeft)) + JMatrix(sigmaCurrent,Lambda(iRight,jRight)) + ...
                JMatrix(sigmaNeighbour,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourLeft,jNeighbourLeft)) + JMatrix(sigmaNeighbour,Lambda(iNeighbourRight,jNeighbourRight));
            EProposed = JMatrix(sigmaNeighbour,Lambda(iBelow,jBelow)) + JMatrix(sigmaNeighbour,Lambda(iLeft,jLeft)) + JMatrix(sigmaNeighbour,Lambda(iRight,jRight)) + ...
                JMatrix(sigmaCurrent,Lambda(iNeighbourAbove,jNeighbourAbove)) + JMatrix(sigmaCurrent,Lambda(iNeighbourLeft,jNeighbourLeft)) + JMatrix(sigmaCurrent,Lambda(iNeighbourRight,jNeighbourRight));
        end
    
        % Move spins
        if EProposed < ECurrent
            Lambda(i,j) = sigmaNeighbour;
            Lambda(iNeighbour,jNeighbour) = sigmaCurrent;
        else
            randVar = rand(1);
            if randVar < exp(-betaTemp*(EProposed - ECurrent))
                Lambda(i,j) = sigmaNeighbour;
                Lambda(iNeighbour,jNeighbour) = sigmaCurrent;
            end
        end
    end
    
    % Plot migration result
    figureIteration(figNum) = currentIteration + migrationSteps;
    figureConcentration(figNum) = round(cCurrent(3) + cCurrent(4),2);
    figure(figNum)
    plotLambda(Lambda,fontSize)
    title(['After migration; cSolvent: ' num2str(round(cCurrent(3) + cCurrent(4),2))])
    pause(0.1)
    figNum = figNum + 1;

    fprintf('Done!\n')
end


%% Save figures
fileName = ['L1_' num2str(L1) '-L2_' num2str(L2) '-beta_' num2str(betaTemp) '-c1_' num2str(c1) '-c2_' num2str(c2) '-cSol_' num2str(cSol) '-c3_' num2str(c3) '-c4_' num2str(c4) '-cStop_' num2str(cStop)];

if saveFig == 1
    fprintf('Exporting figures... \n')
    saveFigures(prefix,fileName,createDisks,domainMigration,figureIteration,figureConcentration)
    fprintf('Done!\n')
end





% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  Functions.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function Lambda = initialseLambda(L1,L2,c1,c2,cSol,c3,c4)
% Mapping
% 1 is AC       (blue)
% 2 is Resin    (yellow)
% 3 is Ethyl    (red)
% 4 is Benzine  (green)

% Intialise spin matrix
Lambda = zeros(L2,L1);
for i = 1:L2
    for j = 1:L1
        rand1 = rand(1);

        if rand1 < c1
            Lambda(i,j) = 1;
        elseif rand1 < c1 + c2
            Lambda(i,j) = 2;
        else
            rand2 = rand(1);
            if rand2 < c3
                Lambda(i,j) = 3;
            else
                Lambda(i,j) = 4;
            end
        end
    end
end
end

% Plot the spin matrix
function plotLambda(Lambda,fontSize)
map = [0 0 1; 1 1 0; 1 0 0; 0 1 0]; % Map 1 to blue, 2 to yellow, 3 to red and 4 to green
imagesc(Lambda)
daspect([1 1 1]) % Fix aspect ratio
colormap(map);
xticks(round(linspace(1,size(Lambda,2),8),0));
yticks(round(linspace(1,size(Lambda,1),8),0));
set(gca,'FontSize',fontSize)
end

% Compute concentrations
function concVec = getConcentrations(Lambda)
num1 = sum(Lambda(:) == 1);
num2 = sum(Lambda(:) == 2);
num3 = sum(Lambda(:) == 3);
num4 = sum(Lambda(:) == 4);

c1_current = num1/numel(Lambda);
c2_current = num2/numel(Lambda);
c3_current = num3/numel(Lambda);
c4_current = num4/numel(Lambda);

concVec = [c1_current c2_current c3_current c4_current];
end

% Find the four nearest-neighbours of a given site
function [iAbove, jAbove, iBelow, jBelow, iLeft, jLeft, iRight, jRight] = findNeighbours(i,j,L1,L2)
if i == 1
    iAbove = L2; % Top bondary
else
    iAbove = i - 1;
end
jAbove = j;

if i == L2
    iBelow = 1; % Bottom boundary
else
    iBelow = i + 1;
end
jBelow = j;

if j == 1
    jLeft = L1; % Left boundary
else
    jLeft = j - 1;
end
iLeft = i;

if j == L1
    jRight = 1; % Right boundary
else
    jRight = j + 1;
end
iRight = i;
end

% Save figures automagically
function saveFigures(prefix,fileName,createDisks,domainMigration,figureIteration,figureConcentration)
figHandles = findall(0,'Type','figure'); % Get all figure handles
% Sort in order of appearance in script
for i = 1:numel(figHandles)
    figNums(i) = figHandles(i).Number;
end
[~,index] = sort(figNums);
figHandles = figHandles(index);

for i = 1:numel(figHandles)
    nameVec{i} = ['iteration_' num2str(figureIteration(i)) '-cSolvent_' num2str(figureConcentration(i))];
end

% Create directory
fileName = regexprep(fileName,'.mat','','emptymatch');
dirName = [prefix fileName '_FIGS'];
dirExists = dir(dirName);
if isempty(dirExists) == 0
    fprintf('Deleting old directory... ')
    rmdir(dirName,'s')
    fprintf('Done!\n')
end
fprintf('Creating directory... ')
mkdir(dirName)
fprintf('Done!\n')

for i = 1:numel(figHandles)
    set(figHandles(i),'Units','pixels');
    set(figHandles(i),'Position', [0 0 600 350])    
    set(figHandles(i),'color','w');
    pause(0.5)
    
    if i == 1
        filename = [dirName '/initial_' nameVec{i} '.png'];
    elseif i == 2 && createDisks == 1
        filename = [dirName '/disks_' nameVec{i} '.png'];
    elseif i == numel(figHandles) && domainMigration == 1
        filename = [dirName '/migrated_' nameVec{i} '.png'];
    else
        filename = [dirName '/' nameVec{i} '.png'];
    end
    frame = getframe(figHandles(i));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'png');

    if i == 1
        filename = [dirName '/initial_' nameVec{i} '.pdf'];
    elseif i == 2 && createDisks == 1
        filename = [dirName '/disks_' nameVec{i} '.pdf'];
    elseif i == numel(figHandles) && domainMigration == 1
        filename = [dirName '/migrated_' nameVec{i} '.pdf'];
    else
        filename = [dirName '/' nameVec{i} '.pdf'];
    end
    set(figHandles(i),'Units','Inches');
    pos = get(figHandles(i),'Position');
    set(figHandles(i),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(figHandles(i),filename,'-dpdf','-r0')

    if i == 1
        filename = [dirName '/initial_' nameVec{i} '.figure'];
    elseif i == 2 && createDisks == 1
        filename = [dirName '/disks_' nameVec{i} '.figure'];
    elseif i == numel(figHandles) && domainMigration == 1
        filename = [dirName '/migrated_' nameVec{i} '.figure'];
    else
        filename = [dirName '/' nameVec{i} '.figure'];
    end
    savefig(figHandles(i),filename)
end
end