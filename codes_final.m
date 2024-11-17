
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the cancer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Download the prostrate cancer data from TCGA portal
% Load the TCGA prostate cancer data
load('TCGAprad.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Normalization and Filtration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract and convert the count data to double
geneExpressionCounts = double(TCGAprad(3:end, 2:end-1));
pseudoReferenceSample = geomean(geneExpressionCounts, 2);
nonZeroPseudoRefSample = pseudoReferenceSample > 0;

% Normalize the counts using size factors
ratios = bsxfun(@rdivide, geneExpressionCounts(nonZeroPseudoRefSample, :), pseudoReferenceSample(nonZeroPseudoRefSample));
sizeFactors = median(ratios, 1);
normalizedCounts = bsxfun(@rdivide, geneExpressionCounts, sizeFactors);

% Extract gene names
geneNames = TCGAprad(3:end, 555);
TCGAprad = TCGAprad(:, 2:end-1);

% Identify normal and tumor sample indices
normalSampleIndices = find(TCGAprad(2, :) == 'Solid Tissue Normal');
tumorSampleIndices = find(TCGAprad(2, :) == 'Primary Tumor');

% Separate normal and tumor samples
normalSamples = [geneNames, normalizedCounts(:, normalSampleIndices)];
tumorSamples = [geneNames, normalizedCounts(:, tumorSampleIndices)];

% Find lowly expressed genes from both samples with a cutoff of 10
lowExpressionCutoff = 10;
lowExpressionData = [];

for i = 1:size(normalSamples, 1)
    lowNormalExpression = length(find(double(normalSamples(i, 2:end)) < lowExpressionCutoff)) / size(normalSamples, 2);
    lowTumorExpression = length(find(double(tumorSamples(i, 2:end)) < lowExpressionCutoff)) / size(tumorSamples, 2);
    lowExpressionData = [lowExpressionData; lowNormalExpression, lowTumorExpression];
end

% Plot the lowly expressed genes
numLowNormal = length(find(lowExpressionData(:, 1) >= 0.7));
numLowTumor = length(find(lowExpressionData(:, 2) >= 0.7));
numLowBoth = length(find(lowExpressionData(:, 1) >= 0.7 & lowExpressionData(:, 2) >= 0.7));
barData = [numLowNormal, numLowTumor, numLowBoth];

bar(barData)
for i = 1:numel(barData)
    text(i, barData(i), num2str(barData(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold')
end

grid on
box on
xticks(1:3)
xticklabels(["Normal", "Tumor", "Both"])
xlabel('Categories')
ylabel('Lowly Expressed Genes')
ylim([0 41000])
set(gca, 'FontWeight', 'bold', 'FontSize', 10, 'LineWidth', 2)

% Remove the lowly expressed genes
lowExpressionIndices = find(lowExpressionData(:, 1) >= 0.7 & lowExpressionData(:, 2) >= 0.7);
normalSamples(lowExpressionIndices, :) = [];
tumorSamples(lowExpressionIndices, :) = [];
geneNames(lowExpressionIndices, :) = [];


% Remove gene names from the matrices and handle outliers
normalSamples = double(normalSamples(:, 2:end));
tumorSamples = double(tumorSamples(:, 2:end));

normalSamples = filloutliers(normalSamples, "nearest", "median");
tumorSamples = filloutliers(tumorSamples, "nearest", "median");

% Merge duplicate genes
duplicateGeneNames = string(tabulate(geneNames(:, 1)));
duplicateGeneNames = sortrows(duplicateGeneNames, -2);
duplicateNames = duplicateGeneNames(double(duplicateGeneNames(:, 2)) > 1, 1);

removeIndices = [];
for i = 1:size(duplicateNames, 1)
    duplicateIndices = find(strcmp(geneNames(:, 1), duplicateNames(i, 1)));
    normalSamples(duplicateIndices(1), :) = mean(double(normalSamples(duplicateIndices, :)), 1);
    tumorSamples(duplicateIndices(1), :) = mean(double(tumorSamples(duplicateIndices, :)), 1);
    removeIndices = [removeIndices; duplicateIndices(2:end)];
end

normalSamples(removeIndices, :) = [];
tumorSamples(removeIndices, :) = [];
geneNames(removeIndices, :) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify the differentially expressed genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate fold change and perform statistical tests
meanNormal = mean(normalSamples, 2);
meanTumor = mean(tumorSamples, 2);
foldChange = meanTumor ./ meanNormal;
pValues = mattest(tumorSamples, normalSamples);
log2FoldChange = log2(foldChange);

% Volcano plot
negLog10PValues = -log10(pValues);

figure;

% Plot non-significant points in gray
scatter(log2FoldChange(pValues > 0.05), negLog10PValues(pValues > 0.05), 10, [0.8 0.8 0.8], 'filled');
hold on

% Plot significant up-regulated points in red
scatter(log2FoldChange(pValues <= 0.05 & log2FoldChange >= 1), negLog10PValues(pValues <= 0.05 & log2FoldChange >= 1), 40, 'r', 'filled');
hold on

% Plot significant down-regulated points in blue
scatter(log2FoldChange(pValues <= 0.05 & log2FoldChange <= -1), negLog10PValues(pValues <= 0.05 & log2FoldChange <= -1), 40, 'b', 'filled');

xlim([-7 7]);
ylim([0 max(negLog10PValues) + 1]);
xlabel('Log2 Fold Change');
ylabel('-Log10 P-value');

yline(-log10(0.05), '--k', 'p = 0.05', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
xline(1, '--r', 'Log2FC = 1', 'FontSize', 16,'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
xline(-1, '--b', 'Log2FC = -1', 'FontSize', 16,'LineWidth', 1, 'LabelHorizontalAlignment', 'left');

grid on;
%box on;
set(gca, 'FontSize', 18, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
legend({'Non-significant', 'Significant Up-regulated', 'Significant Down-regulated'}, 'Location', 'northwest');
set(gcf, 'Position', [100, 100, 800, 600]);


% Identify up- and down-regulated genes
upRegulatedIndices = find(log2FoldChange >= 1 & pValues <= 0.05);
downRegulatedIndices = find(log2FoldChange <= -1 & pValues <= 0.05);

upRegulatedGenes = [geneNames(upRegulatedIndices, 1), log2FoldChange(upRegulatedIndices)];
downRegulatedGenes = [geneNames(downRegulatedIndices, 1), log2FoldChange(downRegulatedIndices)];

countMatrix = [geneNames, normalSamples, tumorSamples];
degMatrix = [upRegulatedGenes; downRegulatedGenes];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction and analysis of the PPI network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load undirected network
load('final_undirected_network.mat')

% Find matching rows in the network for differentially expressed genes
logicalIndex = ismember(final_undirected_network(:, 1), degMatrix) | ismember(final_undirected_network(:, 2), degMatrix);
matchingIndices = find(logicalIndex);
constructedNetwork = final_undirected_network(matchingIndices, :);

% Remove self-loops and duplicate edges
uniqueEdges = {};

for i = 1:size(constructedNetwork, 1)
    edge = sort(constructedNetwork(i, :));
    if ~strcmp(edge{1}, edge{2})
        uniqueEdges = [uniqueEdges; edge];
    end
end

uniqueEdges = unique(uniqueEdges, 'rows', 'stable');

nodes = unique(uniqueEdges(:));
nodeMap = containers.Map(nodes, 1:length(nodes));

adjacencyMatrix = zeros(length(nodes));

for i = 1:size(uniqueEdges, 1)
    idx1 = nodeMap(uniqueEdges{i, 1});
    idx2 = nodeMap(uniqueEdges{i, 2});
    adjacencyMatrix(idx1, idx2) = 1;
    adjacencyMatrix(idx2, idx1) = 1;
end

% Analyze the network
graphNetwork = graph(adjacencyMatrix);

degreeRanks = centrality(graphNetwork, 'degree');
degreeRanks = [nodes, degreeRanks];

betweennessRanks = centrality(graphNetwork, 'betweenness');
betweennessRanks = [nodes, betweennessRanks];

closenessRanks = centrality(graphNetwork, 'closeness');
closenessRanks = [nodes, closenessRanks];

[~, sortedDegreeIndices] = sort(double(degreeRanks(:, 2)), 'descend');
sortedDegreeRanks = degreeRanks(sortedDegreeIndices, :);

[~, sortedBetweennessIndices] = sort(double(betweennessRanks(:, 2)), 'descend');
sortedBetweennessRanks = betweennessRanks(sortedBetweennessIndices, :);

[~, sortedClosenessIndices] = sort(double(closenessRanks(:, 2)), 'descend');
sortedClosenessRanks = closenessRanks(sortedClosenessIndices, :);

topNodes = intersect(sortedDegreeRanks(1:180, 1), intersect(sortedClosenessRanks(1:180, 1), sortedBetweennessRanks(1:180, 1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRISPR analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load CRISPR DATA for Prostrate
load('CRISPRDepMap.mat')

colnames = CRISPRDepMap(1,:);
total_cell_lines = size(CRISPRDepMap,1)-1;
[top_nodes_effect,b,~] = intersect(colnames,topNodes); 

top_nodes_effect(:,2) = b;

for i = 1:size(top_nodes_effect,1)
    j = double(CRISPRDepMap(2:end,double(top_nodes_effect(i,2))));
    top_nodes_effect(i,3) = length(find(j<0))*100/total_cell_lines;
end

[~,i] = sort(double(top_nodes_effect(:,3)),'descend');
top_nodes_effect = top_nodes_effect(i,:);

